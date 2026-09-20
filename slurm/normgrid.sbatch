#!/bin/bash
# Slurm job array for building a normalization table (python/normgrid.py).
# One array task = one chunk of grid points, single core, no communication;
# the whole grid is embarrassingly parallel.  A TEMPLATE: adjust the #SBATCH
# lines (partition/account/time) for your cluster.  Locally tested with a
# stand-in evaluator and the container paths below checked against the
# Dockerfiles, but not run on a live cluster.
#
# Usage (from the repository root, after `mkdir -p logs results`):
#
#   python python/normgrid.py make-spec --band NR --out spec_NR.json     # RECOMMENDED_NODES
#   python python/normgrid.py make-random-spec --band NR --n 200 --out held_NR.json
#   N=$(python python/normgrid.py chunks spec_NR.json --size 100 | wc -l)
#   sbatch --array=0-$((N-1))%200 slurm/normgrid.sbatch spec_NR.json 100 results
#   # same for the held-out validation points (used by `validate` below):
#   N=$(python python/normgrid.py chunks held_NR.json --size 20 | wc -l)
#   sbatch --array=0-$((N-1)) slurm/normgrid.sbatch held_NR.json 20 results
#
#   # when the arrays are done (quick; run on a login node or interactively):
#   python python/normgrid.py merge --spec spec_NR.json --results 'results/res_spec_NR_*.txt' \
#          --out tables/norm_NR_region2-200_4-100.h5
#   python python/normgrid.py validate --table tables/norm_NR_region2-200_4-100.h5 \
#          --heldout-spec held_NR.json --results 'results/res_held_NR_*.txt'
#
# Arguments: SPEC.json  [CHUNK_SIZE=100]  [RESULTS_DIR=results]
#   Pick CHUNK_SIZE so a task runs ~10-30 min (a point takes ~3-10 s of one
#   core), and the %200 throttle to what the cluster is happy with.
#
# Robustness: `normgrid run` is resumable and crash-tolerant (a point whose
# integral errors out is recorded as nan and the worker restarts past it).
# If tasks are killed (time limit, preemption) or fail, just submit the same
# command again -- finished points are skipped, so only unfinished chunks do
# real work.  Failed points are re-tried with a looser tolerance via
#   python python/normgrid.py run --spec SPEC --start A --stop B --out FILE \
#          --retry-failed --epsrel 1e-6
#
# Environment (all optional):
#   BAND_SIF         container image (default $HOME/band.sif, built from
#                    Dockerfile_intel or _llvm as in the README).  Set
#                    BAND_NATIVE=1 to skip the container and use the repo's own
#                    python/ and lib/ (i.e. you built the library on the cluster).
#   BAND_PYTHON      python inside the container (default: the image's conda
#                    'band' env; the login-shell conda activation in the
#                    Dockerfile does not happen under apptainer)

#SBATCH --job-name=normgrid
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --requeue
#SBATCH --output=logs/normgrid_%A_%a.out
#SBATCH --error=logs/normgrid_%A_%a.err

set -euo pipefail

: "${SLURM_ARRAY_TASK_ID:?submit as a job array: sbatch --array=0-N-1 slurm/normgrid.sbatch SPEC [SIZE] [DIR] (see the header)}"
SPEC=${1:?usage: slurm/normgrid.sbatch SPEC.json [CHUNK_SIZE] [RESULTS_DIR]}
SIZE=${2:-100}
OUT_DIR=${3:-results}
mkdir -p "$OUT_DIR"

# One thread per task: with an ifx/flang build every worker would otherwise
# start its own OpenMP pool and the tasks would oversubscribe their core.
export OMP_NUM_THREADS=1

if [[ "${BAND_NATIVE:-0}" == "1" ]]; then
    RUN=()          # (expanded as ${RUN[@]+"${RUN[@]}"}: an empty array is an error under set -u on bash < 4.4)
    PY=${BAND_PYTHON:-python}
    NORMGRID=python/normgrid.py
    export LD_LIBRARY_PATH="$PWD/lib:${LD_LIBRARY_PATH:-}"
else
    RUN=(apptainer exec --bind "$PWD:$PWD" --pwd "$PWD" "${BAND_SIF:-$HOME/band.sif}")
    PY=${BAND_PYTHON:-/opt/anaconda/envs/band/bin/python}
    NORMGRID=/app/python/normgrid.py
    export APPTAINERENV_LD_LIBRARY_PATH="/app/lib:/opt/anaconda/envs/band/lib"
    export SINGULARITYENV_LD_LIBRARY_PATH=$APPTAINERENV_LD_LIBRARY_PATH   # older singularity
fi

# this task's chunk: line (task id + 1) of the chunk list (a bad spec makes
# the next line fail loudly under `set -e`)
CHUNKS=$(${RUN[@]+"${RUN[@]}"} "$PY" "$NORMGRID" chunks "$SPEC" --size "$SIZE")
read -r START STOP <<< "$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" <<< "$CHUNKS")"
if [[ -z "${START:-}" ]]; then
    echo "array task $SLURM_ARRAY_TASK_ID is past the last chunk of $SPEC; nothing to do"
    exit 0
fi

STEM=$(basename "$SPEC" .json)
OUT="$OUT_DIR/res_${STEM}_${START}.txt"
echo "[$(date +%T)] $STEM points [$START, $STOP) -> $OUT on $(hostname)"

# nonzero exit only means some points failed (recorded as nan in $OUT);
# report it but let the task finish so the result file is kept
${RUN[@]+"${RUN[@]}"} "$PY" "$NORMGRID" run --spec "$SPEC" --start "$START" --stop "$STOP" --out "$OUT" || \
    echo "WARNING: some points in [$START, $STOP) failed; see $OUT (retry with --retry-failed --epsrel 1e-6)"
