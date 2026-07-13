This code is useful for dark matter searches where detector output is Ep (total phonon energy) and Eq (total charge energy).  The point of the code is to provide the probability of an (Ep, Eq) pair given a set of detector parameters for both electron recoils (`PpqG`, where the G is for gamma because gammas are the cause of most electron recoils) and neutron recoils (`PpqN`, where the N is for neutron).


# Building and testing with the Fortran Package Manager (`fpm`)

This project uses the Fortran Package Manager (fpm).  You'll need to install that to build this project; please see https://fpm.fortran-lang.org/install/index.html#install for instructions on installing fpm on your system.  Currently (Nov 2025), building from source will install version 0.14 while installing the package via e.g. `conda` will install version 0.12.

The code parallelizes its integration loops with `do concurrent` using locality specifiers, including the Fortran 2023 `reduce` clause, so you need a recent compiler.  The versions below have been verified via the docker containers in this repository:

|Vendor| Version(s)      |  Build/Test Command                                                                                                                        |
|------|-----------------|--------------------------------------------------------------------------------------------------------------------------------------------|
|GNU   | 15.2.0          | `fpm test --compiler gfortran --profile release --flag "-march=native -fopenmp -ftree-parallelize-loops=4 -fcoarray=single -fPIC"`         |
|Intel | 2025.2.1        | `fpm test --compiler ifx --flag "-fpp -O3 -qopenmp -DHAVE_MULTI_IMAGE_SUPPORT=0" --profile release`                                        |
|LLVM  | 22              | `fpm test --compiler flang --profile release --flag "-O3 -fopenmp -fdo-concurrent-to-openmp=host"`                                         |

Notes:
* GNU: gfortran is slow.  It has no option for mapping `do concurrent` onto threads (the auto-parallelizer in `-ftree-parallelize-loops` cannot parallelize these loops because they contain function calls), so the band integrals run **serially** in gfortran builds — measured at 1.0 effective threads vs 17.6 for ifx on the same machine, making per-event likelihood evaluation ~25x slower there (roughly the core count times a ~1.5x per-core gap from the vector math library).  The gfortran build is fine for verifying the physics and the Python interface, but use ifx or flang for production likelihood work (MCMC).
* Intel: the shared library must also be linked against the Intel OpenMP runtime for the Python ctypes interface to work; the Intel container does this by setting `FPM_LDFLAGS="-liomp5"`.  The `-fpp -DHAVE_MULTI_IMAGE_SUPPORT=0` flags are for the Julienne dependency.  `-qopenmp` maps `do concurrent` onto the OpenMP thread pool.
* LLVM: the compiler is invoked as `flang`.  `-fdo-concurrent-to-openmp=host` is what parallelizes the `do concurrent` loops; without it they compile to serial loops.  The LLVM container sets `FPM_LDFLAGS="-fopenmp"` so the shared library links against `libomp`.  Intel and LLVM builds benchmark identically (~6 us per event in vector mode on 18 emulated cores).

## Building the shared library for the python interface
The python wrappers load `lib/libband_distribution.so`, which `fpm install` builds and places under the repository root.  Use the same flags as the test commands above; for Intel and LLVM, `FPM_LDFLAGS` must also be set so the *shared library* links its OpenMP runtime — without it the library builds but fails to load from python with undefined `__kmpc_*` symbols:

```
# Intel
FPM_LDFLAGS="-liomp5" fpm install --prefix=. --compiler ifx --profile release --flag "-fpp -O3 -qopenmp -DHAVE_MULTI_IMAGE_SUPPORT=0"

# LLVM
FPM_LDFLAGS="-fopenmp" fpm install --prefix=. --compiler flang --profile release --flag "-O3 -fopenmp -fdo-concurrent-to-openmp=host"

# GNU (verification only — see the notes above)
fpm install --prefix=. --compiler gfortran --profile release --flag "-march=native -fopenmp -ftree-parallelize-loops=4 -fcoarray=single -fPIC"
```

# Testing the python calls
This code builds a library that may be called within python (this is the original intent of the code).  Build and install the shared library first (see "Building the shared library for the python interface" above).  The python test scripts live in `test/python/` and should be run from the repository root.  To test the python calls, run

```
LD_LIBRARY_PATH=lib python test/python/test_PpqFort.py
```

or if you just want to test the vectorized functions `PpqN_vector` and `PpqG_vector`

```
LD_LIBRARY_PATH=lib python test/python/test_PpqFort_vectorFuncs.py
```

(`LD_LIBRARY_PATH=lib` is needed when running outside the docker containers, which set it in their environment.)

# PDF self-consistency tests
These verify that data sampled from the PDFs produces Pearson chi-square values that follow the theoretical chi2(n_bins − 1) distribution.  Reference plots are committed in `figures/`.

```
python test/python/verify_sample_from_pdf.py                       # sampler moment checks (~2 s)
python test/python/test_chisquare_gaussian.py                      # analytic Gaussian PDF (~2 min)
LD_LIBRARY_PATH=lib python test/python/test_chisquare_ppqn.py      # Fortran PpqN PDF (~2 min first run)
```

`test_chisquare_ppqn.py` takes the same two optional positional arguments as the simulator tests below: `n_throws` (default 100,000) and `n_bins` (default 64).  For a quick smoke test with fewer throws, run

```
LD_LIBRARY_PATH=lib python test/python/test_chisquare_ppqn.py 2000 64
```

It caches its PDF grid evaluation in `ppqn_vertex_grid.npz` (not committed) and reuses it on later runs.

The `test_chisquare_ppqn.py` timing assumes an ifx or flang build; gfortran builds run the band integrals serially (see the compiler notes above) and take correspondingly longer.

# Physics-simulator validation tests
The tests above only check that the PDFs are *self-consistent* (samples drawn from a PDF match that same PDF).  The two tests below are the stronger check: they generate events from an independent physics simulator (`python/generate_events.py`, which draws Er from the recoil spectrum, N from a truncated normal, and applies detector resolution) and compare the binned counts against the Fortran PDFs.  A pass means the Fortran `PpqN` / `PpqG` implementations correctly describe the physics.

```
LD_LIBRARY_PATH=lib python test/python/test_chisquare_nr_simulator.py [n_throws] [n_bins]   # NR band vs PpqN
LD_LIBRARY_PATH=lib python test/python/test_chisquare_er_simulator.py [n_throws] [n_bins]   # ER band vs PpqG
```

Both write their chi-square histograms to `figures/chisquare_nr_simulator.png` and `figures/chisquare_er_simulator.png`.  `n_throws` defaults to 1000 and `n_bins` to 400.  The one-time integration of the PDF over the bins dominates the wall time, so reducing `n_bins` is the way to get a fast smoke test:

```
# smoke test: 100 throws, 64 bins
LD_LIBRARY_PATH=lib python test/python/test_chisquare_nr_simulator.py 100 64   # ~1 min
LD_LIBRARY_PATH=lib python test/python/test_chisquare_er_simulator.py 100 64   # ~1 min

# full validation: 10,000 throws, 400 bins
LD_LIBRARY_PATH=lib python test/python/test_chisquare_nr_simulator.py 10000    # ~4 min
LD_LIBRARY_PATH=lib python test/python/test_chisquare_er_simulator.py 10000    # ~5 min
```

(Timings measured with the ifx build, 18 workers, under x86 emulation on an Apple Silicon Mac; native x86 hardware should be faster.  gfortran builds run the band integrals serially and will be dramatically slower here — use ifx or flang.)

To run inside the Intel docker container, mount the repository's `figures/` directory so the plot survives the container:

```
docker run --rm -v $(pwd)/figures:/app/figures band_distribution_intel \
    bash --login -c "conda activate band && python /app/test/python/test_chisquare_nr_simulator.py 100 64"
```

# Performance notes
`PpqN` / `PpqG` integrate over a window placed around the located peak(s) of the Er integrand rather than sampling the full physical range, evaluating at roughly 6 microseconds per (Ep, Eq) point in vector mode (18 threads under x86 emulation; measured via `PpqN_vector` over a representative grid).  This number assumes an ifx or flang build — gfortran builds run the band integrals serially and are ~25x slower (see the compiler notes at the top).  For likelihood loops (e.g. MCMC):

* Call the vectorized entry points (`PpqN_vector` / `PpqG_vector`) with all events in one call — the parallelism lives there, and per-event scalar calls pay OpenMP fork/join overhead instead.
* **Shuffle the event array once at load time if it is ordered.**  The vector loops split the events into one contiguous chunk per thread (static scheduling), and the loop only finishes when the slowest chunk does.  Per-event cost varies several-fold across the (Ep, Eq) plane — deep-tail events short-circuit in ~2 us while on-band events cost ~10-15 us — so an energy-ordered array hands some threads chunks of expensive events while others idle at the barrier.  Shuffling gives every chunk a similar cost mix and measured 8-15% faster than energy-ordered input.  The result is identical either way, and if your events are already in effectively random order this changes nothing.

# Build the singularity/apptainer container for HPC submissions
There are multiple Dockerfiles, each building the code with a compiler from a different vendor (GNU, Intel, and LLVM).  

|Vendor| Dockerfile name     | Notes |
|------|---------------------|-------|
|GNU   | Dockerfile_gfortran | Slow: gfortran runs the band integrals serially, ~25x slower than ifx/flang (see the compiler notes above) |
|Intel | Dockerfile_intel    | Recommended for work that needs speed |
|LLVM  | Dockerfile_llvm     | Same performance as the Intel build |

Choose which compiler you want, determine the name of the dockerfile, and then issue the following command:

```
docker build -f {dockerfile name} -t band .
```

If you need to troubleshoot the docker build, you can shell into this container with the command

```
docker run -it band
```

Now you have a docker container that contains the fortran binary, but this is not usable on HPC systems.  Run this command to create `band.sif`, an image file that can be used on HPC systems.  The command can be issued in any location (it is not directory dependent).

```
apptainer build band.sif docker-daemon://band:latest
```

# Build the docker container for running Jupyter and interacting with notebooks

```
docker build --rm -f Dockerfile_jupyter -t band_jupyter .
```

## Run the jupyter container
You can issue this command from any directory.  Note the absolute path names for mounting the volume.  This enables your work to persist!  You will need to replace `/mnt/c/Users/canto/Repositories/nrFanoII` with the path to your repository directory.  You should leave `home/jovyan/work/nrFano` the same.  Note that this command refers to the nrFanoII repository, which uses this (band_distribution) repository.

```
docker run -it --rm -p 8888:8888 -v /mnt/c/Users/canto/Repositories/nrFanoII:/home/jovyan/work/nrFano band_jupyter:latest
```

In notebooks, select the **Python (band)** kernel — it runs in the `band` conda environment (built from `environment.yaml`, the same environment used by the other containers), which has the compiled library's runtime dependencies and all the python packages.

# Use the docker container for local development
For local development, you most likely want the files available to you in a way that persists once you close the container.  In this case you need to supply arguments to `docker run` that mount the top-level directory:

```
docker run -it --mount type=bind,src=.,dst=/app --entrypoint=/bin/bash band
```

# Profiling with TAU
TAU is built into `Dockerfile_tau_intel`. This image uses the ifx compiler because that comipler is the easiest to install outside a Docker container.  Build that image first:

```
docker build -f Dockerfile_tau_intel -t tau_intel .
```

Once the image is built, you can run the container and shell into it with. The repository is mounted at `/repo` so that any changes to files persist on the host. The built library remains at `/app`.

```
# for fish shell
docker run -it -v /tmp/.X11-unix:/tmp/.X11-unix -v "$XAUTH/.Xauthority:/root/.Xauthority" -v $PWD:/repo tau_intel /bin/bash

# for bash shell
docker run -it -v /tmp/.X11-unix:/tmp/.X11-unix -v "$XAUTH/.Xauthority:/root/.Xauthority" -v ($pwd):/repo tau_intel /bin/bash
```

Inside the docker prompt, you'll need to set the DISPLAY environment variable.  After this, tools that require a GUI window like `paraprof` will work.

```
export DISPLAY=:0.0
```

Set the TAU environment variables:

```
export TAU_MAKEFILE=/packages/tau2/x86_64/lib/Makefile.tau-pthread
export TAU_OPTIONS="-optCompInst -optVerbose -optNoRevert"
export TAU_THROTTLE=0
```

Compile the profiling driver using `tau_f90.sh`, which instruments every Fortran source file.
The module and submodule must be compiled before the driver (the `-c` step generates the `.mod` files needed downstream):

```
cd /tmp

tau_f90.sh -O3 -g -qopenmp -fpp -I/app/include -c /app/src/PpqFort_m.f90
tau_f90.sh -O3 -g -qopenmp -fpp -I/app/include -c /app/src/PpqFort_s.f90

tau_f90.sh -O3 -g -qopenmp -fpp \
    -I/app/include \
    PpqFort_m.o PpqFort_s.o /repo/app/profile_driver.f90 \
    -L/app/lib -lband_distribution \
    -Wl,-rpath,/app/lib \
    -o /repo/profile_driver
```

Run the instrumented binary directly (no `tau_exec` needed with compiler instrumentation):

```
/repo/profile_driver
```

Running the code will produce `profile.*` files in the current directory (`/tmp` if you followed the steps above), which you can investigate using `pprof` (terminal summary) or `paraprof` (GUI view).

## Re-profiling after a source change

If you edit source files in `/repo` and want to re-profile, rebuild the library from `/repo` first, then recompile the profile driver against the new library.

**Step 1: Rebuild the library**
```
cd /repo
fpm install --compiler ifx --flag "-fpp -O3 -qopenmp -DHAVE_MULTI_IMAGE_SUPPORT=0" --profile release --prefix /tmp/band_new
```

**Step 2: Add the dependency libraries to `LD_LIBRARY_PATH`**

`fpm install` copies the main library but not its dependencies (`libassert.so`, `libjulienne.so`). Point the linker at the fpm build directory:
```
export LD_LIBRARY_PATH=$(find /repo/build -name "libassert.so" -printf "%h\n" | head -1):$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$(find /repo/build -name "libjulienne.so" -printf "%h\n" | head -1):$LD_LIBRARY_PATH
```

**Step 3: Recompile the profile driver with TAU instrumentation**
```
cd /tmp

tau_f90.sh -O3 -g -qopenmp -fpp -I/tmp/band_new/include -c /repo/src/PpqFort_m.f90
tau_f90.sh -O3 -g -qopenmp -fpp -I/tmp/band_new/include -c /repo/src/PpqFort_s.f90

tau_f90.sh -O3 -g -qopenmp -fpp \
    -I/tmp/band_new/include \
    PpqFort_m.o PpqFort_s.o /repo/app/profile_driver.f90 \
    -L/tmp/band_new/lib -lband_distribution \
    -Wl,-rpath,/tmp/band_new/lib \
    -o /repo/profile_driver
```

**Step 4: Run and inspect**
```
/repo/profile_driver
pprof
```

# Documentation
With [ford](https://github.com/Fortran-FOSS-Programmers/ford) installed, run `ford ford.md`.
Then open `doc/html/index.html` in a web browser to see the band_distribution documentation.
