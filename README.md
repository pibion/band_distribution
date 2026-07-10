This code is useful for dark matter searches where detector output is Ep (total phonon energy) and Eq (total charge energy).  The point of the code is to provide the probability of an (Ep, Eq) pair given a set of detector parameters for both electron recoils (`PpqG`, where the G is for gamma because gammas are the cause of most electron recoils) and neutron recoils (`PpqN`, where the N is for neutron).


## Building and testing with the Fortran Package Manager (`fpm`)

This project uses the Fortran Package Manager (fpm).  You'll need to install that to build this project; please see https://fpm.fortran-lang.org/install/index.html#install for instructions on installing fpm on your system.  Currently (Nov 2025), building from source will install version 0.14 while installing the package via e.g. `conda` will install version 0.12.

The code parallelizes its integration loops with `do concurrent` using locality specifiers, including the Fortran 2023 `reduce` clause.  This is a hard requirement — there is no serial fallback — so you need a recent compiler.  The versions below have been verified via the docker containers in this repository:

|Vendor| Version(s)      |  Build/Test Command                                                                                                                        |
|------|-----------------|--------------------------------------------------------------------------------------------------------------------------------------------|
|GNU   | 15.2.0          | `fpm test --compiler gfortran --profile release --flag "-march=native -fopenmp -ftree-parallelize-loops=4 -fcoarray=single -fPIC"`         |
|Intel | 2025.2.1        | `fpm test --compiler ifx --flag "-fpp -O3 -qopenmp -DHAVE_MULTI_IMAGE_SUPPORT=0" --profile release`                                        |
|LLVM  | 22              | `fpm test --compiler flang --profile release --flag "-O3 -fopenmp -fdo-concurrent-to-openmp=host"`                                         |

Notes:
* Intel: the shared library must also be linked against the Intel OpenMP runtime for the Python ctypes interface to work; the Intel container does this by setting `FPM_LDFLAGS="-liomp5"`.  The `-fpp -DHAVE_MULTI_IMAGE_SUPPORT=0` flags are for the Julienne dependency.
* LLVM: the compiler is invoked as `flang` (the `flang-new` name was dropped in LLVM 20).  `-fdo-concurrent-to-openmp=host` is what parallelizes the `do concurrent` loops; without it they compile to serial loops.  The LLVM container sets `FPM_LDFLAGS="-fopenmp"` so the shared library links against `libomp`.

# Testing the python calls
This code builds a library that may be called within python (this is the original intent of the code).  The python test scripts live in `test/python/` and should be run from the repository root.  To test the python calls, run

```
LD_LIBRARY_PATH=lib python test/python/test_PpqFort.py
```

or if you just want to test the vectorized functions `PpqN_vector` and `PpqG_vector`

```
LD_LIBRARY_PATH=lib python test/python/test_PpqFort_vectorFuncs.py
```

(`LD_LIBRARY_PATH=lib` is needed when running outside the docker containers, which set it in their environment.)

## Chi-square validation tests
These verify that data sampled from the PDFs produces Pearson chi-square values that follow the theoretical chi2(n_bins − 1) distribution.  Reference plots are committed in `figures/`.

```
python test/python/verify_sample_from_pdf.py                       # sampler moment checks (~1 min)
python test/python/test_chisquare_gaussian.py                      # analytic Gaussian PDF (~2 min)
LD_LIBRARY_PATH=lib python test/python/test_chisquare_ppqn.py      # Fortran PpqN PDF (~1 hr first run)
```

`test_chisquare_ppqn.py` caches its PDF grid evaluation in `ppqn_vertex_grid.npz` (not committed); the first run takes ~20 minutes to build it, subsequent runs reuse it.  Pass an integer argument to reduce the number of throws for a quick smoke test, e.g. `... test_chisquare_ppqn.py 2000`.

### Physics-simulator validation tests
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
LD_LIBRARY_PATH=lib python test/python/test_chisquare_nr_simulator.py 10000    # ~7 min
LD_LIBRARY_PATH=lib python test/python/test_chisquare_er_simulator.py 10000    # ~11 min
```

(Timings measured with 18 workers under x86 emulation on an Apple Silicon Mac; native x86 hardware should be faster.)

To run inside the Intel docker container, mount the repository's `figures/` directory so the plot survives the container:

```
docker run --rm -v $(pwd)/figures:/app/figures band_distribution_intel \
    bash --login -c "conda activate band && python /app/test/python/test_chisquare_nr_simulator.py 100 64"
```

# Build the singularity/apptainer container for HPC submissions
There are multiple Dockerfiles, each building the code with a compiler from a different vendor (GNU, Intel, and LLVM).  

|Vendor| Dockerfile name     |
|------|---------------------|
|GNU   | Dockerfile_gfortran | 
|Intel | Dockerfile_intel    | 
|LLVM  | Dockerfile_llvm     |

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

# Use the docker container for local development
For local development, you most likely want the files available to you in a way that persists once you close the container.  In this case you need to supply arguments to `docker run` that mount the top-level directory:

```
docker run -it --mount type=bind,src=.,dst=/app --entrypoint=/bin/bash band
```

# Build the docker container for running Jupyter and interacting with notebooks

```
docker build --rm -f Dockerfile_jupyter -t band_jupyter .
```

# Run the jupyter container
You can issue this command from any directory.  Note the absolute path names for mounting the volume.  This enables your work to persist!  You will need to replace `/mnt/c/Users/canto/Repositories/nrFanoII` with the path to your repository directory.  You should leave `home/jovyan/work/nrFano` the same.  Note that this command refers to the nrFanoII repository, which uses this (band_distribution) repository.

```
docker run -it --rm -p 8888:8888 -v /mnt/c/Users/canto/Repositories/nrFanoII:/home/jovyan/work/nrFano band_jupyter:latest
```

In notebooks, select the **Python (band)** kernel — it runs in the `band` conda environment (built from `environment.yaml`, the same environment used by the other containers), which has the compiled library's runtime dependencies and all the python packages.

# Profiling with TAU
TAU is built into `Dockerfile_tau_intel`. Build that image first:

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
