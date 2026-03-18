This code is useful for dark matter searches where detector output is Ep (total phonon energy) and Eq (total charge energy).  The point of the code is to provide the probability of an (Ep, Eq) pair given a set of detector parameters for both electron recoils (`PpqG`, where the G is for gamma because gammas are the cause of most electron recoils) and neutron recoils (`PpqN`, where the N is for neutron).


## Building and testing with the Fortran Package Manager (`fpm`)

This project uses the Fortran Package Manager (fpm).  You'll need to install that to build this project; please see https://fpm.fortran-lang.org/install/index.html#install for instructions on installing fpm on your system.  Currently (Nov 2025), building from source will install version 0.14 while installing the package via e.g. `conda` will install version 0.12.

The commands below should work with `fpm` 0.12.0 and with the compiler versions shown.
With `fpm` releases more recent than 0.12.0, one can replace `flang-new` with `flang`.

|Vendor| Version(s)      |  Build/Test Command                                                                                                                      |
|------|-----------------|------------------------------------------------------------------------------------------------------------------------------------------|
|GNU   | 14.3.0, 15.2.0  | `fpm test --compiler gfortran --profile release --flag "-cpp -march=native -fopenmp -ftree-parallelize-loops=4"`                         |
|      | 13.4.0          | `fpm test --compiler gfortran --profile release --flag "-cpp -march=native -fopenmp -ftree-parallelize-loops=4 -ffree-line-length-none"` |
|Intel | 2025.2.1        | `fpm test --compiler ifx --flag "-fpp -O3 -qopenmp -DHAVE_MULTI_IMAGE_SUPPORT=0" --profile release`                                          |
|LLVM  | 20-22           | `fpm test --compiler flang-new --profile release --flag "-cpp -O3"`                                                                      | 
|      | 19              | `fpm test --compiler flang-new --profile release --flag "-cpp -O3 -mmlir -allow-assumed-rank"`                                           |
|NAG   | 7.2, Build 7235 | `fpm test --compiler nagfor --flag "-fpp -O4"`                                                                                           |

**Caveat:** In the case of LLVM 19-20, the above commands succeed for testing band_distribution's Julienne dependency.
A future pull request could test band_distribution itself with LLVM 19-20 via GitHub Actions.

### Preprocessor macros
* Add `-DCANNOT_DO_CONCURRENT=1` in the `--flag` argument to turn off `do concurrent` code if you have a compiler that does not implement this feature.

### Experimental parallelization
To multithread `do concurrent` on CPUs with LLVM `flang` 21 or later, try the following:
```
fpm test --compiler flang-new --profile release --flag "-O3 -cpp -fopenmp -fdo-concurrent-to-openmp=host"
```

# Testing the python calls
This code builds a library that may be called within python (this is the original intent of the code).  To test the python calls, run

```
python test_PpqFort.py
```

or if you just want to test the vectorized functions `PpqN_vector` and `PpqG_vector`

```
python test_PpqFort_vectorFuncs.py
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
docker run -it --entrypoint /bin/bash band
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
