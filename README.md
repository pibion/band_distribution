This code is useful for dark matter searches where detector output is Ep (total phonon energy) and Eq (total charge energy).  The point of the code is to provide the probability of an (Ep, Eq) pair given a set of detector parameters for both electron recoils (`PpqG`, where the G is for gamma because gammas are the cause of most electron recoils) and neutron recoils (`PpqN`, where the N is for neutron).


## Building and testing with the Fortran Package Manager (`fpm`)

This project uses the Fortran Package Manager (fpm).  You'll need to install that to build this project; please see https://fpm.fortran-lang.org/install/index.html#install for instructions on installing fpm on your system.  Currently (Nov 2025), building from source will install version 0.14 while installing the package via e.g. `conda` will install version 0.12.

The commands below should work with `fpm` 0.12.0 and with the compiler versions shown.
With `fpm` releases more recent than 0.12.0, one can replace `flang-new` with `flang`.

|Vendor| Version(s)      |  Build/Test Command                                          |
|------|-----------------|--------------------------------------------------------------|
|GNU   | 14.3.0, 15.2.0  | `fpm test --compiler gfortran --profile release --flag "-cpp -march=native -fopenmp -ftree-parallelize-loops=4"` |
|      | 13.4.0          | `fpm test --compiler gfortran --profile release --flag "-cpp -march=native -fopenmp -ftree-parallelize-loops=4 -ffree-line-length-none"` |
|Intel | 2025.2.1        | `FOR_COARRAY_NUM_IMAGES=1 fpm test --compiler ifx --flag "-fpp -O3 -coarray" --profile release` |
|LLVM  | 20-22           | `fpm test --compiler flang-new --profile release --flag "-cpp -O3"` |
|      | 19              | `fpm test --compiler flang-new --profile release --flag "-cpp -O3 -mmlir -allow-assumed-rank"` |
|NAG   | 7.2, Build 7235 | `fpm test --compiler nagfor --flag "-fpp -O4"`               |

**Caveat:** In the case of LLVM 19-20, the above commands succeed for testing band_distribution's Julienne dependency.
A future pull request could test band_distribution itself with LLVM 19-20 via GitHub Actions.

### Preprocessor macros
* Add `-DPREFER_DO_CONCURRENT` in the `--flag` argument to switch from array statements or `do` loops to `do concurrent`

### Experimental parallelization
To multithread `do concurrent` on CPUs with LLVM `flang` 21 or later, try the following:
```
fpm test --compiler flang-new --profile release --flag "-O3 -cpp -DPREFER_DO_CONCURRENT -fopenmp -fdo-concurrent-to-openmp=host"
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
|GNU   | Dockerfile          | 
|Intel | Dockerfile_intel    | 
|LLVM  | Dockerfile_llvm     |

Choose which compiler you want, determine the name of the dockerfile, and then issue the following command in the `fortran-python` directory:

```
docker build -f {dockerfile name} -t fano_fort .
```

If you need to troubleshoot the docker build, you can shell into this container with the command

```
docker run -it --entrypoint /bin/bash fano_fort
```

Now you have a docker container that contains the fortran binary, but this is not usable on HPC systems.  Run this command to create `fano_fort.sif`, an image file that can be used on HPC systems.  The command can be issued in any location (it is not directory dependent).

```
apptainer build fano_fort.sif docker-daemon://fano_fort:latest
```

# Use the docker container for local development
For local development, you most likely want the files available to you in a way that persists once you close the container.  In this case you need to supply arguments to `docker run` that mount the top-level directory:

```
docker run -it --mount type=bind,src=.,dst=/app --entrypoint=/bin/bash fano_fort
```

# Build the docker container for running Jupyter and interacting with notebooks
You should issue the following command in the `fortran-python` directory:

```
docker build --rm -f Dockerfile_jupyter -t fano_jupyter .
```

# Run the jupyter container
You can issue this command from any directory.  Note the absolute path names for mounting the volume.  This enables your work to persist!  You will need to replace `/mnt/c/Users/canto/Repositories/nrFanoII` with the path to your repository directory.  You should leave `home/jovyan/work/nrFano` the same.

```
docker run -it --rm -p 8888:8888 -v /mnt/c/Users/canto/Repositories/nrFanoII:/home/jovyan/work/nrFano fano_jupyter:latest
```

# Profiling with TAU
I ran into trouble installing TAU, so I recommend using Docker. I use the container created by E4S (https://e4s.io/):

```
docker pull ecpe4s/e4s-cpu:latest
```

Once the dockerfile is pulled, you can run the container and shell into it with

```
# for fish shell
docker run -it -v /tmp/.X11-unix:/tmp/.X11-unix -v "$XAUTH/.Xauthority:/root/.Xauthority" -v $PWD:/app ecpe4s/e4s-cpu:latest

# for bash shell
docker run -it -v /tmp/.X11-unix:/tmp/.X11-unix -v "$XAUTH/.Xauthority:/root/.Xauthority" -v ($pwd):/app ecpe4s/e4s-cpu:latest
```

Inside the docker prompt, you'll need to set the DISPLAY environment variable.  After this, tools that require a GUI window like `paraprof` will work.

```
export DISPLAY=:0.0
```

And you'll need to move into the mounted directory that contains the fortran files:

```
# this command assumes you're in the root directory
cd app
```

E4S currently uses an older (pre-13) version of gcc, so for the gfortran compiler to work we'll have to set `FCFLAGS`:

```
 FCFLAGS+=/usr/lib/gcc/x86_64-linux-gnu/11/
```

Now we can check that the code compiles normally:

```
gfortran -O2 PpqFort.f90 test_ppq.f90 -o test_ppq
```

And now we can try to instrument with TAU:

```
rm test_ppq
export TAU_MAKEFILE=/spack/opt/spack/linux-ubuntu22.04-x86_64_v3/gcc-11.4.0/tau-2.34.1-qruklusqgoww5pzgc4f2ffcpybmkbzpy/lib/Makefile.tau-papi-ompt-mpi-pthread-python-pdt-openmp
export TAU_OPTIONS="-optCompInst -optVerbose -optNoRevert"
export TAU_THROTTLE=0
tau_f90.sh -O2 PpqFort.f90 test_ppq.f90 -o test_ppq
```

And then run the instrumented code with 

```
mpirun -np 1 ./test_ppq
```

In other circumstances you'd just run `./test_ppq` but in this case that doesn't work because TAU was compiled with MPI support and this is a non-MPI application.

Running the code will produce a file `profile.0.0.0` in the `/app` directory which you can investigate using `pprof` (terminal summary) or `paraprof` (GUI view).

# Documentation
With [ford](https://github.com/Fortran-FOSS-Programmers/ford) installed, run `ford ford.md`.
Then open `doc/html/index.html` in a web browser to see the band_distribution documentation.
