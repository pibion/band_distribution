This code is useful for dark matter searches where detector output is Ep (total phonon energy) and Eq (total charge energy).  The point of the code is to provide the probability of an (Ep, Eq) pair given a set of detector parameters for both electron recoils (`PpqG`, where the G is for gamma because gammas are the cause of most electron recoils) and neutron recoils (`PpqN`, where the N is for neutron).


# Compile with
You can replace the 4 with the number of processes you'd like to parallelize to.

```
gfortran -fPIC -shared -O3 -march=native -ffast-math -fopenmp -ftree-parallelize-loops=4 -o PpqFort.so PpqFort.f90
```

# And then test with
```
python test_PpqFort.py
```

or if you just want to test the vectorized functions `PpqN_vector` and `PpqG_vector`

```
python test_PpqFort_vectorFuncs.py
```

# Build the singularity/apptainer container for HPC submissions
You should issue the following command in the `fortran-python` directory:

```
docker build -f Dockerfile -t fano_fort .
```

Now you have a docker container, but this is not usable on HPC systems.  Run this command to create `fano_fort.sif`, an image file that can be used on HPC systems.  The command can be issued in any location (it is not directory dependent).

```
apptainer build fano_fort.sif docker-daemon://fano_fort:latest

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
I ran into trouble installing TAU, so I recommend using Docker.

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

It's probably wise to make the code compiles normally:

```
gfortran -O2 PpqFort.f90 test_ppq.f90 -o test_ppq
```

Before trying to instrument with TAU:

```
rm test_ppq
tau_exec gfortran -O2 PpqFort.f90 test_ppq.f90 -o test_ppq
```
