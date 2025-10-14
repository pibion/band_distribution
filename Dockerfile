# Use Ubuntu 22.04 as base image
FROM ubuntu:22.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Install necessary packages
RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    ca-certificates \
    curl \
    git \
    nano \
    gfortran \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set up working directory
WORKDIR /app

# Install Anaconda
RUN wget https://repo.anaconda.com/archive/Anaconda3-2024.10-1-Linux-x86_64.sh -O /tmp/anaconda.sh \
    && bash /tmp/anaconda.sh -b -p /opt/anaconda \
    && rm /tmp/anaconda.sh

# Add Anaconda to PATH
ENV PATH="/opt/anaconda/bin:${PATH}"

# Initialize conda for bash shell
RUN echo ". /opt/anaconda/etc/profile.d/conda.sh" >> ~/.bashrc

# Copy environment.yaml file
COPY NR_Fano_dan_env.yaml .

# Create conda environment from the yaml file
RUN conda env create -f NR_Fano_dan_env.yaml
RUN echo "conda activate NR_Fano" >> ~/.bashrc

# Copy the fortran code over
COPY PpqFort.f90 .

# Now compile the fortran code
RUN gfortran -fPIC -shared -O3 -march=native -ffast-math -fopenmp -ftree-parallelize-loops=4 -o PpqFort.so PpqFort.f90

# Set the shell to bash with login to ensure .bashrc is sourced
SHELL ["/bin/bash", "--login", "-c"]
