Reinforcement learning algorithm for adaptive phase estimation
==============================================================

This project aims to show the use of reinforcement learning algorithm in finding a feedback policy for adaptive quantum-enhanced measurement, which has the goal of achieving the scaling in precision that exceeds the conventional techniques. We use adaptive phase estimation that includes phase noise and loss as the test problem.
        
The program is designed to streamline the implementation of population-based optimization algorithms on high-performance computing clusters and support parallel computing using MPI.

Features
--------

  - Easy selection of optimization algorithm and optimization problem

  - Libraries for differential evolution and particle swarm optimization

  - Both uniformly random and clustered initialization of population 

  - Accept-reject criteria to ensure quantum-enhanced results

  - Support for VSL and GPU to provide fast random number generation

Usage
=====

The program is designed to work on HPC clusters and it requires MPI to run. The basic use is as follows:

    $ [mpirun -np NPROC] phase_estimation [config_file]

Arguments:

    config_file              Configuration file name

If it is run without a configuration file, some default values are taken for all parameters; the exact settings are identical to the one in the provided `default.cfg` file. The configuration file is a plain text file with the name of the parameter on the left, followed by an equation sign surrounded by a space on either side, and a value on the right-hand side. For example, the contents of `default.cfg` are as follows:

    pop_size = 20
    N_begin = 4
    N_cut = 5
    N_end = 10
    iter = 100
    iter_begin = 300
    repeat = 10
    output_filename = output.dat
    time_filename = time.dat

If you supply a configuration file, but do not set a specific value to every possible option, the default values are again the ones described in `default.cfg`.

The meaning of the individual parameters:

  - `pop_size`: population size.
  
  - `N_begin`: the starting number of particles.
  
  - `N_cut`: the number of particles where the program use cluster initialization around previous solution.
  
  - `N_end`: the final number of particles.

  - `iter`: number of iterations.
  
  - `iter_begin`:
  
  - `repeat`:
  
  - `output_filename`: the name of the file to write the results to.
  
  - `time_filename`: the name of the file to write the time taken to run the program for each number of variables.

Compilation
-----------

The project contains the support for compilation using Autotools, and has been tested using GNU Compile Chain (GCC) and Intel Compilers. The Intel VSL library and CUDA are automatically detected. An MPI implementation is required to compile and run the code.

If you cloned the Git repository, first run `autogen.sh` in order to create missing files and generate the executable configure from configure.ac. 

Follow the standard POSIX procedure:

    $ ./configure [options]
    $ make
    $ make install

To use the Intel compilers, set the following environment variables:

    export CC=/path/of/intel/compiler/icc
    export CXX=/path/of/intel/compiler/icpc
    export OMPI_CC=/path/of/intel/compiler/icc
    export OMPI_CXX=/path/of/intel/compiler/icpc

In order to use icc and icpc compilers, you have to set these variables
so the mpic++ will invoke icpc instead of the default compiler.

Options for configure

    --prefix=PATH           Set directory prefix for installation
    --with-mpi=MPIROOT      Use MPI root directory.
    --with-mpi-compilers=DIR or --with-mpi-compilers=yes
                              use MPI compiler (mpicxx) found in directory DIR, or
                              in your PATH if =yes
    --with-mpi-libs="LIBS"  MPI libraries [default "-lmpi"]
    --with-mpi-incdir=DIR   MPI include directory [default MPIROOT/include]
    --with-mpi-libdir=DIR   MPI library directory [default MPIROOT/lib]

The above flags allow the identification of the correct MPI library the user wishes to use. The flags are especially useful if MPI is installed in a non-standard location, or when multiple MPI libraries are available.

    --with-cuda=/path/to/cuda           Set path for CUDA

The configure script looks for CUDA in /usr/local/cuda. If your installation is not there, then specify the path with this parameter. If you do not want CUDA enabled, set the parameter to ```--without-cuda```.


  --with-vsl=PATH    prefix where Intel MKL/VSL is installed

Specify the path to the VSL installation with this parameter.

Acknowledgement
---------------
The computational work was enabled in part by support provided by [WestGrid](www.westgrid.ca) and [Calcul Quebec](www.calculquebec.ca) through [Compute Canada Calcul Canada](www.computecanada.ca).
