# Download and installation

## Download

* Latest release: [download](https://github.com/PanPalitta/phase_estimation/releases)
* Development version: [Github](https://github.com/PanPalitta/phase_estimation)
* Feedback and issues: [issue](https://github.com/PanPalitta/phase_estimation/issues)

## Compilation and installation

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
