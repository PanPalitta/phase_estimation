##Reinforcement learning algorithm for adaptive phase estimation

This project aims to show the use of reinforcement learning algorithm in finding a feedback policy for adaptive quantum-enhanced measurement, which has the goal of achieving the scaling in precision that exceeds the conventional techniques. We use adaptive phase estimation that includes phase noise and loss as the test problem.
        
The program is designed to streamline the implementation of population-based optimization algorithms on high-performance computing clusters and support parallel computing using MPI.

**Features:** 

- Easy selection of optimization algorithm and optimization problem

- Libraries for differential evolution and particle swarm optimization

- Both uniformly random and clustered initialization of population 

- Accept-reject criteria to ensure quantum-enhanced results

- Support for VSL and GPU to provide fast random number generation


**Limitation**

-The adaptive phase estimation is only reliable to 100 photons [Other limitation?] 

### Usage:

As the program is designed to work on HPC clusters where interactive input from user might not be allowed, all the parameters are set within the code before compilation.

#### Selecting problem and optimization algorithm 

The library for the optimization algorithm and the problem must be included in the headers of main.cpp. The objects in the corresponding classes can be set in the main function as pointers.

```
		#include "mpi_de.h"
		#include "mpi_pso.h"
``` 

The optimization algorithms are instantiated with default parameters which can be changed within the code. For example, for differential evolution:

```
        ...

        DE(Problem<typeT> *problem_ptr): F(0.1), Cr(0.6) {

        ...
```

The library also contain a function for changing the parameters during runtime:

```
		void DE<typeT>::write_param(double *param_array) {
		...

```

#### Parameters setting 

The beginning of the main function contains a set of parameters including

-the smallest (N\_begin) and the largest number of variables (N\_end). 

-number of variables where the program use cluster initialization around previous solution (N\_cut)

-number of variables where accept-reject criteria starts (data\_end)

-population size (pop\_size)

-number of iterations (iter, iter\_begin)

#### Compilation

The project contains the support for compilation using Autotools, and has been tested using GNU Compile Chain (GCC) and Intel Compilers. Intel VSL library and GPU are automatically detected.

In the first compilation, first run autogen.sh in order to create missing files and generate the executable configure from configure.ac. Then run the executable to generate Makefile from Makefile.in. The code can now be compiled. The sequence of the commands are

```
        ./autogen.sh

        ./configure

        make
``` 

Any job submission command should be directed to run src/phase\_estimation. The solution from the run are in the output file output.dat, and the time taken to run the program for each number of variables is in the file time.dat. 


###Acknowledgement:

The computational work was enabled in part by support provided by WestGrid (www.westgrid.ca) and Calcul Quebec (www.calculquebec.ca) through Compute Canada Calcul Canada (www.computecanada.ca).

[Should we include the funding sources: NSERC, AITF]

###References: 

[To be include: manuscript on arxiv and manual -- what should be in the manual?]