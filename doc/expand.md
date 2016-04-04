# Expanding the library

The intention of this project is to create a library that can be used for solving multiple quantum control problems. The code is therefore designed to ease the process of including new problems and algorithms and make the selection of problems and algorithm as error-free as possible. The following document is a guide on how users can write and include their own problems and optimization algorithms to the existing library, and what is needed to customize and compile the code.

Readers are assumed to be familiar with population-based optimization algorithm, C++, object-oriented programming, class hierarchy and inheritance, and polymorphism.

## Program structure

![Diagram of program](..\Phase_Est_diagram2.png)

### User-specified components

The orange boxes correspond to the components in which the users specifies before the compiling the program. `Phase` class contains the modules for the adaptive phase estimation problem, which can be replaced with other problems. To select a problem of choice, replace `Phase()` by the constructor of the class in `main()` in the following line.
```
problem = new Phase(numvar, gaussian_rng, uniform_rng);
```

The Phase class is accessed through the `Problem` class. The pointer is given to the `OptAlg` class to be used for computing the fitness values and accept-reject criteria.

### Optimization algorithms

The choice of optimization algorithm is specified in the configuration file along with other parameters shown as an orange oval in figure. Otherwise, it can also be coded in `main()` in the following line if necessary.

```
opt = new DE(problem, gaussian_rng, pop_size);
```

### MPI

The MPI library is required for the program to run, as the program is designed to spread the solution candidate evenly on a group of processors. The processors communicates in the following situations.

* Constructing new set of candidates from existing population
* Finding the best candidate in the population or subset of the population
* Selecting the best candidate as solution

## Add a new problem

A new problem should be written as a class derived from `Problem` class. There are five functions that user must include in the new problem.

* `fitness()` is a function intended to be a wrapper for changing conditions in which the fitness function is evaluated.
* `avg_fitness()` is the function for calculating the fitness value.
* `T_condition()` is a function for calculating additional conditions for when the optimization algorithm is set to accept solution after T iteration.
* `error_condition()` is a function for calculating additional conditions for when optimization algorithm is set to accept solution from error bound.
* `boundary()` is used to keep the solution candidate within the boundary of the search space.

This class does not use any MPI functionalities.

## Add new algorithm

New algorithms can be added to the library of optimization algorithms by creating a derived class from the `OptAlg` class. The functions are designed based on swarm intelligence algorithms and evolutionary algorithms which shares the same backbone functions for initializing the population, for selecting the final solution candidate, and so forth. Aspects that are specific to the algorithm, such as how the new candidates are generated and selected, are declared as virtual function in `OptAlg` to allow the functions to be called from the derived class.

The functions, including virtual functions, are listed in OptAlg class document.

## Constructing optimization algorithm

The library provides the module for users to contract the optimization algorithm in `main()`. The basic structure of the algorithm is given in main.cpp, which compile to the following structure.

* Initialize MPI and setting is compute to spread the number of candidates evenly on the processors
* Initialize problem and select the optimization algorithm
* Population is initialized using a user specified method
* The fitness values are computed and the population is prepared for the optimization
* The iterative optimization commences until the solution satisfied the specified criterion
* The program writes the fitness value, the solution, and the computational time as .dat files
* Program terminates

Most of these functionalities are in `OptAlg` class.

For adaptive phase estimation, the program runs many consecutive optimization problems with different number of variables _N_ and the accept-reject criteria changes for different sets of _N_, which is possible by changing conditions given to the optimization algorithm in `main()`. 

## Changing the compilation setting

When a new problem and/or algorithm is included to the library, the following line in `src\Makefile.in`,

`OBJS=main.o candidate.o phase_loss_opt.o io.o problem.o mpi_optalg.o mpi_pso.o mpi_de.o candidate.o rng.o aux_functions.o`,

should be updated to include the new class.