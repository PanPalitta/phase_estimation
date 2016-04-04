# Setting and Usage

## Usage

The program is designed to work on HPC clusters and it requires MPI to run. The basic use is as follows:

    $ [mpirun -np NPROC] phase_estimation [config_file]


## Input and configuration

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

  - `iter`: number of iterations when cluster initialization is used.
  
  - `iter_begin`: number of iteration when uniformly random initialization is used.
  
  - `repeat`: number of time the candidates are compute before the best candidate is selected after the optimization
  
  - `optimization`: choose the heuristic optimization algorithm: de (differential evolution) or pso (particle swarm optimization)
  
  - `output_filename`: the name of the file to write the results to.
  
  - `time_filename`: the name of the file to write the time taken to run the program for each number of variables.

  - `random_seed`: fix a random seed. If it is not specified, the random number generator is initialized with the system time

  - `data_end`: set where accept-reject criterion based on error from expected solution is used

  - `prev_dev`: the boundary for the cluster initialization -- for variables that are initialized from previous solution 

  - `new_dev`: the boundary for the cluster initialization -- for new variable

  - `t_goal`: parameter corresponding to the error which the algorithm will accept to solution


## Output

The program output two files, one containing the policy and fitness value (output.dat) and the other containing the CPU time used in finding the policy (time.dat). The numbers are updated for every _N_ and can be used to track the progress of the optimization.

Example of output:

output.dat

```
#N 	 Sharpness 	 Policy
4	0.854507	4.91056	5.44221	5.73905	5.87313	
5	0.888902	4.92427	5.46181	5.76412	5.88404	5.95493	
6	0.90593	4.9164	5.46353	5.75503	5.89696	5.95538	6.05514	
7	0.920337	4.91471	5.48735	5.72474	5.8765	5.94925	6.06185	6.14271	
8	0.932388	4.88864	5.45983	5.72587	5.88542	5.94729	6.06073	6.12173	6.17251	
9	0.941697	4.87251	5.45756	5.74001	5.90138	5.95748	6.04226	6.13899	6.16975	6.07994	
10	0.941703	4.86374	5.42822	5.74102	5.87777	5.96852	6.04537	6.13999	6.13744	6.0704	6.25615	

```

time.dat

```
#N 	 Time
4	1
5	0
6	1
7	1
8	0
9	0
10	0

```