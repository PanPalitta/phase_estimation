#include "mpi_optalg.h"

/*##############################Function for initializing population##############################*/

void OptAlg::Init_population(int psize) {
    /*! This function first initialize a group of candidates for a processor uniformly over the search space.
    * The function then compute the fitness values for each.
    */
    if(psize <= 0) {
        throw invalid_argument("Population size must be positive.");
        }

    double input[num];
    //store the variables
    pop_size = psize;
    pop = new Candidate[pop_size];

    for(int p = 0; p < pop_size; ++p) {
        //generate candidate uniformly
        for(int i = 0; i < num; ++i) {
            input[i] = double(rand()) / RAND_MAX * (prob->upper_bound[i] - prob->lower_bound[i]) + prob->lower_bound[i];
            }

        //store it in candidate object
        pop[p].init_can(num, num_fit); //First the memories in the object are initialized.
        pop[p].update_cont(input); //The candidates are stored in the contender array.
        Cont_fitness(p); //The fitness values are computed.
        }
    }

void OptAlg::Init_previous(double prev_dev, double new_dev, int psize, double *prev_soln) {
    /*! This function first initialize a group of candidates for a processor from a previous solution.
    * The function then compute the fitness values for each.
    */
    if(psize <= 0) {
        throw invalid_argument("Population size must be positive");
        }
    int dev_cut_off = num - 1;
    double input[num];
    double dev[num];
    //store the variables
    pop_size = psize;
    pop = new Candidate[pop_size];

    prev_soln[num - 1] = prev_soln[num - 2]; //Generate the input array from previous solution.
    dev_gen(dev, prev_dev, new_dev, dev_cut_off); //Generate an array of the width of the Gaussian distribution.

    for(int p = 0; p < pop_size; ++p) {
        //generate candidate
        for(int i = 0; i < num; ++i) {
            input[i] = fabs(gaussian_rng->next_rand(prev_soln[i], dev[i]));
            }
        prob->boundary(input);
        //store it in candidate object
        pop[p].init_can(num, num_fit); //First the memories in the object are initialized.
        pop[p].update_cont(input); //The candidates are stored in the contender array.
        Cont_fitness(p); //The fitness values are computed.
        }

    }

/*##############################Function for calculating fitnesses##############################*/

void OptAlg::Cont_fitness(int p) {
    /*! The mean fitness value of a contender are computed for 2 samples.
    The number of samples is currently fixed.
    */
    double fit1[num_fit];
    double fit2[num_fit];

    prob->avg_fitness(pop[p].contender, prob->num_repeat, fit1);
    prob->avg_fitness(pop[p].contender, prob->num_repeat, fit2);
    for(int i = 0; i < num_fit; i++) {
        fit1[i] += fit2[i];
        fit1[i] = fit1[i] / 2.0;
        }
    pop[p].write_contfit(fit1, 2);
    }

void OptAlg::Best_fitness(int p) {
    /*! A fitness value of the best array is computed, and the mean fitness value is update.
    */
    double fit[num_fit];

    prob->avg_fitness(pop[p].can_best, prob->num_repeat, fit);
    pop[p].write_bestfit(fit);
    }

void OptAlg::update_popfit() {
    /*! The mean fitness values of the population on a processor are updated.
    */
    for(int p = 0; p < pop_size; ++p) {
        Best_fitness(p);
        }
    }

/*##############################Final Selections#################################*/
double OptAlg::Final_select(double *fit, double *solution, double *fitarray) {
    /*! The fitness values are sent to the zeroth processor, where the candidate with the highest fitness value.
    * A copy of that candidate to other processors.
    */
    int indx;
    MPI_Status status;
    int tag = 1;
    double global_fit;

    fit_to_global();//ensuring that global_best contains the solutions

    //Send all fitness value to root (p=0)

    for(int p = 0; p < total_pop; ++p) {
        if(p % nb_proc == 0) {// if p is on root, then just write the fitness to central table
            if(my_rank == 0) {
                fit[p] = pop[int(p / nb_proc)].read_globalfit(0);
                }
            }
        else { //if p is not on the root, then read the fitness and send it to central table in root
            if(my_rank == p % nb_proc) { //send the fitness from non-root
                global_fit = pop[int(p / nb_proc)].read_globalfit(0);
                MPI_Send(&global_fit, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
                }
            else if(my_rank == 0) { //root receive and store fitness value
                MPI_Recv(&fit[p], 1, MPI_DOUBLE, p % nb_proc, tag, MPI_COMM_WORLD, &status);
                }
            else {}
            }
        }

    MPI_Barrier(MPI_COMM_WORLD);

    //find the candidate that with highest fitness and send the index to all processors.
    if(my_rank == 0) {
        indx = find_max(fit);
        //root send the index to all processors.
        for(int p = 1; p < nb_proc; ++p) {
            MPI_Send(&indx, 1, MPI_INT, p, tag, MPI_COMM_WORLD);
            }
        }
    //processors other than root receive the index
    else if(my_rank != 0) {
        MPI_Recv(&indx, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
        }
    else {}

    MPI_Barrier(MPI_COMM_WORLD);

    //read the fitarray and send it to all processor for checking success criteria
    if(my_rank == indx % nb_proc) {

        for(int i = 0; i < num_fit; ++i) {
            fitarray[i] = pop[indx / nb_proc].read_globalfit(i);
            }

        for(int p = 0; p < nb_proc; ++p) {
            if(p != my_rank) {
                MPI_Send(&fitarray[0], num_fit, MPI_DOUBLE, p, tag, MPI_COMM_WORLD);
                }
            }
        }
    else if(my_rank != indx % nb_proc) {
        MPI_Recv(&fitarray[0], num_fit, MPI_DOUBLE, indx % nb_proc, tag, MPI_COMM_WORLD, &status);
        }
    else {}

    MPI_Barrier(MPI_COMM_WORLD);

    //sending the solution back to root //need to check data type
    if(indx % nb_proc == 0) { //if solution is in root, then read it out.
        if(my_rank == 0) {
            pop[int(indx / nb_proc)].read_global(solution);
            }
        else {}
        }
    else { //if solution is not in root, send to root
        if(my_rank == indx % nb_proc) {
            pop[int(indx / nb_proc)].read_global(solution);
            MPI_Send(&solution[0], num, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
            }
        else if(my_rank == 0) {
            MPI_Recv(&solution[0], num, MPI_DOUBLE, indx % nb_proc, tag, MPI_COMM_WORLD, &status);
            }
        else {}
        }

    MPI_Barrier(MPI_COMM_WORLD);

    return fitarray[0];
    }

double OptAlg::avg_Final_select(double* solution, int repeat, double *soln_fit, double *fitarray) {
    /*! Similar to Final_select(), the function finds the solution with the highest fitness in the population.
    * But before doing so, the fitness values were computed for 10 times (the variable is repeat).
    */
    MPI_Status status;
    int tag = 1;
    double final_fit;
    int indx;
    double array[num];
    double fit[pop_size][num_fit];
    //double fitarray[num_fit];

    fit_to_global();//move solution to global_best array in case we're using DE.

    //Need to calculate fitness again for 'repeat' times, independently on each
    for(int p = 0; p < pop_size; ++p) {
        //pop[p].read_global(array);
        fit[p][0] = 0;
        fit[p][1] = 0;
        for(int i = 0; i < repeat; ++i) {
            prob->avg_fitness(pop[p].global_best, prob->num_repeat, fitarray);
            fit[p][0] += fitarray[0];
	    fit[p][1] += fitarray[1];
            }
        fit[p][0] = fit[p][0] / repeat;
        fit[p][1] = fit[p][1] / repeat;
        }

    //filling the fitness table in root
    for(int p = 0; p < total_pop; ++p) {
        if(p % nb_proc == 0) {
            soln_fit[p] = fit[p / nb_proc][0]; //no need for transmitting data
            }
        else {
            if(my_rank == p % nb_proc) {
                MPI_Send(&fit[p / nb_proc][0], 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
                }
            else if(my_rank == 0) {
                MPI_Recv(&soln_fit[p], 1, MPI_DOUBLE, p % nb_proc, tag, MPI_COMM_WORLD, &status);
                }
            else {}
            }
        }

    MPI_Barrier(MPI_COMM_WORLD);

// find the maximum fitness and send out the fitness for success criteria testing
    if(my_rank == 0) {
        indx = find_max(soln_fit);
        final_fit = soln_fit[indx];
        for(int p = 1; p < nb_proc; ++p) {
            MPI_Send(&final_fit, 1, MPI_DOUBLE, p, tag, MPI_COMM_WORLD);
            }
        }
    else {
        MPI_Recv(&final_fit, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
        }

    MPI_Barrier(MPI_COMM_WORLD);

// send index of solution to processors
    if(my_rank == 0) {
        for(int p = 1; p < nb_proc; ++p) {
            MPI_Send(&indx, 1, MPI_INT, p, tag, MPI_COMM_WORLD);
            }
        }
    else {
        MPI_Recv(&indx, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
        }

    MPI_Barrier(MPI_COMM_WORLD);

    //get solution from the processor
    if(my_rank == indx % nb_proc) {
        //pop[indx / nb_proc].read_global(array);
        MPI_Send(pop[indx / nb_proc].global_best, num, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
        }
    else if(my_rank == 0) {
        MPI_Recv(&solution[0], num, MPI_DOUBLE, indx % nb_proc, tag, MPI_COMM_WORLD, &status);
        }
    else {}

    MPI_Barrier(MPI_COMM_WORLD);

   //get fitarray from the processor
    if(my_rank == indx % nb_proc) {
        MPI_Send(&fit[indx / nb_proc][1], 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
        }
    else if(my_rank == 0) {
	fitarray[0]=final_fit;
        MPI_Recv(&fitarray[1], 1, MPI_DOUBLE, indx % nb_proc, tag, MPI_COMM_WORLD, &status);
        }
    else {}

    return final_fit;
    }


int OptAlg::find_max(double *fit) {
    double* soln;
    int indx;

    soln = &fit[0];
    indx = 0;
    for(int p = 1; p < total_pop; ++p) {

        if(*soln < fit[p]) {
            soln = &fit[p];
            indx = p;
            }
        else {}
        }
    return indx;
    }

/*##############################Success Criteria#################################*/

void OptAlg::set_success(int iter, bool goal_in) {
    if (iter <= 0) {
        throw out_of_range("Number of iteration has to be positive.");
        }
    success = 0;
    T = iter;
    goal = goal_in;
    }

bool OptAlg::check_success(int t, double *current_fitarray, double *memory_fitarray, int data_size, double t_goal, bool *mem_ptype, int *numvar, int N_cut, double *memory_forT) {
    /*! The condition for stopping the optimization is checked. If the condition is based on the number of iterations,
    * T_condition() is checked. This condition doesn't change t if the condition is not met, but it reverses numvar
    * so that the optimization process restarts from the initialization. This function can take the memory from several numvar.
    * If the condition is based on whether the fitness value met the error condition, error_condition() is called.
    * This function takes the data from several numvar to generate the error condition according to t_goal.
    */
    bool type;

    if(goal == 0) {

        if(t >= T) {
            return prob->T_condition(current_fitarray, numvar, N_cut, mem_ptype, memory_forT); //This is wrong
            }
        else {
            return 0;
            }
        }
    else {
        return prob->error_condition(current_fitarray, memory_fitarray, data_size, t_goal);
        }

    }


void OptAlg::dev_gen(double *dev_array, double prev_dev, double new_dev, int cut_off) {
    for(int i = 0; i < num; ++i) {
        if(i < cut_off) {
            dev_array[i] = prev_dev;
            }
        else {
            dev_array[i] = new_dev;
            }
        }
    }
