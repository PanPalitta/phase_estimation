#include "mpi_optalg.h"

void PSO::read_param(double *param_array) {
    param_array[0] = w;
    param_array[1] = phi1;
    param_array[2] = phi2;
    param_array[3] = v_max;
    }

void PSO::write_param(double *param_array) {
    w = param_array[0];
    phi1 = param_array[1];
    phi2 = param_array[2];
    v_max = param_array[3];
    }

void PSO::find_global() {
    /*! The global best position is find for a neighborhood of three. Because of the circular topology,
    * the neighbors are on other processors and the fitness values from the best array are sent back to the candidate.
    * The candidate then communicates with the processor containing the highest fitness values of the three
    * and ask for the rest of the information.
    */
    MPI_Status status;
    int tag = 1;
    int ptr, prev, forw;
    double prev_fit, forw_fit, fit;
    double fitarray[this->prob->num_fit];
    double array[this->num];

    for(int p = 0; p < total_pop; ++p) {

        //find index of the candidate in the neighborhood
        find_index(&prev, &forw, p);

        //neighbor send fitness to the processor that contains candidate p
        if(my_rank == prev % nb_proc) {
            prev_fit = this->pop[prev / nb_proc].read_bestfit(0);
            MPI_Send(&prev_fit, 1, MPI_DOUBLE, p % nb_proc, tag, MPI_COMM_WORLD);
            }
        else if(my_rank == p % nb_proc) {
            MPI_Recv(&prev_fit, 1, MPI_DOUBLE, prev % nb_proc, tag, MPI_COMM_WORLD, &status);
            }
        else {}

        if(my_rank == forw % nb_proc) {
            forw_fit = this->pop[forw / nb_proc].read_bestfit(0);
            MPI_Send(&forw_fit, 1, MPI_DOUBLE, p % nb_proc, tag, MPI_COMM_WORLD);
            }
        else if(my_rank == p % nb_proc) {
            MPI_Recv(&forw_fit, 1, MPI_DOUBLE, forw % nb_proc, tag, MPI_COMM_WORLD, &status);
            }
        else {}

        MPI_Barrier(MPI_COMM_WORLD);

        //compare fitness at the processor and store the candidate with best fitness in ptr
        if(my_rank == p % nb_proc) {
            fit = this->pop[p / nb_proc].read_bestfit(0); //read fitness of p
            ptr = find_fitness(prev, prev_fit, forw, forw_fit, p, fit);
            //send ptr to prev and forw
            MPI_Send(&ptr, 1, MPI_INT, prev % nb_proc, tag, MPI_COMM_WORLD);
            }
        else if(my_rank == prev % nb_proc) {
            MPI_Recv(&ptr, 1, MPI_INT, p % nb_proc, tag, MPI_COMM_WORLD, &status);
            }
        else {}

        if(my_rank == p % nb_proc) {
            MPI_Send(&ptr, 1, MPI_INT, forw % nb_proc, tag, MPI_COMM_WORLD);
            }
        else if(my_rank == forw % nb_proc) {
            MPI_Recv(&ptr, 1, MPI_INT, p % nb_proc, tag, MPI_COMM_WORLD, &status);
            }
        else {}

        MPI_Barrier(MPI_COMM_WORLD);

        //updating global best
        if(ptr == p) { //global best already in processor's memory. No need for MPI
            if(my_rank == ptr % nb_proc) {
                this->pop[p / nb_proc].put_to_global();
                }
            else {}
            }
        else if(ptr == prev || ptr == forw) {
            if(my_rank == ptr % nb_proc) {
                this->pop[ptr / nb_proc].read_best(array);
                MPI_Send(&array[0], this->num, MPI_DOUBLE, p % nb_proc, tag, MPI_COMM_WORLD);
                }
            else if(my_rank == p % nb_proc) {
                MPI_Recv(&array[0], this->num, MPI_DOUBLE, ptr % nb_proc, tag, MPI_COMM_WORLD, &status);
                this->pop[p / nb_proc].update_global(array);
                }
            //sending the fitarray
            if(my_rank == ptr % nb_proc) {
                for(int i = 0; i < this->num_fit; ++i) {
                    fitarray[i] = this->pop[ptr / nb_proc].read_globalfit(i);
                    }
                MPI_Send(&fitarray[0], this->num_fit, MPI_DOUBLE, p % nb_proc, tag, MPI_COMM_WORLD);
                }
            else if(my_rank == p % nb_proc) {
                MPI_Recv(&fitarray[0], this->num_fit, MPI_DOUBLE, ptr % nb_proc, tag, MPI_COMM_WORLD, &status);
                this->pop[p / nb_proc].write_globalfit(fitarray);
                }
            else {}
            }
        else {}
        //end updating global best

        }//p loop
    }
inline void PSO::find_index(int * prev, int * forw, int p) {
    /*! Finding the indexes of candidate in the neighborhood. We use the circular topology in order to determine the neighborhood.
    */
    if(p == 0) {
        *prev = total_pop - 1;
        }
    else {
        *prev = p - 1;
        }
    *forw = (p + 1) % total_pop;
    }

inline int PSO::find_fitness(int prev, double prev_fit, int forw, double forw_fit, int p, double fit) {
    /*! Find the candidate with the best fitness value from a group of three.
    */
    int ptr = prev;
    if(prev_fit <= fit) {
        ptr = p;
        if(fit < forw_fit) {
            ptr = forw;
            }
        else {}   //ptr=p
        }
    else if(prev_fit > fit) {
        if(prev_fit < forw_fit) {
            ptr = forw;
            }
        else {}   //ptr still points to prev
        }
    else {}
    return ptr;
    }

void PSO::put_to_best() {
    /*! This function in PSO record the first position as personal best, initialize the velocity and find global best position
    * in preparation for the iterative steps.
    */
    double array[this->num];

    for(int p = 0; p < this->pop_size; ++p) {
        this->pop[p].update_best();
        this->pop[p].init_velocity();
        //generating velocity
        for(int i = 0; i < this->num; ++i) {
            array[i] = double(rand()) / RAND_MAX * v_max;
            }
        this->pop[p].update_vel(array);
        }

    find_global();

    }

void PSO::selection() {
    /*!In PSO, selection chooses whether the position stored in contender is a new personal best,
    * and, if so, store it in best array. The global best position is then searched and updated.
    */
    for(int p = 0; p < this->pop_size; ++p) {
        if(this->pop[p].read_bestfit(0) < this->pop[p].read_contfit(0)) {
            this->pop[p].update_best();
            }
        else {}
        }
    find_global();
    }

void PSO::combination() {
    /*! Combination() generates the position and velocity for the next time step using the rule with inertia term.
    * The new position is stored in the contender and computed for mean fitness value.
    */
    double global_pos[this->num];
    double personal_pos[this->num];
    double pos[this->num];
    double vel[this->num];
    double new_pos[this->num];

    for(int p = 0; p < this->pop_size; ++p) {
        this->pop[p].read_global(global_pos);
        this->pop[p].read_best(personal_pos);
        this->pop[p].read_cont(pos);
        this->pop[p].read_vel(vel);
        for(int i = 0; i < this->num; ++i) {
            new_pos[i] = pos[i] + vel[i];
            vel[i] = w * vel[i] + phi1 * double(rand()) / RAND_MAX * (personal_pos[i] - pos[i]) + phi2 * double(rand()) / RAND_MAX * (global_pos[i] - pos[i]);
            if(vel[i] > v_max) {
                vel[i] = v_max;
                }
            else if(vel[i] < -v_max) {
                vel[i] = -v_max;
                }
            else {}
            }
        this->prob->boundary(new_pos);
        this->pop[p].update_cont(new_pos);
        this->pop[p].update_vel(vel);
        this->Cont_fitness(p);
        }
    }
