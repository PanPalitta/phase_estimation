#include "mpi_optalg.h"

void DE::put_to_best() {
    for(int p = 0; p < this->pop_size; ++p) { //pop_size here is the number of candidates assigned for a processor from the initialization.
        this->pop[p].update_best();
        }
    }

void DE::write_param(double *param_array) {
    F = param_array[0];
    Cr = param_array[1];
    }

void DE::read_param(double *param_array) {
    param_array[0] = F;
    param_array[1] = Cr;
    }

void DE::fit_to_global() {
    for(int p = 0; p < this->pop_size; ++p) {
        this->pop[p].put_to_global();
        }
    }

void DE::selection() {
    /*! Selection() reads the mean fitness value of best array and contender and decides which one is to be the member of the population.
    In this version, we simple select the one with the higher fitness value.
    */
    for(int p = 0; p < this->pop_size; ++p) {
        if(this->pop[p].read_bestfit(0) < this->pop[p].read_contfit(0)) {
            this->pop[p].update_best();
            }
        }
    }

void DE::combination() {
    /*! Combination() generates the next generation of candidates. First, all the candidate give their best array to the zeroth processor
    * who then generates the new candidate and send it back to the parent candidate. The parent store the new array into the contender
    * and compute the mean fitness value.
    */
    MPI_Status status;
    int tag = 1;
    double coin;
    int fam_size = 3;
    double all_soln[total_pop][this->num];
    int fam[fam_size];
    double input[this->num];

    //get the solution from all processor to root ***POTENTIAL BOTTLENECK***
    for(int p = 0; p < total_pop; ++p) {
        if(p % nb_proc != 0) { //if candidate is not in root, send the solution to root.
            if(my_rank == p % nb_proc) {
                this->pop[int(p / nb_proc)].read_best(input);
                MPI_Send(&input, this->num, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
                }
            else if(my_rank == 0) {
                MPI_Recv(&input, this->num, MPI_DOUBLE, p % nb_proc, tag, MPI_COMM_WORLD, &status);
                for(int i = 0; i < this->num; ++i) {
                    all_soln[p][i] = input[i];
                    }
                }
            }
        else if(p % nb_proc == 0) { //if candidate is in root, read the solution into the memory
            if(my_rank == 0) {
                this->pop[int(p / nb_proc)].read_best(input);
                for(int i = 0; i < this->num; ++i) {
                    all_soln[p][i] = input[i];
                    }
                }
            }
        }// p loop

    MPI_Barrier(MPI_COMM_WORLD);

    //generate the new candidate
    for(int p = 0; p < total_pop; ++p) {
        if(my_rank == 0) {
            family_gen(fam, p, fam_size);
            //create donor
            for(int i = 0; i < this->num; ++i) {
                //create donor
                coin = double(rand()) / RAND_MAX;
                if(coin <= Cr || i == rand() % this->num) {
                    input[i] = all_soln[fam[0]][i] + F * (all_soln[fam[1]][i] - all_soln[fam[2]][i]);
                    }
                else {
                    input[i] = all_soln[p][i];
                    }
                }
            this->prob->boundary(input);

            }//my_rank==0

        //update of candidate
        if(p % nb_proc == 0) { // if p candidate is in root, update the candidate
            if(my_rank == 0) {
                this->pop[int(p / nb_proc)].update_cont(&input[0]);
                }
            }
        else { // if p candidate is not on root, send the new candidate from root to processor
            if(my_rank == 0) {
                MPI_Send(&input[0], this->num, MPI_DOUBLE, p % nb_proc, tag, MPI_COMM_WORLD);
                }
            else if(my_rank == p % nb_proc) { //processor that contains p candidate updates the candidate
                MPI_Recv(&input[0], this->num, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
                this->pop[int(p / nb_proc)].update_cont(&input[0]);
                }
            else {}
            }

        }//p loop

    MPI_Barrier(MPI_COMM_WORLD);

    //all candidate calculate fitness of contender
    for(int p = 0; p < this->pop_size; ++p) {
        this->Cont_fitness(p);   //compute the fitness
        }

    MPI_Barrier(MPI_COMM_WORLD);
    }

void DE::family_gen(int* fam, int p, int fam_size) {
    /*! Three candidates are chosen at random from the population are selected to generate the next generation.
    * The candidates selected for the family are not allow to be the parent or the same individual more than once.
    */
    if(total_pop <= fam_size + 1) {
        for(int f = 0; f < fam_size; ++f) {
            fam[f] = rand() % total_pop;
            }
        }
    else {
        do fam[0] = rand() % total_pop;
        while(fam[0] == p);
        do fam[1] = rand() % total_pop;
        while(fam[1] == p || fam[1] == fam[0]);
        do fam[2] = rand() % total_pop;
        while(fam[2] == p || fam[2] == fam[0] || fam[2] == fam[1]);
        }
    }
