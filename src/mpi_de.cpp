#include "mpi_optalg.h"

void DE::put_to_best(int my_rank, int total_pop, int nb_proc) {
    int p;
    for(p = 0; p < this->pop_size; ++p) { //pop_size here is the number of candidates assigned for a processor from the initialization.
        this->pop[p].update_best();
        }
    }

void DE::write_param(double *param_array) {
    F = param_array[1];
    Cr = param_array[2];
    }

void DE::combination(int my_rank, int total_pop, int nb_proc) {
    MPI_Status status;
    int tag = 1;
    int i, f, p;
    double coin;
    int fam_size = 3;
    double all_soln[total_pop][this->num];
    int fam[fam_size];
    double input[this->num];
    double can[this->num];

    srand(0 + my_rank);

    //get the solution from all processor to root ***POTENTIAL BOTTLENECK***
    for(p = 0; p < total_pop; ++p) {
        if(p % nb_proc != 0) { //if candidate is not in root, send the solution to root.
            if(my_rank == p % nb_proc) {
                this->pop[int(p / nb_proc)].read_best(input);
                MPI_Send(&input, this->num, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
                }
            else if(my_rank == 0) {
                MPI_Recv(&input, this->num, MPI_DOUBLE, p % nb_proc, tag, MPI_COMM_WORLD, &status);
                for(i = 0; i < this->num; ++i) {
                    all_soln[p][i] = input[i];
                    }
                }
            }
        else if(p % nb_proc == 0) { //if candidate is in root, read the solution into the memory
            if(my_rank == 0) {
                this->pop[int(p / nb_proc)].read_best(input);
                for(i = 0; i < this->num; ++i) {
                    all_soln[p][i] = input[i];
                    }
                }
            }
        }// p loop

    MPI_Barrier(MPI_COMM_WORLD);

    //select family members for the candidate on the processor other that root
    for(p = 0; p < total_pop; ++p) {
        if(my_rank == 0) {
            //srand(time(NULL)+p);
            //temp=rand();//flush out the first sampling which is not random
            if(total_pop <= fam_size + 1) {
                for(f = 0; f < fam_size; ++f) {
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

            //cout<<"family="<<fam[0]<<","<<fam[1]<<","<<fam[2]<<endl;

            //create donor
            for(i = 0; i < this->num; ++i) {
                //create donor
                input[i] = all_soln[fam[0]][i] + F * (all_soln[fam[1]][i] - all_soln[fam[2]][i]);
                //cross-over
                coin = double(rand()) / RAND_MAX;
                if(coin <= Cr || i == rand() % this->num) {
                    can[i] = input[i];
                    }
                //if(coin<=Cr) {
                //    can[i]=input[i];
                //}
                else {
                    can[i] = all_soln[p][i];
                    }
                }
            //send the new candidate back
            MPI_Send(&can, this->num, MPI_DOUBLE, p % nb_proc, tag, MPI_COMM_WORLD);
            }//my_rank==0

        if(my_rank == p % nb_proc) { //receive and store the new candidates in contender
            MPI_Recv(&can, this->num, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
            //this->prob->normalize(input);
            this->prob->boundary(can);
            this->pop[int(p / nb_proc)].update_cont(can);
            }

        }//p loop

    MPI_Barrier(MPI_COMM_WORLD);

    for(p = 0; p < this->pop_size; ++p) {
        this->Cont_fitness(p);   //compute the fitness
        }

    MPI_Barrier(MPI_COMM_WORLD);
    }

void DE::selection(int my_rank, int total_pop, int nb_proc) {
    int p;
    for(p = 0; p < this->pop_size; ++p) {
        if(this->pop[p].read_bestfit(0) < this->pop[p].read_contfit(0)) {
            this->pop[p].update_best();
            //cout<<"can select:"<<this->pop[p].read_bestfit()<<endl;
            }
        else {
            /*do nothing cout<<endl;*/
            }
        }
    };

void DE::fit_to_global() {
    int p;
    for(p = 0; p < this->pop_size; ++p) {
        this->pop[p].put_to_global();
        //cout<<p<<":"<<this->pop[p].read_globalfit(0)<<","<<this->pop[p].read_globalfit(1)<<endl;
        }
    }
