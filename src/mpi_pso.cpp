#include "mpi_pso.h"

void PSO::write_param(double *param_array) {
    w = param_array[1];
    phi1 = param_array[2];
    phi2 = param_array[3];
    v_max = param_array[4];
    }

void PSO::put_to_best(int my_rank, int total_pop, int nb_proc) {
    int i, p;
    double array[this->num];

    for(p = 0; p < this->pop_size; ++p) {
        this->pop[p].update_best();
        this->pop[p].init_velocity();
        for(i = 0; i < this->num; ++i) {
            array[i] = double(rand()) / RAND_MAX * v_max;
            }
        this->pop[p].update_vel(array);
        }

    find_global(my_rank, total_pop, nb_proc);

    }

void PSO::find_global(int my_rank, int total_pop, int nb_proc) {
    MPI_Status status;
    int tag = 1;
    int p, ptr, prev, forw, i;
    double prev_fit, forw_fit, fit;
    double fitarray[this->prob->num_fit];
    double array[this->num];

    for(p = 0; p < total_pop; ++p) {

        //find index of the candidate in the neighborhood
        if(p == 0) {
            prev = total_pop - 1;
            }
        else {
            prev = p - 1;
            }
        forw = (p + 1) % total_pop;

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

        //fitness values sent to p
        if(my_rank == p % nb_proc) {
            fit = this->pop[p / nb_proc].read_bestfit(0); //read fitness of p
            ptr = prev;
            if(prev_fit <= fit) {
                ptr = p;
                if(fit < forw_fit) {
                    ptr = forw;
                    }
                else {}   //ptr=p
                }
            else if(prev_fit > fit) {
                //ptr=prev
                if(prev_fit < forw_fit) {
                    ptr = forw;
                    }
                else {}   //ptr still points to prev
                }
            else {
                cout << "Error!!!" << endl;
                }
            }

        //send ptr
        if(my_rank == p % nb_proc) {
            MPI_Send(&ptr, 1, MPI_INT, prev % nb_proc, tag, MPI_COMM_WORLD);
            }
        else if(my_rank == prev % nb_proc) {
            MPI_Recv(&ptr, 1, MPI_INT, p % nb_proc, tag, MPI_COMM_WORLD, &status);
            //cout<<",prev ptr="<<ptr;
            }
        else {}



        if(my_rank == p % nb_proc) {
            MPI_Send(&ptr, 1, MPI_INT, forw % nb_proc, tag, MPI_COMM_WORLD);
            }
        else if(my_rank == forw % nb_proc) {
            MPI_Recv(&ptr, 1, MPI_INT, p % nb_proc, tag, MPI_COMM_WORLD, &status);
            //cout<<",forw ptr="<<ptr<<endl;
            }
        else {}

        MPI_Barrier(MPI_COMM_WORLD);

        //updating global best
        if(ptr == p) { //Global best already in precessor's memory. No need for MPI
            if(my_rank == ptr % nb_proc) {
                this->pop[p / nb_proc].read_best(array);
                this->pop[p / nb_proc].update_global(array);
                for(i = 0; i < this->prob->num_fit; i++) {
                    fitarray[i] = this->pop[p / nb_proc].read_bestfit(i);
                    }
                this->pop[p / nb_proc].write_globalfit(fitarray);
                }
            else {}
            //cout<<"ptr=p:"<<ptr<<","<<p<<endl;
            }
        else if(ptr == prev || ptr == forw) {
            //cout<<"ptr:"<<ptr<<",prev="<<prev<<",forw="<<forw<<endl;
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
                for(i = 0; i < this->prob->num_fit; i++) {
                    fitarray[i] = this->pop[ptr / nb_proc].read_globalfit(i);
                    }
                MPI_Send(&fitarray[0], this->prob->num_fit, MPI_DOUBLE, p % nb_proc, tag, MPI_COMM_WORLD);
                }
            else if(my_rank == p % nb_proc) {
                MPI_Recv(&fitarray[0], this->prob->num_fit, MPI_DOUBLE, ptr % nb_proc, tag, MPI_COMM_WORLD, &status);
                this->pop[p / nb_proc].write_globalfit(fitarray);
                }
            else {}
            }
        else {
            //cout<<"ptr:"<<ptr<<",prev="<<prev<<",forw="<<forw<<endl;
            //cout<<"Can't find global. Code Error"<<endl;
            }
        //end updating global best

        }//p loop
    }

void PSO::combination(int my_rank, int total_pop, int nb_proc) {
    int p, i;
    double global_pos[this->num];
    double personal_pos[this->num];
    double pos[this->num];
    double vel[this->num];
    double new_pos[this->num];

    for(p = 0; p < this->pop_size; ++p) {
        srand(time(NULL) + p);
        i = rand(); //flushing out the deterministic results
        this->pop[p].read_global(global_pos);
        this->pop[p].read_best(personal_pos);
        this->pop[p].read_cont(pos);
        this->pop[p].read_vel(vel);
        for(i = 0; i < this->num; ++i) {
            new_pos[i] = pos[i] + vel[i];
            vel[i] = vel[i] + phi1 * double(rand()) / RAND_MAX * (personal_pos[i] - pos[i]) + phi2 * double(rand()) / RAND_MAX * (global_pos[i] - pos[i]);
            if(vel[i] > v_max) {
                vel[i] = v_max;
                }
            else if(vel[i] < -v_max) {
                vel[i] = -v_max;
                }
            else {
                /*keep vel[i]*/
                }
            }
        //this->prob->normalize(new_pos);
        this->prob->modulo(new_pos);
        this->pop[p].update_cont(new_pos);
        this->pop[p].update_vel(vel);
        this->Cont_fitness(p);
        }
    }

void PSO::selection(int my_rank, int total_pop, int nb_proc) {
    int p;
    for(p = 0; p < this->pop_size; ++p) {
        if(this->pop[p].read_bestfit(0) < this->pop[p].read_contfit(0)) {
            this->pop[p].update_best();
            }
        else {
            /*do nothing cout<<endl;*/
            }
        }
    find_global(my_rank, total_pop, nb_proc);

    }
