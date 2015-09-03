#include <cmath>

#include "mpi_optalg.h"

void OptAlg::Init_population(int psize) {
    int i, p;
    double input[num];
    //store the variables
    pop_size = psize;
    pop = new Candidate<double>[pop_size];

    //srand(0);

    for(p = 0; p < pop_size; ++p) {
        //generate candidate
        for(i = 0; i < num; ++i) {
            input[i] = double(rand()) / RAND_MAX * (prob->upper_bound[i] - prob->lower_bound[i]) + prob->lower_bound[i];
            }
        //store it in candidate object
        pop[p].init_can(num, prob->num_fit);
        pop[p].update_cont(input);
        Cont_fitness(p);
        }

    }

void OptAlg::Init_previous(double prev_dev, double new_dev, int psize, double *prev_soln) {
    int i, p;
    int dev_cut_off = num - 1;
    double input[num];
    double dev[num];
    //store the variables
    pop_size = psize;
    pop = new Candidate<double>[pop_size];

    prev_soln[num - 1] = prev_soln[num - 2];
    dev_gen(dev, prev_dev, new_dev, dev_cut_off);

    //srand(0);

    for(p = 0; p < pop_size; ++p) {
        //generate candidate

        for(i = 0; i < num; ++i) {
            input[i] = abs(rand_Gaussian(prev_soln[i], dev[i]));
        }

        //prob->normalize(input);
        prob->modulo(input);
        //store it in candidate object
        pop[p].init_can(num, prob->num_fit);

        pop[p].update_cont(input);
        Cont_fitness(p);
    }

}

void OptAlg::Cont_fitness(int p) {
    int i;
    double fit1[prob->num_fit];
    double fit2[prob->num_fit];
    double soln[num];

    pop[p].read_cont(soln);
    prob->avg_fitness(soln, prob->num_repeat, fit1);
    prob->avg_fitness(soln, prob->num_repeat, fit2);
    for(i = 0; i < prob->num_fit; i++) {
        fit1[i] += fit2[i];
        }
    pop[p].write_contfit(fit1, 2);

    }

void OptAlg::Best_fitness(int p) {
    double fit[prob->num_fit];
    double soln[num];

    pop[p].read_best(soln);
    prob->avg_fitness(soln, prob->num_repeat, fit);
    pop[p].write_bestfit(fit);
    }

void OptAlg::update_popfit() {
    int p;
    for(p = 0; p < pop_size; ++p) {
        Best_fitness(p);
        }
    }
/*##############################Success Criteria#################################*/

void OptAlg::set_success(int iter, bool goal_in) {
    success = 0;
    T = iter;
    goal = goal_in;
    }

bool OptAlg::check_success(int t, int D, double fit, double slope, double intercept) {
    double fit_goal;
    double fit_del = 0;

    if(goal == 0) {
        if(t >= T) {
            return 1;
            }
        else {
            return 0;
            }
        }
    else {
        fit_goal = 1 / sqrt(pow(10.0, intercept) * pow(double(D), slope + fit_del / 2) + 1);
        if(t >= T) {
            return 1;
            }
        else {
            if(fit >= fit_goal) {
                return 1;
                }
            else {
                return 0;
                }
            }
        }
    }

void OptAlg::linear_fit(int data_size, double *x, double *y, double *slope, double *intercept, double *mean_x) {
    int v;
    double sum_x = 0;
    double sum_xx = 0;
    double sum_y = 0;
    double sum_xy = 0;

    for(v = 0; v < data_size; ++v) {
        sum_x = sum_x + x[v];
        sum_xx = sum_xx + x[v] * x[v];
        sum_y = sum_y + y[v];
        sum_xy = sum_xy + x[v] * y[v];
        *mean_x = *mean_x + x[v];
        }

    *mean_x = *mean_x / double(data_size);

    *slope = (sum_xy - sum_y * sum_x / double(data_size)) / (sum_xx - sum_x * sum_x / double(data_size));
    *intercept = sum_y / double(data_size) - *slope * sum_x / double(data_size);

    /*error in slope*/
    /*
    double res_y,res_x;
    double fit_del;
    res_y=0;
    res_x=0;
    for(v=0;v<data_size;v++){
    res_y=res_y+pow(y[v]-intercept-slope*x[v],2);
    res_x=res_x+pow((x[v]-mean_x),2);
    }
    fit_del=sqrt(res_y/res_x/double(data_size-2));

     //extrapolation from two data points
    m=(y[numvar-1]-y[numvar])/(x[numvar-1]-x[numvar]);
    b=y[numvar]-m*x[numvar];*/
    //printf("m=%lf, fit_del=%lf, b=%lf\n",m,fit_del,b);

    }

double OptAlg::error_interval(double *x, double *y, double mean_x, int data_size, double *SSres, double slope, double intercept) {
    int i;
    double SSx = 0;

    for(i = 0; i < data_size; ++i) {
        *SSres = *SSres + pow(y[i] - slope * x[i] - intercept, 2);
        SSx = SSx + (x[i] - mean_x) * (x[i] - mean_x);
        }

    return sqrt(*SSres / double(data_size - 2) * (1 / data_size + (pow(x[data_size - 1] - mean_x, 2) / SSx)));
    }

double OptAlg::error_update(int data_size, double *SSres, double *mean_x, double slope, double intercept, double *y, double *x) {
    int i;
    double SSx = 0;

    *mean_x = (*mean_x * data_size + x[data_size]) / double(data_size + 1);
    *SSres = *SSres + pow(y[data_size] - slope * x[data_size] - intercept, 2);

    for(i = 0; i < data_size + 1; ++i) {
        SSx = SSx + (x[i] - *mean_x) * (x[i] - *mean_x);
        }

    return sqrt(*SSres / double(data_size - 1) * (1 / data_size + (pow(x[data_size] - *mean_x, 2) / SSx)));
    }


/*quantile calculation*/

inline int OptAlg::sgn(double x) {
    if(x < 0) {
        return -1;
        }
    else if(x == 0) {
        return 0;
        }
    else {
        return 1;
        }
    }

inline double OptAlg::inv_erf(double x) {
    double a = 0.140012;
    double lnx = log(1 - x * x);
    double temp = sqrt(pow(2.0 / (M_PI * a) + lnx / 2.0, 2) - lnx / a);

    return sgn(x) * sqrt(temp - 2.0 / (M_PI * a) - lnx / 2.0);
    }

double OptAlg::quantile(double p) { //p is percentile
    return sqrt(2) * inv_erf(2 * p - 1);
    }
/*##############################Policy Type#################################*/
bool OptAlg::check_policy(double error, double sharp) {
    double sd = sqrt(1 / (sharp * sharp) - 1);
    if(error >= M_PI - sd) {
        return 1;
        }
    else {
        return 0;
        }
    }



/*##############################Final Selections#################################*/
double OptAlg::Final_select(int my_rank, int total_pop, int nb_proc, double *fit, double *solution, double *fitarray) {
    double *soln;
    int p, i;
    int indx;
    MPI_Status status;
    int tag = 1;
    double global_fit;

    fit_to_global();//ensuring that global_best contains the solutions

    //roots needed all the global fitness

    for(p = 0; p < total_pop; ++p) {
        if(p % nb_proc == 0) {
            if(my_rank == 0) {
                fit[p] = pop[int(p / nb_proc)].read_globalfit(0);
                }
            }
        else {
            if(my_rank == p % nb_proc) {
                global_fit = pop[int(p / nb_proc)].read_globalfit(0);
                MPI_Send(&global_fit, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);

                }
            else if(my_rank == 0) {
                MPI_Recv(&fit[p], 1, MPI_DOUBLE, p % nb_proc, tag, MPI_COMM_WORLD, &status);
                }
            else {}
            }
        }


    MPI_Barrier(MPI_COMM_WORLD);

    //find the candidate that is the solution and send the index to all processor.
    if(my_rank == 0) {
        soln = &fit[0];
        indx = 0;
        for(p = 1; p < total_pop; ++p) {

            if(*soln < fit[p]) {
                soln = &fit[p];
                indx = p;
                }
            else {}
            }


        for(p = 1; p < nb_proc; ++p) {
            MPI_Send(&indx, 1, MPI_INT, p, tag, MPI_COMM_WORLD);
            }
        }
    else if(my_rank != 0) {
        MPI_Recv(&indx, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
        }
    else {}

    MPI_Barrier(MPI_COMM_WORLD);

    //read the fitarray and send it to all processor
    if(my_rank == indx % nb_proc) {

        for(i = 0; i < prob->num_fit; i++) {
            fitarray[i] = pop[indx / nb_proc].read_globalfit(i);
            }

        for(p = 0; p < nb_proc; p++) {
            if(p != my_rank) {
                MPI_Send(&fitarray[0], prob->num_fit, MPI_DOUBLE, p, tag, MPI_COMM_WORLD);
                }
            }
        }
    else if(my_rank != indx % nb_proc) {
        MPI_Recv(&fitarray[0], prob->num_fit, MPI_DOUBLE, indx % nb_proc, tag, MPI_COMM_WORLD, &status);
        }
    else {}

    MPI_Barrier(MPI_COMM_WORLD);

    //sending the solution back to root //need to check data type
    if(indx % nb_proc == 0) {
        if(my_rank == 0) {
            pop[int(indx / nb_proc)].read_global(solution);
            }
        else {}
        }
    else {
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

    return global_fit;
    }

double OptAlg::avg_Final_select(double* solution, int repeat, int my_rank, int total_pop, int nb_proc, double *soln_fit) {
    MPI_Status status;
    int tag = 1;
    double *soln;
    double final_fit;
    int p, i, indx;
    double array[num];
    double fit[pop_size];
    double fitarray[prob->num_fit];

    fit_to_global();//move solution to global_best array in case we're using DE.

    //Need to calculate fitness again for 'repeat' times, independently on each
    for(p = 0; p < pop_size; ++p) {
        pop[p].read_global(array);
        fit[p] = 0;
        for(i = 0; i < repeat; ++i) {
            prob->avg_fitness(array, prob->num_repeat, fitarray);
            fit[p] += fitarray[0];
            }
        fit[p] = fit[p] / repeat;
        }

    //filling the fitness table in root
    for(p = 0; p < total_pop; ++p) {
        if(p % nb_proc == 0) {
            soln_fit[p] = fit[p / nb_proc]; //no need for transmitting data
            }
        else {
            if(my_rank == p % nb_proc) {
                MPI_Send(&fit[p / nb_proc], 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
                }
            else if(my_rank == 0) {
                MPI_Recv(&soln_fit[p], 1, MPI_DOUBLE, p % nb_proc, tag, MPI_COMM_WORLD, &status);
                }
            else {}
            }

        }

    MPI_Barrier(MPI_COMM_WORLD);

    if(my_rank == 0) {
        soln = &soln_fit[0];
        indx = 0;
        for(p = 1; p < total_pop; ++p) {
            if(*soln < soln_fit[p]) {
                soln = &soln_fit[p];
                indx = p;
                }
            else {}
            }
        final_fit = *soln;
        for(p = 1; p < nb_proc; ++p) {
            MPI_Send(&final_fit, 1, MPI_DOUBLE, p, tag, MPI_COMM_WORLD);
            }
        }
    else {
        MPI_Recv(&final_fit, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
        }

    MPI_Barrier(MPI_COMM_WORLD);

    if(my_rank == 0) {
        for(p = 1; p < nb_proc; ++p) {
            MPI_Send(&indx, 1, MPI_INT, p, tag, MPI_COMM_WORLD);
            }
        }
    else {
        MPI_Recv(&indx, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
        }

    MPI_Barrier(MPI_COMM_WORLD);

    //get solution from the processor
    if(my_rank == indx % nb_proc) {
        pop[indx / nb_proc].read_global(array);
        MPI_Send(&array[0], num, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
        }
    else if(my_rank == 0) {
        MPI_Recv(&array[0], num, MPI_DOUBLE, indx % nb_proc, tag, MPI_COMM_WORLD, &status);
        }
    else {}

    MPI_Barrier(MPI_COMM_WORLD);

    if(my_rank == 0) {
        for(i = 0; i < num; ++i) {
            solution[i] = array[i];
            }
        }

    return final_fit;
    }

/*####################Auxilliary functions####################*/

double OptAlg::rand_Gaussian(double mean, /*the average theta*/
                                    double dev /*deviation for distribution*/
                                   ) {
    /*creating random number using Box-Muller Method/Transformation*/
    //double Z0;//,Z1;
    //double U1,U2; /*uniformly distributed random number input*/
    //double r;
    /*create input between [-1,1]*/
    /*do {
        U1=2.0*double(rand())/RAND_MAX-1.0;
        U2=2.0*double(rand())/RAND_MAX-1.0;
        r=U1*U1+U2*U2;
    } while(r==0.0||r>=1.0);*/
    /*using Box-Muller Transformation*/
    //Z0=U1*sqrt(-2*log(r)/r);
    //return Z0*dev+mean;

    //Approximating the Gaussian distribution with uniform distribution
    return (double(rand()) / RAND_MAX - 0.5) * 2 * dev + mean;
    }/*end of rand_Gaussian*/

void OptAlg::dev_gen(double *dev_array, double prev_dev, double new_dev, int cut_off) {
    int i;
    for(i = 0; i < num; ++i) {
        if(i < cut_off) {
            dev_array[i] = prev_dev;
            }
        else {
            dev_array[i] = new_dev;
            }
        }
    }
