#include "phase_loss_opt.h"
#include "mpi_optalg.h"
#include "io.h"

using namespace std;

int main(int argc, char **argv) {

    /*mpi handlers*/
    int my_rank;
    int nb_proc;
    time_t start_time;
    MPI_Status status;
    int tag = 1;

    /*variables*/
    int numvar;
    double *solution;//the type of this array must correspond to that of the solution of the problem.
    double *fitarray;
    int p, t, T = 0;
    double final_fit;
    double *soln_fit;
    int *can_per_proc;
    double *x;
    double *y;
    bool mem_ptype[2] = {false, false};

    /*parameter settings*/
    int pop_size, N_begin, N_cut, N_end, iter, iter_begin, repeat, seed;
    string output_filename, time_filename;
    char const *config_filename;
    if (argc > 1) {
        config_filename = argv[1];
        }
    else {
        config_filename = NULL;
        }
    read_config_file(config_filename, &pop_size, &N_begin, &N_cut, &N_end, &iter,
                     &iter_begin, &repeat, &seed, &output_filename,
                     &time_filename);

    int T_cut_off = N_cut;
    double prev_dev = 0.01 * M_PI;
    double new_dev = 0.25 * M_PI;
    int data_start = N_begin;
    int data_end = 94;
    double t_goal = 0.98; //probability for calculating quantile
    int data_size = data_end - data_start;
    double slope = 0.0, intercept = 0.0;
    double mean_x = 0.0, SSres = 0.0;
    double TSSres, Tmean_x;
    double error;

    /*start mpi*/
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
    int num_repeat = 10 * N_end * N_end;
    //Maximum number of uniform random numbers we will use in one go
    int n_urandom_numbers = num_repeat + 2 * num_repeat * N_end;
    //Maximum number of Gaussian random numbers we will use in one go
    int n_grandom_numbers = 3 * num_repeat * N_end;
    Rng *rng = new Rng(n_urandom_numbers, n_grandom_numbers, seed, my_rank);
    soln_fit = new double[pop_size]; //create an array to store global fitness from each candidate.
    solution = new double[N_end];

    x = new double[N_end - data_start];
    y = new double[N_end - data_start];
    //calculating number of candidate per processor and stores the number in an array.
    can_per_proc = new int[nb_proc];
    for (p = 0; p < nb_proc; ++p) {
        can_per_proc[p] = 0; //make sure the array started with zeroes.
        }
    for (p = 0; p < pop_size; ++p) {
        can_per_proc[p % nb_proc] += 1;
        }


    if (my_rank == 0) {
        output_header(output_filename.c_str(), time_filename.c_str());
        }

    for (numvar = N_begin; numvar <= N_end; ++numvar) {
        if (my_rank == 0) {
            cout << numvar << endl;
            }
        t = 0;

        try {
            Problem* problem = new Phase(numvar, rng);
            }
        catch(invalid_argument) {
            numvar = 4;
            }
        Problem* problem = new Phase(numvar, rng);
        OptAlg* opt = new DE(problem);

        fitarray = new double[problem->num_fit];

        start_time = time(NULL);

        if (numvar < N_cut) {
            try {
                opt->Init_population(can_per_proc[my_rank]);
                }
            catch(invalid_argument) {
                can_per_proc[my_rank] = 1;
                }
            opt->Init_population(can_per_proc[my_rank]);
            }
        else {
            if (my_rank == 0) {
                for(p = 1; p < nb_proc; ++p) {
                    MPI_Send(&solution[0], numvar, MPI_DOUBLE, p, tag, MPI_COMM_WORLD);
                    }
                }
            else {
                MPI_Recv(&solution[0], numvar, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
                }
            try {
                opt->Init_previous(prev_dev, new_dev, can_per_proc[my_rank], solution);
                //each processor initialize the candidates.
                }
            catch(invalid_argument) {
                can_per_proc[my_rank] = 1;
                }
            opt->Init_previous(prev_dev, new_dev, can_per_proc[my_rank], solution);
            }

        opt->put_to_best(my_rank, pop_size, nb_proc);

        //setting the success criterion
        if (numvar < T_cut_off) {
            try {
                opt->set_success(iter_begin, 0);
                }
            catch(out_of_range) {
                iter_begin = 100;
                }
            opt->set_success(iter_begin, 0);
            T = iter_begin;
            }
        else if (numvar >= data_end) {
            opt->set_success(1000, 1);  //stop after perform 1000 iteration or exceeding the goal.
            }
        else {
            try {
                opt->set_success(iter, 0);
                }
            catch(out_of_range) {
                iter = 100;
                }
            opt->set_success(iter, 0);
            T = iter;
            }

        do {
            ++t;
            opt->update_popfit();
            opt->combination(my_rank, pop_size, nb_proc); //this is where a lot of comm goes on between processors
            opt->selection(my_rank, pop_size, nb_proc);

            //root check for success

            final_fit = opt->Final_select(my_rank, pop_size, nb_proc, soln_fit, solution, fitarray); //again, communicate to find the best solution that exist so far

            if(numvar >= data_end) {
                y[numvar - data_start] = log10(pow(final_fit, -2) - 1);
                TSSres = SSres;
                Tmean_x = mean_x;
                //error=opt->error_update(data_size,&TSSres,&Tmean_x,slope,intercept,y,x);
                //cout<<numvar<<":error="<<error<<","<<abs(y[numvar-data_start]-slope*x[numvar-data_start]-intercept)<<endl;
                }

            //checking policy type
            if(numvar == N_cut - 1) {
                if(t == T - 1) {
                    try {
                        opt->policy_type = opt->check_policy(fitarray[1], fitarray[0]);
                        }
                    catch(invalid_argument) {
                        fitarray[0] = 0.999999;
                        }
                    opt->policy_type = opt->check_policy(fitarray[1], fitarray[0]);
                    if (my_rank == 0) {
                        cout << fitarray[1] << endl;
                        }
                    mem_ptype[0] = opt->policy_type;
                    if (my_rank == 0) {
                        cout << "type[0]=" << mem_ptype[0] << endl;
                        }
                    if(mem_ptype[0] == 1) {
                        numvar = numvar - 1;
                        break;
                        }
                    }
                }
            else if(numvar == N_cut) {
                if(t == T / 3) {
                    try {
                        opt->policy_type = opt->check_policy(fitarray[1], fitarray[0]);
                        }
                    catch(invalid_argument) {
                        fitarray[0] = 0.999999;
                        }
                    opt->policy_type = opt->check_policy(fitarray[1], fitarray[0]);
                    if (my_rank == 0) {
                        cout << fitarray[1] << endl;
                        }
                    mem_ptype[1] = opt->policy_type;
                    if (my_rank == 0) {
                        cout << "type[1]=" << mem_ptype[1] << endl;
                        }
                    if(mem_ptype[0] | mem_ptype[1]) {
                        //the policy is bad
                        //reset the policy found in numvar=N_cut-1
                        if (my_rank == 0) {
                            cout << "numvar=" << numvar << " is set to";
                            }
                        numvar = N_cut - 2;
                        if (my_rank == 0) {
                            cout << numvar << endl;
                            }
                        break;
                        }
                    }
                }


            opt->success = opt->check_success(t, numvar, final_fit, slope, intercept);

            }
        while (opt->success == 0);

        final_fit = opt->avg_Final_select(solution, repeat, my_rank, pop_size, nb_proc, soln_fit);

        if (my_rank == 0) {
            if ((numvar == 4) || (!mem_ptype[0] && !mem_ptype[1])) {
                output_result(numvar, final_fit, solution, start_time,
                              output_filename.c_str(), time_filename.c_str());
                }
            }

        //collect data for linear fit
        x[numvar - data_start] = log10(numvar); //collect x data
        if (numvar >= data_start && numvar < data_end) {
            y[numvar - data_start] = log10(pow(final_fit, -2) - 1);
            }

        if (numvar == data_end - 1) {
            opt->linear_fit(data_size, x, y, &slope, &intercept, &mean_x);
            //There should be a catch for the exception here for linear_fit and error_interval,
            //but there is no easy way to keep the program running correctly if the exception is thrown.
            t_goal = (t_goal + 1) / 2;
            t_goal = opt->quantile(t_goal);
            error = opt->error_interval(x, y, mean_x, data_size, &SSres, slope, intercept);
            error = error * t_goal;
            if (my_rank == 0) {
                cout << slope << "," << intercept << endl;
                }
            //cout<<numvar<<":error="<<error<<endl;
            }
        else if (numvar >= data_end) {
            SSres = TSSres;
            mean_x = Tmean_x;
            ++data_size;
            }

        //testing loss
        if (my_rank == 0) {
            problem->fitness(solution, fitarray);
            final_fit = fitarray[0];
            if (my_rank == 0) {
                cout << numvar << "\t" << final_fit << endl;
                }
            }

        delete opt;
        delete problem;
        }
    delete rng;
    delete [] solution;
    delete [] soln_fit;
    delete [] x;
    delete [] y;
    delete [] can_per_proc;

    MPI_Finalize();
    if (my_rank == 0) {
        cout << "done" << endl;
        }

    return 0;
    }
