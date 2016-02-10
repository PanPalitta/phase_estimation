#include "phase_loss_opt.h"
#include "mpi_optalg.h"
#include "io.h"
#include "aux_functions.h"

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
    bool mem_ptype[2] = {false, false};

    /*parameter settings*/
    int pop_size, N_begin, N_cut, N_end, iter, iter_begin, repeat, seed, data_end;
    double prev_dev, new_dev, t_goal;
    string output_filename, time_filename, optimization;
    char const *config_filename;
    if (argc > 1) {
        config_filename = argv[1];
        }
    else {
        config_filename = NULL;
        }
    read_config_file(config_filename, &pop_size, &N_begin, &N_cut, &N_end, &iter,
                     &iter_begin, &repeat, &seed, &output_filename,
                     &time_filename, &optimization, &prev_dev, &new_dev, &t_goal, &data_end);

    prev_dev = prev_dev * M_PI;
    new_dev = new_dev * M_PI;

    int data_start = N_begin;
    int data_size = data_end - data_start;

    /*start mpi*/
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
    int num_repeat = 10 * N_end * N_end;
    //Maximum number of uniform random numbers we will use in one go
    int n_urandom_numbers = num_repeat + 2 * num_repeat * N_end;
    //Maximum number of Gaussian random numbers we will use in one go
    int n_grandom_numbers = 3 * num_repeat * N_end;
    Rng *gaussian_rng = new Rng(true, n_grandom_numbers, seed, my_rank);
    Rng *uniform_rng = new Rng(false, n_urandom_numbers, seed, my_rank);
    soln_fit = new double[pop_size]; //create an array to store global fitness from each candidate.
    solution = new double[N_end];

    double memory_fitarray[2][N_end - data_start + 1];

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

    numvar = N_begin;

    do {
        if (my_rank == 0) {
            cout << numvar << endl;
            }
        t = 0;

        Problem* problem;
        try {
            problem = new Phase(numvar, gaussian_rng, uniform_rng);
            }
        catch(invalid_argument) {
            numvar = 4;
            problem = new Phase(numvar, gaussian_rng, uniform_rng);
            }
        OptAlg* opt;
        if (optimization == "de") {
            opt = new DE(problem, gaussian_rng, pop_size);
            }
        else if (optimization == "pso") {
            opt = new PSO(problem, gaussian_rng, pop_size);
            }
        else {
            throw runtime_error("Unknown optimization algorithm");
            }

        fitarray = new double[problem->num_fit];

        start_time = time(NULL);

        if (numvar < N_cut) {
            try {
                opt->Init_population(can_per_proc[my_rank]);
                }
            catch(invalid_argument) {
                cout << "Population size at processor" << my_rank << "is <=0." << endl;
                terminate();
                }
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
                cout << "Population size at processor" << my_rank << "is <=0." << endl;
                terminate();
                }
            }

        opt->put_to_best();

        //setting the success criterion
        if (numvar < N_cut) {
            try {
                opt->set_success(iter_begin, 0);
                }
            catch(out_of_range) {
                iter_begin = 100;
                opt->set_success(iter_begin, 0);
                }
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
                opt->set_success(iter, 0);
                }
            T = iter;
            }

        do {
            ++t;

            opt->update_popfit();
            opt->combination(); //this is where a lot of comm goes on between processors
            opt->selection();

            //root check for success
            final_fit = opt->Final_select(soln_fit, solution, fitarray); //again, communicate to find the best solution that exist so far

            opt->success = opt->check_success(t, fitarray, &memory_fitarray[0][0], data_size, t_goal, mem_ptype, &numvar, N_cut);
            }
        while (opt->success == 0);

        final_fit = opt->avg_Final_select(solution, repeat, soln_fit);

        if (my_rank == 0) {
            if ((numvar == 4) || (!mem_ptype[0] && !mem_ptype[1])) {
                output_result(numvar, final_fit, solution, start_time,
                              output_filename.c_str(), time_filename.c_str());
                }
            }

        //collect data for linear fit
        if(numvar >= data_start && numvar < data_end) {
            memory_fitarray[0][numvar - data_start] = log10(numvar);
            memory_fitarray[1][numvar - data_start] = log10(pow(final_fit, -2) - 1);
            }
        else if(numvar >= data_end) {
            ++data_size;
            }
        else {}

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

        ++numvar;
        }
    while(numvar <= N_end);

    delete gaussian_rng;
    delete uniform_rng;
    delete [] solution;
    delete [] soln_fit;
    delete [] memory_fitarray;
    delete [] can_per_proc;

    MPI_Finalize();
    if (my_rank == 0) {
        cout << "done" << endl;
        }

    return 0;
    }
