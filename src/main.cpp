/*! \brief The main function contains the user's constructed code for performing optimization on the problem.
            In this particular example, the main function contains the code for optimizing feedback policy of an adaptive phase estimation scheme.
*/

#include "phase_loss_opt.h" //The header file for the specific problem
#include "mpi_optalg.h" //The header file for the optimization algorithms.'mpi.h' is included in this header.
#include "io.h" //The header file for user-specified parameters
#include "iostream"

using namespace std;

int main(int argc, char **argv) {

    /*mpi handlers*/
    int my_rank; //processor ID
    int nb_proc; //number of processors
    time_t start_time;
    MPI_Status status;
    int tag = 1;
    int *can_per_proc; //array containing the number of candidate solutions on each of the processors

    /*variables for the problem*/
    int numvar; //number of variables to be optimized
    double *solution; //array for containing the solution
    double *fitarray; //array containing the fitness values
    int p, t, T = 0;
    double final_fit; //the optimum fitness value
    double *soln_fit; //array for containing the fitness values for all candidate solutions
    bool mem_ptype[2] = {false, false}; //array for storing the type of policy -- specific to adaptive phase estimation
    double memory_forT[2]; //array for storing the fitness values for T_condition

    /*parameter settings*/
    int pop_size, iter, iter_begin, repeat, seed, data_end;
    int N_begin, N_cut, N_end; //variables for setting begin and ending of numvar
    double prev_dev, new_dev;
    double t_goal; //the acceptable level of error: in this case, it is the confidence interval
    string output_filename, time_filename, optimization;
    char const *config_filename;
    /*reading the parameters from io*/
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

    /*specifying data collection over many numvar*/
    int data_start = N_begin;
    int data_size = data_end - data_start;

    /*start mpi*/
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);

    int num_repeat = 10 * N_end * N_end; //The size of the samples used in the the learning algorithm
    soln_fit = new double[pop_size]; //an array to store global fitness from each candidate.
    solution = new double[N_end]; //an array to store solution
    double memory_fitarray[N_end - data_start + 1][2]; //memory for storing data from many numvar's to be used in accept/reject criteria

    /*Initializing the RNG*/
    //Maximum number of uniform random numbers we will use in one go
    int n_urandom_numbers = num_repeat + 2 * num_repeat * N_end;
    //Maximum number of Gaussian random numbers we will use in one go
    int n_grandom_numbers = 3 * num_repeat * N_end;
    Rng *gaussian_rng = new Rng(true, n_grandom_numbers, seed, my_rank);
    Rng *uniform_rng = new Rng(false, n_urandom_numbers, seed, my_rank);

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
    /*beginning the optimization algorithm*/
    do {
        t = 0;

        //instantiate the particular problem using pointer from Problem class
        Problem* problem;
        try {
            problem = new Phase(numvar, gaussian_rng, uniform_rng);
            }
        catch(invalid_argument) {
            numvar = 4;
            problem = new Phase(numvar, gaussian_rng, uniform_rng);
            }
        //instantiating the particular algorithm using pointer from OptAlg class
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

        /*start timer*/
        start_time = time(NULL);

        //Initializa population
        if (numvar < N_cut) {
            //If the user chose to initialize population uniformly over the search space
            try {
                opt->Init_population(can_per_proc[my_rank]);
                }
            catch(invalid_argument) {
                //cout << "Population size at processor" << my_rank << "is <=0." << endl;
                terminate();
                }
            }
        else {
            //If the user chose to initialize population using a Guassian distribution over a point in the search space.
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

        //Copy the candidates that are initialized to personal best
        opt->put_to_best();

        //setting the success criterion
        if (numvar < N_cut) {
            //The optimization terminate after T step
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
            //The optimization terminate after the set level of error is reached
            opt->set_success(1000, 1);  //stop after perform 1000 iteration or exceeding the goal.
            }
        else {
            //The optimization terminate after T step
            try {
                opt->set_success(iter, 0);
                }
            catch(out_of_range) {
                iter = 100;
                opt->set_success(iter, 0);
                }
            T = iter;
            }
        /*Start iterative optimization*/
        do {
            ++t;

            opt->update_popfit(); //recalculated the fitness values and update the mean of the fitness values for the population
            opt->combination(); //create contenders(offsprings). This is where a lot of comm goes on between processors
            opt->selection(); //select candidates or contenders for the next step

            final_fit = opt->Final_select(soln_fit, solution, fitarray); //communicate to find the best solution that exist so far

            //check if optimization is successful. This function includes accept-reject criteria.
            opt->success = opt->check_success(t, fitarray, &memory_fitarray[0][0], data_size, t_goal, mem_ptype, &numvar, N_cut, memory_forT);

            }
        while (opt->success == 0);

        //repeat calculation of fitness value to find the best one in the population
        final_fit = opt->avg_Final_select(solution, repeat, soln_fit, fitarray);

        if (my_rank == 0) {
            //if the policy is of the correct type output the solution
            if ((numvar >= N_begin) && (!mem_ptype[0] && !mem_ptype[1])) {
//		if ((numvar >= N_begin)) {
                output_result(numvar, problem->num_fit, fitarray, solution, start_time,
                              output_filename.c_str(), time_filename.c_str());
                }
            }

        //collect data to use in error calculation
        if(numvar >= data_start && numvar < data_end) {
            memory_fitarray[numvar - data_start][0] = log10(numvar);
            memory_fitarray[numvar - data_start][1] = log10(pow(final_fit, -2) - 1);
            }
        else if(numvar >= data_end) {
            ++data_size;
            }
        else {}

        //testing how the solution perform under lossy condition
        if (my_rank == 0) {
            problem->fitness(solution, fitarray);
            //final_fit = fitarray[0];
            cout << numvar << "\t" << fitarray[0] << "\t" <<fitarray[1] << endl;
            }

        /*delete the objects for algorithm and problem to free memory*/
        delete opt;
        delete problem;

        ++numvar;
        }
    while(numvar <= N_end);

    delete gaussian_rng;
    delete uniform_rng;
    delete [] solution;
    delete [] soln_fit;
    delete [] can_per_proc;

    MPI_Finalize();
    if (my_rank == 0) {
        cout << "done" << endl;
        }

    return 0;
    }
