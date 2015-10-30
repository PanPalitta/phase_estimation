#include "unittest++/UnitTest++.h"
#include "phase_loss_opt.h"
#include "mpi_optalg.h"
#include "rng.h"
#include "problem.h"
#include "candidate.h"
#include "io.h"

#include <stdexcept>
#include <cstdlib>
#include <ctime>
#include <complex>

using namespace std;
typedef complex<double> dcmplx;

// run all tests
int main(int argc, char **argv) {
    return UnitTest::RunAllTests();
    }
	
//from phase_loss_opt
SUITE(phase_mainfunc) {
    TEST(constructor) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = -5;
        bool err = 0;
        Problem* problem;
        CHECK_THROW(problem = new Phase(numvar, gaussian_rng, uniform_rng), invalid_argument);
        }
    TEST(boundary) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 10;
        Phase* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        double Phi = 4 * M_PI / (numvar + 1);
        double array[numvar];
        bool err = 0;

        for(int i = 0; i < numvar; ++i) {
            array[i] = (i - numvar / 2) * Phi;
            }
        problem->boundary(array);
        for(int i = 0; i < numvar; ++i) {
            if(array[i] > problem->upper_bound[i]) {
                err = 1;
                }
            else if(array[i] < problem->lower_bound[i]) {
                err = 1;
                }
            else {}
            }

        CHECK(err == 0);
        }
    TEST(fitness) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 10;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        double can[numvar + 1];
        int K = numvar * numvar * 10;
        double fit[problem->num_fit];

        for(int i = 0; i < numvar; ++i) {
            can[i] = 1;
            }
        problem->avg_fitness(can, K, fit);
        CHECK(fit[0] < 1.0);
        }

    }
//from candidate
SUITE(candidate_initialize_memory) {
    TEST(numvar_invalid) {
        int numvar = -5;
        int num_fit = 2;
        Candidate can1;
        CHECK_THROW(can1.init_can(numvar, num_fit), out_of_range);
        }
    TEST(numfit_invalid) {
        int numvar = 4;
        int num_fit = -2;
        Candidate can1;
        CHECK_THROW(can1.init_can(numvar, num_fit), invalid_argument);
        }
    }
SUITE(candidate_reads) {
    TEST(read_vel) {
        int numvar = 4;
        int num_fit = 2;
        Candidate can1;
        can1.init_can(numvar, num_fit);
        can1.init_velocity();
        double input[numvar];
        double output[numvar];
        double total = 0;
        for(int i = 0; i < numvar; i++) {
            input[i] = 1;
            }
        can1.update_vel(input);
        can1.read_vel(output);
        for(int i = 0; i < numvar; i++) {
            total += output[i];
            }
        CHECK(total == numvar);
        }
    TEST(read_best) {
        int numvar = 4;
        int num_fit = 2;
        Candidate can1;
        can1.init_can(numvar, num_fit);
        double input[numvar];
        double output[numvar];
        double total = 0;
        for(int i = 0; i < numvar; i++) {
            input[i] = 1;
            }
        can1.update_cont(input);
        can1.update_best();
        can1.read_best(output);
        for(int i = 0; i < numvar; i++) {
            total += output[i];
            }
        CHECK(total == numvar);
        }
    }
	

//from mpi_optalg
SUITE(Init_pop) {
    TEST(check_initpop) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 10;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new OptAlg(problem, gaussian_rng);

        int pop_size = -5;
        CHECK_THROW(alg->Init_population(pop_size), invalid_argument);
        }
    TEST(check_previnit) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 5;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new OptAlg(problem, gaussian_rng);

        double prev_dev = 0.1;
        double new_dev = 0.5;
        int pop_size = 0;
        double soln[numvar - 1];
        for(int i = 0; i < numvar - 1; ++i) {
            soln[i] = 0;
            }
        CHECK_THROW(alg->Init_previous(prev_dev, new_dev, pop_size, soln), invalid_argument);
        }
    }

SUITE(success_criteria) {
    TEST(check_iter) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new OptAlg(problem, gaussian_rng);

        int iter = 0;
        double goal_in = 0.1;
        CHECK_THROW(alg->set_success(iter, goal_in), out_of_range);
        }
    TEST(check_fitsuccess) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new OptAlg(problem, gaussian_rng);

        int iter = 10;
        double goal_in = 0.1;
        alg->set_success(iter, goal_in);
        bool succ = 0;
        int t = 5;
        double slope = 1;
        double intercept = 0;
        double fit = 1 / sqrt(pow(10.0, intercept) * pow(double(numvar), slope / 2) + 1);//success due to fit value
        succ = alg->check_success(t, numvar, fit, slope, intercept);
        CHECK(succ == 1);
        }
    TEST(check_tsuccess) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new OptAlg(problem, gaussian_rng);

        int iter = 10;
        double goal_in = 0.1;
        alg->set_success(iter, goal_in);
        bool succ = 0;
        int t = 20; //success due to time up
        double slope = 1;
        double intercept = 0;
        double fit = 0.5 / sqrt(pow(10.0, intercept) * pow(double(numvar), slope / 2) + 1);
        succ = alg->check_success(t, numvar, fit, slope, intercept);
        CHECK(succ == 1);
        }
    }

SUITE(alg_aux_functions) {
    TEST(dev_gen) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new OptAlg(problem, gaussian_rng);

        int dev_cut_off = numvar - 1;
        double dev[numvar];
        double prev_dev = 0.1;
        double new_dev = 0.5;
        alg->dev_gen(dev, prev_dev, new_dev, dev_cut_off);
        CHECK(dev[numvar - 1] == new_dev);
        }
    }

SUITE(policy_criteria) {
    TEST(check_sharp) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new OptAlg(problem, gaussian_rng);

        double error = 4.0;
        double sharp = 1.0;
        bool succ;
        CHECK_THROW(alg->check_policy(error, sharp), invalid_argument);
        }

    TEST(success_policy) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new OptAlg(problem, gaussian_rng);

        double error = 4.0;
        double sharp = 0.9;
        bool succ;
        succ = alg->check_policy(error, sharp);
        CHECK(succ == 1);
        }
    TEST(success_policy2) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new OptAlg(problem, gaussian_rng);

        double error = 2.5;
        double sharp = 0.9;
        bool succ;
        succ = alg->check_policy(error, sharp);
        CHECK(succ == 0);
        }

    }

SUITE(Linear_Regression) {
    TEST(quantile) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new OptAlg(problem, gaussian_rng);

        double input = 0.975;
        double result;
        result = alg->quantile(input);
        CHECK(fabs(result - 1.957) <= 0.001);
        }
    TEST(linear_regression) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new OptAlg(problem, gaussian_rng);

        double slope, intercept, mean_x;
        double x[numvar];
        double y[numvar];
        double m = 1.0;
        double b = 1.0;
        for(int i = 0; i < numvar; i++) {
            x[i] = i;
            y[i] = m * x[i] + b;
            }
        alg->linear_fit(numvar, x, y, &slope, &intercept, &mean_x);
        CHECK((slope == m) && (intercept == b));
        }
    TEST(linear_InvalidSize) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new OptAlg(problem, gaussian_rng);

        double slope, intercept, mean_x;
        double x[numvar];
        double y[numvar];
        double m = 1.0;
        double b = 1.0;
        for(int i = 0; i < numvar; i++) {
            x[i] = i;
            y[i] = m * x[i] + b;
            }
        CHECK_THROW(alg->linear_fit(-1, x, y, &slope, &intercept, &mean_x), invalid_argument);
        }

    TEST(error_interval) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new OptAlg(problem, gaussian_rng);

        double slope, intercept, mean_x, SSres;
        double interval;
        double x[numvar];
        double y[numvar];
        double m = 1.0;
        double b = 1.0;
        for(int i = 0; i < numvar; i++) {
            x[i] = i;
            y[i] = m * x[i] + b;
            }
        alg->linear_fit(numvar, x, y, &slope, &intercept, &mean_x);
        interval = alg->error_interval(x, y, mean_x, numvar, &SSres, slope, intercept);
        CHECK(interval <= 0.000000000000000001);
        }
    TEST(interval_InvalidSize) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new OptAlg(problem, gaussian_rng);

        double slope, intercept, mean_x, SSres;
        double interval;
        double x[numvar];
        double y[numvar];
        double m = 1.0;
        double b = 1.0;
        for(int i = 0; i < numvar; i++) {
            x[i] = i;
            y[i] = m * x[i] + b;
            }
        alg->linear_fit(numvar, x, y, &slope, &intercept, &mean_x);
        CHECK_THROW(interval = alg->error_interval(x, y, mean_x, -1, &SSres, slope, intercept), invalid_argument);
        }

    TEST(error_update) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new OptAlg(problem, gaussian_rng);

        double slope, intercept, mean_x, SSres;
        double interval;
        double x[numvar];
        double y[numvar];
        double m = 1.0;
        double b = 1.0;
        for(int i = 0; i < numvar; i++) {
            x[i] = i;
            y[i] = m * x[i] + b;
            }
        alg->linear_fit(numvar - 1, x, y, &slope, &intercept, &mean_x);
        interval = alg->error_interval(x, y, mean_x, numvar - 1, &SSres, slope, intercept);
        interval = alg->error_update(numvar - 1, &SSres, &mean_x, slope, intercept, y, x);
        CHECK(interval <= 0.000000000000000001);
        }

    TEST(update_InvalidSize) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new OptAlg(problem, gaussian_rng);

        double slope, intercept, mean_x, SSres;
        double interval;
        double x[numvar];
        double y[numvar];
        double m = 1.0;
        double b = 1.0;
        for(int i = 0; i < numvar; i++) {
            x[i] = i;
            y[i] = m * x[i] + b;
            }
        alg->linear_fit(numvar - 1, x, y, &slope, &intercept, &mean_x);
        interval = alg->error_interval(x, y, mean_x, numvar - 1, &SSres, slope, intercept);
        CHECK_THROW(interval = alg->error_update(- 1, &SSres, &mean_x, slope, intercept, y, x), invalid_argument);
        }
    }

SUITE(Selections) {
    TEST(find_max) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new OptAlg(problem, gaussian_rng);

        double fit[numvar];
        int indx;
        for(int i = 0; i < numvar; ++i) {
            fit[i] = numvar * numvar / 4 - (i - numvar / 2) * (i - numvar / 2);
            }
        indx = alg->find_max(fit, numvar);
        CHECK(indx == numvar / 2);
        }
//I don't know how to test Final_select or avg_Final_select
//Finding the highest fitness has been changed to to find_max function.
    }

//from mpi_de + throw checking
SUITE(de_constructor) {
    TEST(check_num) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new DE(problem, gaussian_rng);
        CHECK(alg->num == numvar);
        }
    }
SUITE(de_functions) {
    TEST(check_readparam) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new DE(problem, gaussian_rng);

        double out_param[2];
        alg->read_param(out_param);
        CHECK(out_param[0] == 0.1);
        }
    TEST(check_writeparam) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new DE(problem, gaussian_rng);

        double in_param[2] = {1.3, 2.1};
        double out_param[2];
        alg->write_param(in_param);
        alg->read_param(out_param);
        CHECK( out_param[0] == in_param[0]);
        }
    }

SUITE(check_catch) {
    TEST(catch_numvar) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = -10;
        try {
            Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
            OptAlg* alg = new DE(problem, gaussian_rng);
            }
        catch(invalid_argument) {
            numvar = 4;
            }
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new DE(problem, gaussian_rng);
        CHECK(alg->num == 4);
        }
    TEST(catch_quantile) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new DE(problem, gaussian_rng);

        double p = 1.0;
        double result = 0;
        result = alg->quantile(p);
        CHECK(result != 0);
        }
    TEST(catch_policy) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new DE(problem, gaussian_rng);

        double error = 0.1;
        double sharp = 1.0;
        bool out = true;
        try {
            out = alg->check_policy(error, sharp);
            }
        catch(invalid_argument) {
            sharp = 0.999999;
            }
        out = alg->check_policy(error, sharp);
        CHECK(out == false);
        }
    TEST(catch_linearfit) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new DE(problem, gaussian_rng);

        int data_size = 0;
        double x[numvar];
        double y[numvar];
        double slope, intercept, mean_x;
        try {
            alg->linear_fit(data_size, x, y, &slope, &intercept, &mean_x);
            }
        catch(invalid_argument) {
            data_size = 100;
            }
        CHECK(data_size == 100);
        }
    TEST(catch_errorinterval) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new DE(problem, gaussian_rng);

        int data_size = 0;
        double x[numvar];
        double y[numvar];
        double slope, intercept, mean_x, SSres;
        try {
            alg->error_interval(x, y, mean_x, data_size, &SSres, slope, intercept);
            }
        catch(invalid_argument) {
            data_size = 100;
            }
        CHECK(data_size == 100);
        }
    TEST(catch_errorupdate) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *gaussian_rng = new Rng(true, xn, seed, rank);
        Rng *uniform_rng = new Rng(false, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, gaussian_rng, uniform_rng);
        OptAlg* alg = new DE(problem, gaussian_rng);

        int data_size = -10;
        double x[numvar];
        double y[numvar];
        double slope, intercept, mean_x, SSres;
        try {
            alg->error_update(data_size, &SSres, &mean_x, slope, intercept, y, x);
            }
        catch(invalid_argument) {
            data_size = 100;
            }
        CHECK(data_size == 100);
        }
    }

