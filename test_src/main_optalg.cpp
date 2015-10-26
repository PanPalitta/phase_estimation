#include <C:\UnitTest++-1.3\src\UnitTest++.h> //change to directory for UnitTest++
#include <stdexcept>
#include "rng.h"
#include "problem.h"
#include "phase_loss_opt.h"
#include "mpi_optalg.h"

// run all tests
int main(int argc, char **argv) {
    return UnitTest::RunAllTests();
    }

SUITE(Init_pop) {
    TEST(check_initpop) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 10;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new OptAlg(problem);

        int pop_size = -5;
        CHECK_THROW(alg->Init_population(pop_size), invalid_argument);
        }
    TEST(check_previnit) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 5;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new OptAlg(problem);

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

SUITE(success) {
    TEST(check_iter) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new OptAlg(problem);

        int iter = 0;
        double goal_in = 0.1;
        CHECK_THROW(alg->set_success(iter, goal_in), out_of_range);
        }
    TEST(check_fitsuccess) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new OptAlg(problem);

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
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new OptAlg(problem);

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

SUITE(aux_functions) {
    TEST(dev_gen) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new OptAlg(problem);

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
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new OptAlg(problem);

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
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new OptAlg(problem);

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
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new OptAlg(problem);

        double error = 2.5;
        double sharp = 0.9;
        bool succ;
        succ = alg->check_policy(error, sharp);
        CHECK(succ == 0);
        }

    }

SUITE(Linear_Regression) {
    TEST(inverf_input) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new OptAlg(problem);

        double input = 1.0;
        CHECK_THROW(alg->inv_erf(input), invalid_argument);
        }

    TEST(quantile) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new OptAlg(problem);

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
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new OptAlg(problem);

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
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new OptAlg(problem);

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
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new OptAlg(problem);

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
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new OptAlg(problem);

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
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new OptAlg(problem);

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
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new OptAlg(problem);

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
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new OptAlg(problem);

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
