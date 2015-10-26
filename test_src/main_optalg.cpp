#include <C:\UnitTest++-1.3\src\UnitTest++.h>
#include <stdexcept>
#include "rng.h"
#include "problem.h"
#include "phase_loss_opt.h"
#include "mpi_optalg.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// To add a test, simply put the following code in the a .cpp file of your choice:
//
// =================================
// Simple Test
// =================================
//
//  TEST(YourTestName)
//  {
//  }
//
// The TEST macro contains enough machinery to turn this slightly odd-looking syntax into legal C++, and automatically register the test in a global list.
// This test list forms the basis of what is executed by RunAllTests().
//
// If you want to re-use a set of test data for more than one test, or provide setup/teardown for tests,
// you can use the TEST_FIXTURE macro instead. The macro requires that you pass it a class name that it will instantiate, so any setup and teardown code should be in its constructor and destructor.
//
//  struct SomeFixture
//  {
//    SomeFixture() { /* some setup */ }
//    ~SomeFixture() { /* some teardown */ }
//
//    int testData;
//  };
//
//  TEST_FIXTURE(SomeFixture, YourTestName)
//  {
//    int temp = testData;
//  }
//
// =================================
// Test Suites
// =================================
//
// Tests can be grouped into suites, using the SUITE macro. A suite serves as a namespace for test names, so that the same test name can be used in two difference contexts.
//
//  SUITE(YourSuiteName)
//  {
//    TEST(YourTestName)
//    {
//    }
//
//    TEST(YourOtherTestName)
//    {
//    }
//  }
//
// This will place the tests into a C++ namespace called YourSuiteName, and make the suite name available to UnitTest++.
// RunAllTests() can be called for a specific suite name, so you can use this to build named groups of tests to be run together.
// Note how members of the fixture are used as if they are a part of the test, since the macro-generated test class derives from the provided fixture class.
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// run all tests
int main(int argc, char **argv) {
    return UnitTest::RunAllTests();
    }

SUITE(setting_up) {
    TEST(check_num) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 10;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new OptAlg(problem);
        CHECK(alg->num == numvar);
        }
    TEST(check_numfit) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 10;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new OptAlg(problem);
        CHECK(alg->num_fit == problem->num_fit);
        }
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
    TEST(check_popsize) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 10;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new OptAlg(problem);

        int pop_size = 10;
        alg->Init_population(pop_size);
        CHECK(alg->pop_size == pop_size);
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
    TEST(check_canprev) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 5;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new OptAlg(problem);

        double total = 0.0;

        double prev_dev = 0.05;
        double new_dev = 0.05;
        int pop_size = 5;
        double soln[numvar - 1];
        for(int i = 0; i < numvar - 1; ++i) {
            soln[i] = 0;
            }
        alg->Init_previous(prev_dev, new_dev, pop_size, soln);
        for(int p = 0; p < pop_size; ++p) {
            for(int i = 0; i <= numvar; ++i) {
                total += alg->pop[p].contender[i];
                }
            }
        total = total / (pop_size * numvar);
        CHECK(fabs(total) <= (new_dev + prev_dev) / 2);
        }
    }

SUITE(Calling_fitness) {
    TEST(cont_fit) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new OptAlg(problem);
        int pop_size = 5;
        alg->Init_population(pop_size);
        alg->Cont_fitness(1);
        CHECK(alg->pop[1].times == 2);
        }
    TEST(best_fit) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new OptAlg(problem);
        int pop_size = 5;
        alg->Init_population(pop_size);
        for(int p = 0; p < pop_size; ++p) {
            alg->pop[p].update_best();
            }
        alg->Best_fitness(1);
        CHECK(alg->pop[1].best_times == 3);
        }
    TEST(update_bestfit) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new OptAlg(problem);
        int pop_size = 5;
        alg->Init_population(pop_size);
        for(int p = 0; p < pop_size; ++p) {
            alg->pop[p].update_best();
            }
        alg->update_popfit();
        CHECK(alg->pop[0].best_times == 3);
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
