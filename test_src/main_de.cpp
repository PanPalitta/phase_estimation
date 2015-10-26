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

SUITE(de_constructor) {
    TEST(check_num) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new DE(problem);
        CHECK(alg->num == numvar);
        }
    TEST(check_numfit) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new DE(problem);
        CHECK(alg->num_fit == problem->num_fit);
        }
    }
SUITE(de_functions) {
    TEST(check_puttobest) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new DE(problem);

        int pop_size = 5;
        int my_rank = 0;
        int nb_proc = 0;
        int total_pop = pop_size;
        alg->Init_population(pop_size);
        alg->put_to_best(my_rank, total_pop, nb_proc);
        CHECK(alg->pop[pop_size - 1].best_fit[1] == alg->pop[pop_size - 1].cont_fit[1]);
        }
    TEST(check_readparam) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new DE(problem);

        double out_param[2];
        alg->read_param(out_param);
        CHECK(out_param[0] == 0.1);
        }
    TEST(check_writeparam) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new DE(problem);

        double in_param[2] = {1.3, 2.1};
        double out_param[2];
        alg->write_param(in_param);
        alg->read_param(out_param);
        CHECK( out_param[0] == in_param[0]);
        }
    TEST(check_puttoglobal) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new DE(problem);

        int pop_size = 5;
        int my_rank = 0;
        int nb_proc = 0;
        int total_pop = pop_size;
        alg->Init_population(pop_size);
        alg->put_to_best(my_rank, total_pop, nb_proc);
        alg->fit_to_global();
        CHECK(alg->pop[0].global_fit[0] == alg->pop[0].cont_fit[0]);
        }
    TEST(check_selection) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new DE(problem);

        int pop_size = 5;
        int my_rank = 0;
        int nb_proc = 0;
        int total_pop = pop_size;
        alg->Init_population(pop_size);
        alg->put_to_best(my_rank, total_pop, nb_proc);
        for(int p = 0; p < pop_size; ++p) {
            alg->pop[p].cont_fit[0] = 1.0;
            }
        alg->selection(my_rank, total_pop, nb_proc);
        CHECK(alg->pop[0].best_fit[0] == alg->pop[0].cont_fit[0]);
        }
    }

SUITE(combination_refactoring) {
    TEST(family_smallpop) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        DE* alg = new DE(problem);

        int pop_size = 2;
        int fam_size = 1;
        int fam[fam_size];
        alg->Init_population(pop_size);
        alg->family_gen(fam, 1, fam_size, pop_size);
        CHECK(fam[0] < pop_size);
        }
    TEST(family_bigpop) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        DE* alg = new DE(problem);

        int p = 1;
        int pop_size = 5;
        int fam_size = 1;
        int fam[fam_size];
        alg->Init_population(pop_size);
        alg->family_gen(fam, p, fam_size, pop_size);
        CHECK(fam[0] != p);
        }
//Not sure how to test the combination function itself.
    }

SUITE(check_catch) {
    TEST(catch_numvar) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = -10;
        try {
            Problem* problem = new Phase(numvar, rng);
            OptAlg* alg = new DE(problem);
            }
        catch(invalid_argument) {
            numvar = 4;
            }
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new DE(problem);
        CHECK(alg->num == 4);
        }
    TEST(catch_inverf) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new DE(problem);

        double p = 1.0;
        double result;
        try {
            result = alg->inv_erf(p);
            }
        catch(invalid_argument) {
            p = 0.999999;
            }
        CHECK(p == 0.999999);
        }
    TEST(catch_quantile) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new DE(problem);

        double p = 1.0;
        double result = 0;
        result = alg->quantile(p);
        CHECK(result != 0);
        }
    TEST(catch_initpop) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new DE(problem);

        int pop_size = -10;
        try {
            alg->Init_population(pop_size);
            }
        catch(invalid_argument) {
            if(numvar >= 50) {
                pop_size = 50;
                }
            else {
                pop_size = numvar;
                }
            }
        alg->Init_population(pop_size);
        CHECK(alg->pop_size > 0);
        }
    TEST(catch_initprev) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new DE(problem);

        int pop_size = -10;
        double prev_dev = 0.2;
        double new_dev = 0.2;
        double soln[numvar];
        try {
            alg->Init_previous(prev_dev, new_dev, pop_size, soln);
            }
        catch(invalid_argument) {
            if(numvar >= 50) {
                pop_size = 50;
                }
            else {
                pop_size = numvar;
                }
            }
        alg->Init_previous(prev_dev, new_dev, pop_size, soln);
        CHECK(alg->pop_size > 0);
        }
    TEST(catch_setsuccess) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new DE(problem);

        int iter = 0;
        bool goal = false;
        try {
            alg->set_success(iter, goal);
            }
        catch (out_of_range) {
            iter = 100;
            }
        alg->set_success(iter, goal);
        CHECK(alg->T == iter);
        }
    TEST(catch_policy) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new DE(problem);

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
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new DE(problem);

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
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new DE(problem);

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
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 4;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new DE(problem);

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
