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

SUITE(pso_construct) {
    TEST(check_num) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new PSO(problem);
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
        OptAlg* alg = new PSO(problem);
        CHECK(alg->num_fit == problem->num_fit);
        }
    }
SUITE(test_param) {
    TEST(read_param) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new PSO(problem);

        double total = 0;
        double input[4] = {0.8, 0.6, 1.0, 0.2};
        double output[4];
        alg->read_param(output);
        for(int i = 0; i < 4; ++i) {
            total += input[i] - output[i];
            }

        CHECK(total == 0);
        }

    TEST(write_param) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        OptAlg* alg = new PSO(problem);

        double total = 0;
        double input[4] = {1.7, 2.6, 5.6, 0.9};
        double output[4];
        alg->write_param(input);
        alg->read_param(output);
        for(int i = 0; i < 4; ++i) {
            total += input[i] - output[i];
            }
        CHECK(total == 0);
        }

    }
SUITE(pso_findglobal) {
    TEST(index_p0) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        PSO* alg = new PSO(problem);

        int pop_size = 5;
        int prev, forw;
        int p = 0;
        alg->find_index(&prev, &forw, p, pop_size);
        CHECK(prev == pop_size - 1);
        }
    TEST(index_popsize) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        PSO* alg = new PSO(problem);

        int pop_size = 5;
        int prev, forw;
        int p = pop_size - 1;
        alg->find_index(&prev, &forw, p, pop_size);
        CHECK(forw == 0);
        }
    TEST(index_fitness) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        PSO* alg = new PSO(problem);

        int pop_size = 5;
        int prev, forw, ptr;
        int p = pop_size - 1;
        alg->find_index(&prev, &forw, p, pop_size);
        double prev_fit = 0.1 * (p + 1);
        double fit = 0.1 * p;
        double forw_fit = 0.1 * (p - 1);
        ptr = alg->find_fitness(prev, prev_fit, forw, forw_fit, p, fit);
        CHECK(ptr == prev);
        }
    /*    TEST(find_global) {
            int xn = 10;
            int xu = 10;
            int seed = 0;
            int rank = 0;
            Rng *rng = new Rng(xn, xu, seed, rank);
            int numvar = 8;
            Problem* problem = new Phase(numvar, rng);
            PSO* alg = new PSO(problem);

            int pop_size = 5;
            int my_rank = 0;
            int nb_proc = 1;
            double fit[alg->num_fit];
            alg->Init_population(pop_size);
            for(int i = 0; i < pop_size; ++i) {
                fit[0] = 0.1 * i;
                fit[1] = i;
                alg->pop[i].write_contfit(fit, 1);
                }
            alg->put_to_best(my_rank, pop_size, nb_proc);
            alg->find_global(my_rank, pop_size, nb_proc);
            CHECK(alg->pop[2].read_globalfit(0) == 0.1 * (0));
            }*/
    }
SUITE(pso_functions) {
    TEST(put_to_best) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        PSO* alg = new PSO(problem);

        int pop_size = 5;
        int my_rank = 0;
        int nb_proc = 1;
        double fit[alg->num_fit];
        alg->Init_population(pop_size);
        for(int i = 0; i < pop_size; ++i) {
            fit[0] = i * 0.1;
            fit[1] = i;
            alg->pop[i].write_contfit(fit, 1);
            }
        alg->put_to_best(my_rank, pop_size, nb_proc);
        CHECK(alg->pop[pop_size - 1].read_bestfit(0) == 0.1 * (pop_size - 1));
        }
    TEST(selection) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        PSO* alg = new PSO(problem);

        int pop_size = 5;
        int my_rank = 0;
        int nb_proc = 1;
        double fit[alg->num_fit];
        alg->Init_population(pop_size);
        for(int i = 0; i < pop_size; ++i) {
            fit[0] = i * 0.1;
            fit[1] = i;
            alg->pop[i].write_contfit(fit, 1);
            }
        alg->put_to_best(my_rank, pop_size, nb_proc);
        for(int i = 0; i < pop_size; ++i) {
            fit[0] = (pop_size - 1) * 0.1;
            fit[1] = i;
            alg->pop[i].write_contfit(fit, 1);
            }
        alg->selection(my_rank, pop_size, nb_proc);
        CHECK(alg->pop[0].read_bestfit(0) == 0.1 * (pop_size - 1));
        }
    TEST(combination) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 8;
        Problem* problem = new Phase(numvar, rng);
        PSO* alg = new PSO(problem);

        int pop_size = 5;
        int my_rank = 0;
        int nb_proc = 1;
        double fit[alg->num_fit];
        double velo[alg->num];
        double out[alg->num];
        double total = 0;
        alg->Init_population(pop_size);
        for(int i = 0; i < pop_size; ++i) {
            fit[0] = i * 0.1;
            fit[1] = i;
            alg->pop[i].write_contfit(fit, 1);
            }
        alg->put_to_best(my_rank, pop_size, nb_proc);
        for(int p = 0; p < pop_size; ++p) {
            alg->pop[p].put_to_global();
            }
        alg->pop[0].read_vel(velo);
        alg->combination(my_rank, pop_size, nb_proc);
        alg->pop[0].read_vel(out);
        for(int i = 0; i < pop_size; ++i) {
            total += out[i] - velo[i];
            }
        CHECK(total == 0);
        }
    }
