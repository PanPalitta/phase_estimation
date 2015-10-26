#include "C:\UnitTest++-1.3\src\UnitTest++.h"
#include "candidate.h"
#include <stdexcept>

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

//Any operation on an uninitialized array will cause the test framework to fail.
//For candidate class, the delete[] in the destructor must be commented out in order to allow this test to work.

SUITE(variables) {
    TEST(check_num) {
        int input = 1;
        Candidate<double> can1;
        can1.num = input;
        CHECK(can1.num == input);
        }
    TEST(check_numfit) {
        int input = 2;
        Candidate<double> can1;
        can1.num_fit = input;
        CHECK(can1.num_fit == input);
        }
    TEST(check_time) {
        int input = 3;
        Candidate<double> can1;
        can1.times = input;
        CHECK(can1.times == input);
        }
    TEST(check_besttimes) {
        int input = 4;
        Candidate<double> can1;
        can1.best_times = input;
        CHECK(can1.best_times == input);
        }
    TEST(check_globaltimes) {
        int input = 5;
        Candidate<double> can1;
        can1.global_times = input;
        CHECK(can1.global_times == input);
        }
    TEST(check_contfit) {
        double input = 1.0;
        Candidate<double> can1;
        can1.num_fit = 2;
        can1.cont_fit = new double[can1.num_fit];
        can1.cont_fit[0] = input;
        CHECK(can1.cont_fit[0] == input);
        }
    TEST(check_bestfit) {
        double input = 2.0;
        Candidate<double> can1;
        can1.num_fit = 2;
        can1.best_fit = new double[can1.num_fit];
        can1.best_fit[0] = input;
        CHECK(can1.best_fit[0] == input);
        }
    TEST(check_globalfit) {
        double input = 2.0;
        Candidate<double> can1;
        can1.num_fit = 2;
        can1.global_fit = new double[can1.num_fit];
        can1.global_fit[0] = input;
        CHECK(can1.global_fit[0] == input);
        }
    }

SUITE(initialize_memory) {
    TEST(numvar_invalid) {
        int numvar = -5;
        int num_fit = 2;
        Candidate<double> can1;
        CHECK_THROW(can1.init_can(numvar, num_fit), out_of_range);
        }
    TEST(numfit_invalid) {
        int numvar = 4;
        int num_fit = -2;
        Candidate<double> can1;
        CHECK_THROW(can1.init_can(numvar, num_fit), invalid_argument);
        }
    TEST(check_num) {
        int numvar = 4;
        int num_fit = 2;
        Candidate<double> can1;
        can1.init_can(numvar, num_fit);
        CHECK(can1.num == numvar);
        }
    TEST(check_numfit) {
        int numvar = 4;
        int num_fit = 2;
        Candidate<double> can1;
        can1.init_can(numvar, num_fit);
        CHECK(can1.num_fit == num_fit);
        }
    }
SUITE(updates) {
    TEST(update_cont) {
        int numvar = 4;
        int num_fit = 2;
        Candidate<double> can1;
        can1.init_can(numvar, num_fit);
        double input[numvar];
        double total = 0;
        for(int i = 0; i < numvar; i++) {
            input[i] = 1;
            }
        can1.update_cont(input);
        for(int i = 0; i < numvar; i++) {
            total += can1.contender[i];
            }
        CHECK(total == numvar);
        }
    TEST(update_velocity) {
        int numvar = 4;
        int num_fit = 2;
        Candidate<double> can1;
        can1.init_can(numvar, num_fit);
        can1.init_velocity();
        double input[numvar];
        double total = 0;
        for(int i = 0; i < numvar; i++) {
            input[i] = 1;
            }
        can1.update_vel(input);
        for(int i = 0; i < numvar; i++) {
            total += can1.velocity[i];
            }
        CHECK(total == numvar);
        }
    TEST(update_best) {
        int numvar = 4;
        int num_fit = 2;
        Candidate<double> can1;
        can1.init_can(numvar, num_fit);
        double input[numvar];
        double total = 0;
        for(int i = 0; i < numvar; i++) {
            input[i] = 1;
            }
        can1.update_cont(input);
        can1.update_best();
        for(int i = 0; i < numvar; i++) {
            total += can1.can_best[i];
            }
        CHECK(total == numvar);
        }
    TEST(update_global) {
        int numvar = 4;
        int num_fit = 2;
        Candidate<double> can1;
        can1.init_can(numvar, num_fit);
        double input[numvar];
        double total = 0;
        for(int i = 0; i < numvar; i++) {
            input[i] = 1;
            }
        can1.update_global(input);
        for(int i = 0; i < numvar; i++) {
            total += can1.global_best[i];
            }
        CHECK(total == numvar);
        }
    TEST(put_to_global) {
        int numvar = 4;
        int num_fit = 2;
        Candidate<double> can1;
        can1.init_can(numvar, num_fit);
        double input[numvar];
        double total = 0;
        for(int i = 0; i < numvar; i++) {
            input[i] = 1;
            }
        can1.update_cont(input);
        can1.update_best();
        can1.put_to_global();
        for(int i = 0; i < numvar; i++) {
            total += can1.global_best[i];
            }
        CHECK(total == numvar);
        }
    }

SUITE(reads) {
    TEST(read_cont) {
        int numvar = 4;
        int num_fit = 2;
        Candidate<double> can1;
        can1.init_can(numvar, num_fit);
        double input[numvar];
        double output[numvar];
        double total = 0;
        for(int i = 0; i < numvar; i++) {
            input[i] = 1;
            }
        can1.update_cont(input);
        can1.read_cont(output);
        for(int i = 0; i < numvar; i++) {
            total += output[i];
            }
        CHECK(total == numvar);
        }
    TEST(read_vel) {
        int numvar = 4;
        int num_fit = 2;
        Candidate<double> can1;
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
        Candidate<double> can1;
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
    TEST(read_global) {
        int numvar = 4;
        int num_fit = 2;
        Candidate<double> can1;
        can1.init_can(numvar, num_fit);
        double input[numvar];
        double output[numvar];
        double total = 0;
        for(int i = 0; i < numvar; i++) {
            input[i] = 1;
            }
        can1.update_global(input);
        can1.read_global(output);
        for(int i = 0; i < numvar; i++) {
            total += output[i];
            }
        CHECK(total == numvar);
        }
    }

SUITE(fits) {
    TEST(write_contfit) {
        int numvar = 4;
        int num_fit = 2;
        Candidate<double> can1;
        can1.init_can(numvar, num_fit);

        double fit[num_fit];
        int tt = numvar;
        double total = 0;

        for(int i = 0; i < num_fit; ++i) {
            fit[i] = numvar;
            }
        can1.write_contfit(fit, tt);
        for(int i = 0; i < num_fit; ++i) {
            total += can1.cont_fit[i];
            }
        CHECK(total == num_fit);
        }
    TEST(write_globalfit) {
        int numvar = 4;
        int num_fit = 2;
        Candidate<double> can1;
        can1.init_can(numvar, num_fit);

        double fit[num_fit];
        int tt = numvar;
        double total = 0;

        for(int i = 0; i < num_fit; ++i) {
            fit[i] = 1;
            }
        can1.write_globalfit(fit);
        for(int i = 0; i < num_fit; ++i) {
            total += can1.global_fit[i];
            }
        CHECK(total == num_fit);
        }
    TEST(write_bestfit) {
        int numvar = 4;
        int num_fit = 2;
        Candidate<double> can1;
        can1.init_can(numvar, num_fit);

        double fit[num_fit];
        int tt = numvar;
        double total = 0;

        for(int i = 0; i < num_fit; ++i) {
            fit[i] = 2;
            can1.best_fit[i] = 2;
            }
        can1.best_times = 2;
        can1.write_bestfit(fit);
        for(int i = 0; i < num_fit; ++i) {
            total += can1.best_fit[i];
            }
        CHECK(total == 2 * num_fit);
        }
    }
