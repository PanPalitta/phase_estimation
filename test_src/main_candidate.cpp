#include "C:\UnitTest++-1.3\src\UnitTest++.h" //enter directory to UnitTest++ header here 
#include "candidate.h"
#include <stdexcept>

// run all tests
int main(int argc, char **argv) {
    return UnitTest::RunAllTests();
    }

//Any operation on an uninitialized array will cause the test framework to fail.
//For candidate class, the delete[] in the destructor must be commented out in order to allow this test to work.
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
