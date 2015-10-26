#include <C:\UnitTest++-1.3\src\UnitTest++.h> //change to UnitTest++ path on the machine.
#include <stdexcept>
#include "rng.h"
#include "problem.h"
#include "phase_loss_opt.h"
#include "mpi_optalg.h"

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
