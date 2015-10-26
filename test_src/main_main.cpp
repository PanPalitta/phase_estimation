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

SUITE(test_loss) {
    TEST(fitness_under_loss) {
        int xn = 10;
        int xu = 10;
        int seed = 0;
        int rank = 0;
        Rng *rng = new Rng(xn, xu, seed, rank);
        int numvar = 5;
        Problem* problem = new Phase(numvar, rng);

        double solution[numvar];
        double fitarray[problem->num_fit];
        for(int i = 0; i < numvar; ++i) {
            solution[i] = 0;
            }
        for(int i = 0; i < problem->num_fit; ++i) {
            fitarray[i] = 0;
            }
        problem->fitness(solution, fitarray);
        CHECK(fitarray[1] != 0);
        }
    }
