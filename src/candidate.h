#ifndef CANDIDATE_H
#define CANDIDATE_H
#include <stdexcept>
#include <cstring>

using namespace std;

/*! \brief Candidate class contains arrays to be used by optimization algorithm to store data of a solution candidate.
*/

class Candidate {

        friend class OptAlg;

        /** OptAlg class is friend of Candidate class and can access the arrays of solution candidates directly.
                *   However, the specific algorithms are not.
                */

    public:
        Candidate() {
            is_velocity_initialized = false;
            is_candidate_initialized = false;
            };
        ~Candidate();

        void init_can(int numvar, int fit_size);/*!< A function for allocating memories for storing information for a candidate solution.*/
        void init_velocity();/*!< A function for allocating memory for additional information: velocity.*/

        void update_cont(double *input); /*!< A function for copying onto the contender array from input.*/
        void update_vel(double *input); /*!< A function for copying onto the velocity array from input.*/
        void update_best(); /*!< A function for putting the contender array to the can_best array.This is done by swapping pointers.*/
        void update_global(double *input);  /*!< A function for copying onto the global_best array from input.*/
        void put_to_global(); /*!< A function for copying can_best array to global_best array.*/

        void read_cont(double *output); /*!< A function for reading the contender array to output array.*/
        void read_vel(double *output); /*!< A function for reading the velocity array to output array.*/
        void read_best(double *output); /*!< A function for reading the can_best array to output array.*/
        void read_global(double *output); /*!< A function for reading the global_best array to output array.*/

        void write_contfit(double *fit, int tt); /*!< A function for copying an array of fitness values to cont_fit array.*/
        void write_bestfit(double *fit); /*!< A function for copying an array of fitness values to best_fit array.*/
        void write_globalfit(double *fit); /*!< A function for copying an array of fitness values to global_fit array.*/

        double read_contfit(int i /*!< index of fitness value to be read*/) { /*! A function for reading a fitness value.*/
            return cont_fit[i];
            }
        double read_bestfit(int i /*!< index of fitness value to be read*/) {/*! A function for reading a fitness value.*/
            return best_fit[i];
            }
        double read_globalfit(int i /*!< index of fitness value to be read*/) {/*! A function for reading a fitness value.*/
            return global_fit[i];
            }
        int read_bestt() {/*! A function for reading the number of times the best_fit has been averaged over.*/
            return best_times;
            }


    private:

        int num; /*number of variables*/
        int num_fit; /*number of fitness values (i.e., total number of contraints and objetive functions)*/
        double *best_fit; /*pointer to array for storing fitness values calculated from the data in contender array*/
        double *cont_fit; /*pointer to array for storing fitness values calculated from the data in can_best array*/
        double *global_fit; /*pointer to array for storing fitness values calculated from the data in global_best array*/
        int times; /*the number of time the fitess value is computed for contender*/
        int best_times; /*the number of time the fitess value is computed for can_best*/
        int global_times; /*the number of time the fitess value is computed for global_best*/
        bool is_velocity_initialized; /*indicating whether addition memory needs to be allocated*/
        bool is_candidate_initialized; /*indicating whether memories are allocated*/

        //memory arrays
        double *can_best; /*the pointer to array that contains solution that would survive several iterations*/
        double *contender; /*the pointer to array that serves mostly as temporary memory that does not survive an interation*/
        double *velocity; /*the pointer to additional array*/
        double *global_best; /*the pointer to array that contains solution to be compared for the final solution*/
    };

#endif // CANDIDATE_H
