#ifndef CANDIDATE_H
#define CANDIDATE_H
#include <stdexcept>
#include <cstring>

using namespace std;

class Candidate {

        friend class OptAlg;

    public:
        Candidate() {
            is_velocity_initialized = false;
            is_candidate_initialized = false;
            };
        ~Candidate();

        void init_can(int numvar, int fit_size);
        void init_velocity();

        void update_cont(double *input);
        void update_vel(double *input);
        void update_best();
        void update_global(double *input);
        void put_to_global();

        void read_cont(double *output);
        void read_vel(double *output);
        void read_best(double *output);
        void read_global(double *output);

        void write_contfit(double *fit, int tt);
        void write_bestfit(double *fit);
        void write_globalfit(double *fit);

        double read_contfit(int i) {
            return cont_fit[i];
            }
        double read_bestfit(int i) {
            return best_fit[i];
            }
        double read_globalfit(int i) {
            return global_fit[i];
            }
        int read_bestt() {
            return best_times;
            }


    private:

        int num, num_fit;
        double *best_fit, *cont_fit, *global_fit;
        int times, best_times, global_times;//number of samples used to calculate average best_fit
        bool is_velocity_initialized, is_candidate_initialized;

        //memory arrays
        double *can_best;
        double *contender;
        double *velocity;
        double *global_best;

    };

#endif // CANDIDATE_H
