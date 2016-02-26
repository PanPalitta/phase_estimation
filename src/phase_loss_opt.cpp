#include <iostream>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

using namespace std;

#include "phase_loss_opt.h"

Phase::Phase(const int numvar, Rng *gaussian_rng, Rng *uniform_rng):
    gaussian_rng(gaussian_rng), uniform_rng(uniform_rng) {
    if(numvar <= 0) {
        throw invalid_argument("numvar<=0. Instantiating Phase fails.");
        }

    int i;
    lower = 0;
    upper = 2 * M_PI;
    loss = 0.0;
    num = numvar;
    num_fit = 2;
    num_repeat = 10 * num * num;
    lower_bound = new double[num];
    upper_bound = new double[num];
    for(i = 0; i < num; ++i) {
        lower_bound[i] = lower;
        upper_bound[i] = upper;
        }
    sqrt_cache = new double[num + 1];
    for(i = 0; i < num + 1; ++i) {
        sqrt_cache[i] = sqrt(i);
        }
    input_state = new dcmplx[num + 1];
    sqrtfac_mat = new double[num + 1];
    overfac_mat = new double[num + 1];
    state = new dcmplx[num + 1];
    update0 = new dcmplx[num + 1];
    update1 = new dcmplx[num + 1];
    }

Phase::~Phase() {
    delete[] state;
    delete[] update0;
    delete[] update1;
    delete[] input_state;
    delete[] sqrtfac_mat;
    delete[] overfac_mat;
    }

void Phase::fitness(double *soln, double *fitarray) {
    loss = 0.2;
    avg_fitness(soln, num_repeat, fitarray);
    loss = 0.0;
    }

void Phase::avg_fitness(double *soln, const int K, double *fitarray) {
    dcmplx sharp(0.0, 0.0);
    bool dect_result;
    double PHI, phi, coin, PHI_in;
    int m, k, d;
    double error = 0.0;

    WK_state();
    for(k = 0; k < K; ++k) {
        phi = uniform_rng->next_rand(0.0, 1.0) * (upper - lower) + lower;
        PHI = 0;
        //copy input state: the optimal solution across all compilers is memcpy:
        //nadeausoftware.com/articles/2012/05/c_c_tip_how_copy_memory_quickly
        memcpy(state, input_state, (num + 1)*sizeof(dcmplx));
        //measurement
        d = 0;
        for (m = 0; m < num; ++m) {
            // This loop is the most critical part of the entire program. It
            // executes K*num=10*num^3 times on each call of avg_fitness. All
            // optimization should focus on this loop.

            //randomly decide whether loss occurs
            coin = uniform_rng->next_rand(0.0, 1.0);

            if(coin <= loss) {
                state_loss(num - m); //update only the state using loss function
                }
            else {
                PHI_in = gaussian_rng->next_rand(PHI, THETA_DEV);
                PHI_in = mod_2PI(PHI_in);//noisy PHI
                dect_result = noise_outcome(phi, PHI_in, num - m);
                if (dect_result == 0) {
                    PHI = PHI - soln[d++];
                    }
                else {
                    PHI = PHI + soln[d++];
                    }
                PHI = mod_2PI(PHI);
                }
            }
        sharp.real() = sharp.real() + cos(phi - PHI);
        sharp.imag() = sharp.imag() + sin(phi - PHI);
        error += abs(phi - PHI);
        }
    fitarray[0] = abs(sharp) / double(K);
    fitarray[1] = error / double(K);

    }

void Phase::T_condition(double *fitarray, int *numvar, int N_cut, bool *mem_ptype) {
    bool type;

    if(*numvar == N_cut - 1) {
        try {
            type = check_policy(fitarray[1], fitarray[0]);
            }
        catch(invalid_argument) {
            fitarray[0] = 0.999999;
            type = check_policy(fitarray[1], fitarray[0]);
            }
        mem_ptype[0] = type;
        if(mem_ptype[0] == 1) {
            *numvar = *numvar - 1;
            }
        }
    else if(*numvar == N_cut) {
        try {
            type = check_policy(fitarray[1], fitarray[0]);
            }
        catch(invalid_argument) {
            fitarray[0] = 0.999999;
            type = check_policy(fitarray[1], fitarray[0]);
            }
        mem_ptype[1] = type;
        if(mem_ptype[0] | mem_ptype[1]) {
            //the policy is bad
            //reset the policy found in numvar=N_cut-1
            *numvar = N_cut - 2;
            }
        }

    }

bool Phase::error_condition(double *memory_fitarray, int data_size, double t_goal) {

    double slope, intercept;
    double mean_x, SSres;
    double tn2;
    double error, error_goal;
    double x[data_size + 1];
    double y[data_size + 1];


    //split into x-y arrays

    for(int i = 0; i < data_size + 1; ++i) {
        x[i] = memory_fitarray[2 * i];
        y[i] = memory_fitarray[2 * i + 1];
        }

    linear_fit(data_size + 1, x, y, &slope, &intercept, &mean_x);

    error_goal = error_interval(x, y, mean_x, data_size + 1, &SSres, slope, intercept);

    t_goal = (t_goal + 1) / 2;
    tn2 = quantile(t_goal);
    error_goal = error_goal * tn2;

    error = y[data_size + 1] - x[data_size + 1] * slope - intercept;

    //cout<<data_size+1<<": error_goal="<<error_goal<<", error"<<error<<endl;

    if(error <= error_goal) {
        return 1;
        }
    else {
        return 0;
        }
    }



void Phase::boundary(double *can1) {
    for(int i = 0; i < num; ++i) {
        if(can1[i] < lower_bound[i]) {
            can1[i] = lower_bound[i];
            }
        else if (can1[i] > upper_bound[i]) {
            can1[i] = upper_bound[i];
            }
        }
    }

/*private functions*/
/*########### Generating Input State ###########*/
/*Generation function*/
void Phase::WK_state() {
    const double beta = M_PI / 2;
    const double cosN = pow(cos(beta / 2), num);
    tan_beta = tan(beta / 2);

    int n, k;
    double n_part[num + 1];
    double k_part[num + 1];
    double s_part;
    double temp;

    sqrtfac(sqrtfac_mat); //initializing
    one_over_fac(overfac_mat); //initializing

    /*factors in d_matrix calculated*/
    //calculating n_parts and k_parts in
    for (n = 0; n <= num; ++n) {
        n_part[n] = pow(-1.0, n) * sqrtfac_mat[n] * sqrtfac_mat[num - n] * cosN * pow(tan_beta, n);
        }
    for (k = 0; k <= num; ++k) {
        k_part[k] = 1 / pow(-1.0, k) * sqrtfac_mat[k] * sqrtfac_mat[num - k] * pow(1 / tan_beta, k);
        }

    for (n = 0; n <= num; ++n) { //we have N+1 state b/c we include n=0 and n=N.
        temp = 0;
        input_state[n].real() = 0;
        input_state[n].imag() = 0;
        for (k = 0; k <= num; ++k) {
            s_part = cal_spart(n, k, num);
            temp = s_part * k_part[k] * sin((k + 1) * M_PI / (num + 2));
            input_state[n].real() = input_state[n].real() + temp * cos(M_PI / 2.0 * (k - n));
            input_state[n].imag() = input_state[n].imag() + temp * sin(M_PI / 2.0 * (k - n));
            }//end k
        // FIXED:This correction corrects the probability overshoot.
        input_state[n].real() = input_state[n].real() * n_part[n] / sqrt(num / 2.0 + 1);
        input_state[n].imag() = input_state[n].imag() * n_part[n] / sqrt(num / 2.0 + 1);
        }//end n
    }

inline double Phase::cal_spart(const int n, const int k, const int N) {
    //calculating the Wigner d-matrix element
    int s;
    int s_min;
    int s_max;
    double s_part = 0;

    //find lower limit
    if(n - k >= 0) {
        s_min = 0;
        }
    else if(n - k < 0) {
        s_min = k - n;
        }
    else {}
    //find upper limit
    if(k <= N - n) {
        s_max = k;
        }
    else if(k > N - n) {
        s_max = N - n;
        }
    else {}

    //calculating the s_part
    for (s = s_min; s <= s_max; ++s) {
        s_part = s_part + pow(-1.0, s) * overfac_mat[k - s] * overfac_mat[s] * overfac_mat[n - k + s] * overfac_mat[N - n - s] * pow(tan_beta, 2 * s);
        }
    return s_part;
    }

inline void Phase::one_over_fac(double *over_mat) { //calculate one over factorial matrix
    over_mat[0] = 1;
    for (int i = 1; i <= num; ++i) {
        over_mat[i] = over_mat[i - 1] / i;
        }//end i
    }

inline void Phase::sqrtfac(double *fac_mat) { //calculate sqrt of factorial matrix
    //check array size
    fac_mat[0] = 1;
    for (int i = 1; i <= num; ++i) {
        fac_mat[i] = fac_mat[i - 1] * sqrt_cache[i];
        }//end i
    }

/*########### Measurement Functions ###########*/

inline bool Phase::noise_outcome(const double phi, const double PHI, const int N) {
    //N is the number of photons currently available, not equal to 'num'
    const double theta = (phi - PHI) / 2.0;
    const double cos_theta = cos(theta) / sqrt_cache[N];
    const double sin_theta = sin(theta) / sqrt_cache[N];
    //noise in operator: currently not in use
    const double oper_n0 = gaussian_rng->next_rand(0.0, DEV_N);//n_x
    const double oper_n2 = gaussian_rng->next_rand(0.0, DEV_N);//n_z
    const double oper_n1 = sqrt(1.0 - (oper_n0 * oper_n0 + oper_n2 * oper_n2));
    const dcmplx U00(sin_theta * oper_n1, -oper_n0 * sin_theta);
    const dcmplx U01(cos_theta, oper_n2 * sin_theta);
    const dcmplx U10(cos_theta, -oper_n2 * sin_theta);
    const dcmplx U11(sin_theta * oper_n1, sin_theta * oper_n0);
    int n;
    double prob = 0.0;
    for (n = 0; n < N; ++n) {
        //if C_0 is measured
        update0[n] = state[n + 1] * U00 * sqrt_cache[n + 1] + state[n] * U01 * sqrt_cache[N - n];
        prob += abs(update0[n] * conj(update0[n]));
        }

    if (uniform_rng->next_rand(0.0, 1.0) <= prob) {
        //measurement outcome is 0
        state[N] = 0;
        prob = 1.0 / sqrt(prob);
        for(n = 0; n < N; ++n) {
            state[n] = update0[n] * prob;
            }
        return 0;
        }
    else {
        //measurement outcome is 1
        prob = 0;
        for(n = 0; n < N; ++n) {
            state[n] = state[n + 1] * U10 * sqrt_cache[n + 1] - state[n] * U11 * sqrt_cache[N - n];
            prob += abs(state[n] * conj(state[n]));
            }
        state[N] = 0;
        prob = 1.0 / sqrt(prob);
        for(n = 0; n < N; ++n) {
            state[n] *= prob;
            }
        return 1;
        }
    }

inline void Phase::state_loss(const int N) {
    //state update if the photon is loss
    double total = 0;
    double factor = 1 / sqrt(2 * N);
    for(int n = 0; n < N; ++n) {
        state[n] = (state[n] * sqrt_cache[N - n] + state[n + 1] * sqrt_cache[n + 1]) * factor;
        total += state[n].real() * state[n].real() + state[n].imag() * state[n].imag();
        }
    state[N] = 0;
    //Necessary state renormalization
    for(int n = 0; n < N; ++n) {
        state[n] = state[n] / sqrt(total);
        }
    }

inline double Phase::mod_2PI(double PHI) {
    while(PHI >= 2 * M_PI) {
        PHI = PHI - 2 * M_PI;
        }
    while (PHI < 0) {
        PHI = PHI + 2 * M_PI;
        }
    return PHI;
    }

inline bool Phase::check_policy(double error, double sharp) {
    if (sharp == 1.0) {
        throw invalid_argument("sharpness cannot be one.");
        }
    double sd = sqrt(1 / (sharp * sharp) - 1);
    if(error >= M_PI - sd) {
        return 1;
        }
    else {
        return 0;
        }
    }
