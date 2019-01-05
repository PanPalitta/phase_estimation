#include <iostream>
#include <cstdlib>
#include <cstring>
#include <stdexcept>

using namespace std;

#include "phase_loss_opt.h"

Phase::Phase(const int numvar, Rng *gaussian_rng, Rng *uniform_rng):
    gaussian_rng(gaussian_rng), uniform_rng(uniform_rng) {
    /*! The instantiation function of the class initialize necessary variables that the optimization algorithm uses.
    * Memories used by the problem is also allocated in this stage.
    */
    if(numvar <= 0) {
        throw invalid_argument("numvar<=0. Instantiating Phase fails.");
        }

    int i;

    //Initializing the conditions the simulation uses.
    lower = 0;
    upper = 2 * M_PI;
    loss = 0.0;

    //Initializing the numbers used by the optimization algorithm.
    num = numvar;
    num_fit = 2;
    num_repeat = 10 * num * num;

    //Initializing the lower bound and upper bound for the optimized variables.
    lower_bound = new double[num];
    upper_bound = new double[num];
    for(i = 0; i < num; ++i) {
        lower_bound[i] = lower;
        upper_bound[i] = upper;
        }
    //Initializing memories used by the simulation.
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
    /*freeing the memory.
    * This is a crucial step if the main function runs a loop with multiple instantiation of the phase problem
    as the computer might run out of memory.
    */
    delete[] state;
    delete[] update0;
    delete[] update1;
    delete[] input_state;
    delete[] sqrtfac_mat;
    delete[] overfac_mat;
    }

void Phase::fitness(double *soln, double *fitarray) {
    /*In this particular problem, this function serves as a wrapper to change the loss rate,
    * so that we can test the policy we found in lossless interferometry with a chosen level of loss.
    */
    loss = 0.2; //This loss rate can be changed by user.
    avg_fitness(soln, num_repeat, fitarray);
    loss = 0.0; //Change back in case the optimization process has to be redone.
    }

void Phase::avg_fitness(double *soln, const int K, double *fitarray) {
    /* A function calculates the fitness values (reported in fitarray) of a solution (soln) over a sample size of K.
    * This function simulates the adaptive phase interferometry and so is the most computationally expensive part of the program.
    */
    dcmplx sharp(0.0, 0.0); //variable to store the sharpness function, which is the first fitness value in the array.
    double error = 0.0; //variable to store the bias of the estimate, which is the second fitness value in the array.
    bool dect_result; //variable to store which path a photon comes out at any step.
    double PHI, phi, coin, PHI_in;
    int m, k, d;

    WK_state(); //Generate the WK state.

    for(k = 0; k < K; ++k) {
        phi = uniform_rng->next_rand(0.0, 1.0) * (upper - lower) + lower;
        PHI = 0;
        //copy input state: the optimal solution across all compilers is memcpy:
        //nadeausoftware.com/articles/2012/05/c_c_tip_how_copy_memory_quickly
        memcpy(state, input_state, (num + 1)*sizeof(dcmplx));
        //Begining the measurement of one sample
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
                //PHI_in = rand_Gaussian(PHI, THETA_DEV); //select if the noise is normally distributed.
		//PHI_in = rand_Hat(PHI, THETA_DEV); //select if the noise is distributed as a Hat function.
		//PHI_in = Lognormal(MU,THETA,PHI); //select if the noise is distributed as a lognormal function.
		//PHI_in = rand_RTN(PHI,Ps,THETA_DEV);//select if the noise is random telegraph.
		PHI_in = rand_skewed(PHI, THETA_DEV, RATIO); //select if the noise is skewed normal.
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
        //store fitness values
        sharp.real(sharp.real() + cos(phi - PHI));
        sharp.imag(sharp.imag() + sin(phi - PHI));
        }
    //find the averages and return
    fitarray[0] = abs(sharp) / double(K); //Calculate the dispersion
    fitarray[1] = atan2(sharp.imag(),sharp.real()); //Calculate the mean direction
    }

bool Phase::T_condition(double *fitarray, int *numvar, int N_cut, bool *mem_ptype, double *memory_forT) {
    /*This function contains the conditions that has to be checked after a step T elapses
    * before the algorithm decides to accept or reject the solution.
    * In particular this function is called when time step is used as the main condition to end the optimization.
    * For this particular problem, it resets the number of variables so the optimization starts over.
    */
	//Let's start this condition from scratch.
    bool type;

    //The conditions are checked only if the algorithm is going to change the way it initializes the population.
    //Policy type 1 is the one with pi bias and is susceptible to loss. CThe algorithm will run until type 0 is found.
    if(*numvar == N_cut - 1) {
        try {
            type = check_policy(fitarray[1], fitarray[0]);
            }
        catch(invalid_argument) {
            fitarray[0] = 0.999999;
            type = check_policy(fitarray[1], fitarray[0]);
            }
        mem_ptype[0] = type;
        if(type == 1) {
            *numvar = *numvar - 1;
            }
        else if(type == 0) {
            memory_forT[0] = fitarray[0];
            memory_forT[1] = fitarray[1];
            }
        else {}
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
//        if(mem_ptype[0] | type) {
        if(check_policy(memory_forT[1], memory_forT[0]) | type) {
            //the policy is bad
            //reset the policy found in numvar=N_cut-1
            *numvar = N_cut - 2;
            }
        }
    return 1;
    }

bool Phase::error_condition(double *current_fitarray, double *memory_fitarray, int data_size, double goal) {
    /*This function contains the function to compute error for when the algorithm is set to use error as the accept-reject condition.
    *It allows for the information of previous optimization to be used to compute the condition as well as using the latest result.
    *In this particular problem, we use the error that corresponds to the confidence interval of 0.98
    * from the linear relation between logN and logV_H.
    */

    double slope, intercept;
    double mean_x;
    double tn2;
    double error, error_goal;
    double x[data_size + 1];
    double y[data_size + 1];
    bool out;

    memory_fitarray[2 * (data_size)] = log10(num);
    memory_fitarray[2 * (data_size) + 1] = log10(pow(current_fitarray[0], -2) - 1);

    //split into x-y arrays

    for(int i = 0; i <= data_size; ++i) {
        x[i] = memory_fitarray[2 * i];
        y[i] = memory_fitarray[2 * i + 1];
        }
	
    //Computing the linear equation using previous data
    linear_fit(data_size, x, y, &slope, &intercept, &mean_x);
    //Compute the goal for error from confidence interval of 0.98
    error_goal = error_interval(x, y, mean_x, data_size + 1, slope, intercept);

    goal = (goal + 1) / 2;
    tn2 = quantile(goal);
    error_goal = error_goal * tn2;

    //Compute the distance between the data and the prediction from linear equation
    error = y[data_size]-x[data_size] * slope - intercept;

    //Check if error is smaller than the goal
    if(error <= error_goal) {
        out = 1;
        }
    else {
        out = 0;
        }
    return out;
    }



void Phase::boundary(double *can1) {
    /*This function cuts the value of the variables so they are within the boundary of the problem.
    *This can be replaced by one or more method, such as periodic boundary or normalization.
    */
    for(int i = 0; i < num; ++i) {
        if(can1[i] < lower_bound[i]) {
            can1[i] = lower_bound[i];
            }
        else if (can1[i] > upper_bound[i]) {
            can1[i] = upper_bound[i];
            }
        }
    }

/*private functions: specific to the simulation and not used by the optimization algorithm*/

/*########### Generating Input State ###########*/
/*Generation function*/
void Phase::WK_state() {
    /* This is the main function for generating the input state.
    * We call it here the Wiseman-Killip state, but the more general term would be the sine state.
    * This algorithm for generating WK state is the most precise one to date (measured by how much the normalization differs from zero),
    * but the rounding error still appears as early as N=80 and is detrimental to the shape of the state for N>100.
    * This is a limitation caused by the use of double precision float.
    * For N>100, either switch to quadruple precision, which will cause a large overhead
    * (about 15 for the simulation, but even more for this function), or another method has to be used to approximate the state.
    */
    const double beta = M_PI / 2;
    const double cosN = pow(cos(beta / 2), num);
    tan_beta = tan(beta / 2);

    int n, k;
    double n_part[num + 1];
    double k_part[num + 1];
    double s_part;
    double temp;

    //Preparing constant arrays to save time.

    //initializing array of factorial numbers
    sqrtfac(sqrtfac_mat);
    one_over_fac(overfac_mat);

    /*factors in d_matrix calculated*/
    //calculating n_parts and k_parts in
    for (n = 0; n <= num; ++n) {
        n_part[n] = pow(-1.0, n) * sqrtfac_mat[n] * sqrtfac_mat[num - n] * cosN * pow(tan_beta, n);
        }
    for (k = 0; k <= num; ++k) {
        k_part[k] = 1 / pow(-1.0, k) * sqrtfac_mat[k] * sqrtfac_mat[num - k] * pow(1 / tan_beta, k);
        }

    //Compute the state
    for (n = 0; n <= num; ++n) { //we have N+1 state b/c we include n=0 and n=N.
        temp = 0;
        input_state[n].real(0);
        input_state[n].imag(0);
        for (k = 0; k <= num; ++k) {
            s_part = cal_spart(n, k, num);
            temp = s_part * k_part[k] * sin((k + 1) * M_PI / (num + 2));
            input_state[n].real(input_state[n].real() + temp * cos(M_PI / 2.0 * (k - n)));
            input_state[n].imag(input_state[n].imag() + temp * sin(M_PI / 2.0 * (k - n)));
            }//end k
        input_state[n].real(input_state[n].real() * n_part[n] / sqrt(num / 2.0 + 1));
        input_state[n].imag(input_state[n].imag() * n_part[n] / sqrt(num / 2.0 + 1));
        }//end n
    }

inline double Phase::cal_spart(const int n, const int k, const int N) {
    /*This function calculates the sum in an element of a Wigner d-matrix
    */

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

inline void Phase::one_over_fac(double *over_mat) {
    /*This function calculates one over factorial matrix.*/
    over_mat[0] = 1;
    for (int i = 1; i <= num; ++i) {
        over_mat[i] = over_mat[i - 1] / i;
        }//end i
    }

inline void Phase::sqrtfac(double *fac_mat) {
    /*This function calculates sqrare root of factorial matrix.*/
    //check array size
    fac_mat[0] = 1;
    for (int i = 1; i <= num; ++i) {
        fac_mat[i] = fac_mat[i - 1] * sqrt_cache[i];
        }//end i
    }

/*########### Measurement Functions ###########*/
/*The following functions are involved in the simulation of the noisy interferometer.*/

inline bool Phase::noise_outcome(const double phi, const double PHI, const int N) {
    /*This function computes the output path of a photon from a noisy interferometer
    by computing the probablity of photon coming out of either path.
    This simulation allows for noise in the unitary operation, but we only consider the noise in the phase shift.
    */
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
    /*This function updates the state when one of the photon is loss.*/
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
    /*This function compute the modulo of the phase.
    */
    while(PHI >= 2 * M_PI) {
        PHI = PHI - 2 * M_PI;
        }
    while (PHI < 0) {
        PHI = PHI + 2 * M_PI;
        }
    return PHI;
    }

inline bool Phase::check_policy(double error, double sharp) {
    /*This function takes the bias and output the policy type.
    *A policy is considered to be type zero if its error falls within the uncertainty of the scheme.
    *This is the desirable type as its estimate no bias and the policy can be used when loss is presence.
    *In the error falls outside the uncertainty,
    *it is very likely that the estimate has a pi bias and is the type of policy that fails when there is loss.
    */
    if (sharp == 1.0) {
        throw invalid_argument("sharpness cannot be one.");
        }
    double sd = sqrt(1 / (sharp * sharp) - 1);
    if(fabs(error)>=fabs(M_PI-sd)){
        return 1;
        }
    else {
        return 0;
        }
    }

double Phase::rand_Gaussian(double mean, /*the average theta*/
				double dev /*deviation for distribution*/
				){
	/*creating random number using Box-Muller Method/Transformation*/
	double Z0;//,Z1;
	double U1,U2; /*uniformly distributed random number input*/
	double r;
	
	/*create input between [-1,1]*/
	do{
	U1=2.0*double(rand())/RAND_MAX-1.0;
	U2=2.0*double(rand())/RAND_MAX-1.0;
	r=U1*U1+U2*U2;
	}while(r==0.0||r>=1.0);
	/*using Box-Muller Transformation*/
	Z0=U1*sqrt(-2*log(r)/r);
	
	return Z0*dev+mean;
	}/*end of rand_Gaussian*/

double Phase::rand_Hat(double PHI, double dev){
	return gaussian_rng->next_rand(PHI, dev);
}

inline double Phase::inv_erf(double x) {
    if(x == 1) {
        throw invalid_argument("Input leads to error.");
        }
    double a = 0.140012;
    double lnx = log(1 - x * x);
    double temp = sqrt(pow(2.0 / (M_PI * a) + lnx / 2.0, 2) - lnx / a);

    return sgn(x) * sqrt(temp - 2.0 / (M_PI * a) - lnx / 2.0);
    }

inline int Phase::sgn(double x) {
    /** Sign function of x: https://en.wikipedia.org/wiki/Sign_function
    */
    if(x < 0) {
        return -1;
        }
    else if(x == 0) {
        return 0;
        }
    else {
        return 1;
        }
    }

double Phase::Lognormal(double mu, double sigma, double peak) {
    double mode = exp(mu - sigma * sigma); //where the peak is, so we know how much to move the distribution later.
    double diff = mode - peak;
    double ans;

    double u = (double)(rand()) / ((double)(RAND_MAX));
    double p = 2 * u - 1;

    try{
	ans=exp(sqrt(2) * sigma * inv_erf(p) + mu) - diff;
	}
    catch(invalid_argument){
	u = (double)(rand()) / ((double)(RAND_MAX));
	p = 2 * u - 1;
	ans=exp(sqrt(2) * sigma * inv_erf(p) + mu) - diff;
	}	

    return ans;
    }

double Phase::rand_RTN(double PHI,double ps,double dev){
	double coin = double(rand())/double(RAND_MAX);
	double ans;
		if (coin<=ps) {
		    coin = double(rand())/double(RAND_MAX);
		    if (coin<0.5){
		        ans = PHI + dev; 
			}
		    else{
			ans = PHI - dev; 
			}
		}
		else {
		    ans=PHI;
		}
	return ans;
}

double Phase::rand_skewed(double mean, double dev, double ratio) {
    double u1, u2;
    double k1, k2;

    double alpha = ratio*dev;

    k1 = (1 + alpha) / sqrt(2 * (1 + alpha * alpha));
    k2 = (1 - alpha) / sqrt(2 * (1 + alpha * alpha));

    u1 = rand_Gaussian(0, 1);
    u2 = rand_Gaussian(0, 1);

    return (k1 * max(u1, u2) + k2 * min(u1, u2))*dev+mean;
    }