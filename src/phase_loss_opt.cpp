/*NOTES on using fitness functions for phase estimation problem with loss*/
/* Both avg_fitness() and fitness() contain the same code.
 * The policies are learned with out loss (loss in avg_fitness() set to zero).
 * Then the policies are selected based on its average fitness value over
 * the lossy case (loss in fitness() set to other than zero) which is called
 * through avg_Final_select() in OptAlg class.
 * */
#include <iostream>
#include <cstdlib>
#include <cstring>

#include "phase_loss_opt.h"
#ifdef CUDA
#include "rng_gpu.h"
#endif

Phase::Phase(const int numvar) {
    int i;
    lower = 0;
    upper = 2*M_PI;
    num = numvar;
    num_repeat = 10*num*num;
    lower_bound = new double[num];
    upper_bound = new double[num];
    for(i=0; i<num; ++i) {
        lower_bound[i] = lower;
        upper_bound[i] = upper;
    }
    sqrt_cache = new double[num+1];
    for(i=0; i<num+1; ++i) {
        sqrt_cache[i] = sqrt(i);
    }
    input_state = new dcmplx[num+1];
    sqrtfac_mat = new double[num+1];
    overfac_mat = new double[num+1];
    state = new dcmplx[num+1];
    update0 = new dcmplx[num+1];
    update1 = new dcmplx[num+1];
    //Maximum number of uniform random numbers we will use in one go
    n_urandom_numbers = num_repeat+2*num_repeat*num;
    urandom_numbers = new double[n_urandom_numbers];
    index_urandom_numbers = 0;
    //Maximum number of Gaussian random numbers we will use in one go
    n_grandom_numbers = 3*num_repeat*num;
    grandom_numbers = new double[n_grandom_numbers];
    index_grandom_numbers = 0;
#ifdef CUDA
    gpu_cache_alloc(n_urandom_numbers, n_grandom_numbers);
#endif
}

Phase::~Phase() {
    delete[] state;
    delete[] update0;
    delete[] update1;
    delete[] input_state;
    delete[] sqrtfac_mat;
    delete[] overfac_mat;
    delete[] grandom_numbers;
    delete[] urandom_numbers;
#ifdef CUDA
    gpu_cache_free();
#endif
}

double Phase::fitness(double *soln) {
    const double loss = 0.2;//loss level
    const int K = 10*num*num;
    dcmplx sharp(0.0, 0.0);
    bool dect_result;
    double PHI, phi, coin, PHI_in;
    int m, k, d;

    WK_state();
    for(k=0; k<K; ++k) {
        phi = double(rand())/RAND_MAX*(upper-lower)+lower;
        PHI = 0;
        //copy input state: the optimal solution across all compilers is memcpy:
        //nadeausoftware.com/articles/2012/05/c_c_tip_how_copy_memory_quickly
        memcpy(state, input_state, (num+1)*sizeof(dcmplx));
        //measurement
        d = 0;
        for (m=0; m<num; ++m) {
            //randomly decide whether loss occurs
            coin = double(rand())/RAND_MAX;
            if(coin<=loss) {
                state_loss(num-m);//update only the state using loss function
            } else {
                PHI_in = rand_Gaussian(PHI, THETA_DEV);
                PHI_in = mod_2PI(PHI_in);//noisy PHI
                dect_result = noise_outcome(phi, PHI_in, num-m);
                //dect_result=outcome(phi,PHI_in,num-m);
                if(dect_result == 0) {
                    PHI = PHI-soln[d++];
                } else {
                    PHI = PHI+soln[d++];
                }
                PHI = mod_2PI(PHI);
            }
        }
        sharp.real()=sharp.real()+cos(phi-PHI);
        sharp.imag()=sharp.imag()+sin(phi-PHI);
    }
    return abs(sharp)/double(K);
}

double Phase::avg_fitness(double *soln, const int K) {
    const double loss = 0.0;//loss level
    dcmplx sharp(0.0, 0.0);
    bool dect_result;
    double PHI, phi, coin, PHI_in;
    int m, k, d;
    // Random numbers are generated in advance to take advantage of
    // vectorization
    index_urandom_numbers = 0;
    index_grandom_numbers = 0;
#ifndef CUDA
    init_urandom_number_cache(K+2*K*num);
    init_grandom_number_cache(3*K*num);
#else
    gpu_cache_init(urandom_numbers, K+2*K*num, grandom_numbers, 3*K*num);
#endif
    WK_state();
    //cout << "HERE " << K << "\n";
    for(k=0; k<K; ++k) {
        //cout << k << " ";
        phi = next_urand()*(upper-lower)+lower;
        PHI = 0;
        //copy input state: the optimal solution across all compilers is memcpy:
        //nadeausoftware.com/articles/2012/05/c_c_tip_how_copy_memory_quickly
        memcpy(state, input_state, (num+1)*sizeof(dcmplx));
        //measurement
        d = 0;
        for (m=0; m<num; ++m) {
            // This loop is the most critical part of the entire program. It
            // executes K*num=10*num^3 times on each call of avg_fitness. All
            // optimization should focus on this loop.

            //randomly decide whether loss occurs
            coin = next_urand();
            if(coin<=loss) {
                state_loss(num-m);//update only the state using loss function
            } else {
                PHI_in = next_grand(PHI, THETA_DEV);
                PHI_in = mod_2PI(PHI_in);//noisy PHI
                dect_result = noise_outcome(phi, PHI_in, num-m);
                //dect_result = outcome(phi,PHI_in,num-m);
                if (dect_result==0) {
                    PHI = PHI-soln[d++];
                } else {
                    PHI = PHI+soln[d++];
                }
                PHI = mod_2PI(PHI);
            }
        }
        sharp.real()=sharp.real()+cos(phi-PHI);
        sharp.imag()=sharp.imag()+sin(phi-PHI);
    }
  return abs(sharp)/double(K);
}

/*private functions*/
/*########### Generating Input State ###########*/
/*Aux functions*/
inline void Phase::one_over_fac(double *over_mat) { //calculate one over factorial matrix
    int i;
    over_mat[0] = 1;
    for (i=1; i<=num; ++i) {
        over_mat[i] = over_mat[i-1]/i;
    }//end i
}

inline void Phase::sqrtfac(double *fac_mat) { //calculate sqrt of factorial matrix
    int i;
    fac_mat[0] = 1;
    for (i=1; i<=num; ++i) {
        fac_mat[i] = fac_mat[i-1]*sqrt_cache[i];
    }//end i
}

inline double Phase::cal_spart(const int n, const int k, const int N) {
    //calculating the Wigner d-matrix element
    int s;
    int s_min;
    int s_max;
    double s_part=0;

    //find lower limit
    if(n-k>=0) {
        s_min=0;
    }
    else if(n-k<0) {
        s_min=k-n;
    }
    else {}
    //find upper limit
    if(k<=N-n) {
        s_max=k;
    }
    else if(k>N-n) {
        s_max=N-n;
    }
    else {}

    //calculating the s_part
    for (s=s_min; s<=s_max; ++s) {
        s_part=s_part+pow(-1.0,s)*overfac_mat[k-s]*overfac_mat[s]*overfac_mat[n-k+s]*overfac_mat[N-n-s]*pow(tan_beta,2*s);
    }
    return s_part;
}


/*Generation function*/
void Phase::WK_state() {
    const double beta=M_PI/2;
    const double cosN=pow(cos(beta/2), num);
    tan_beta=tan(beta/2);

    int n,k;
    double n_part[num+1];
    double k_part[num+1];
    double s_part;
    double temp;
    //double total=0;

    sqrtfac(sqrtfac_mat); //initializing
    one_over_fac(overfac_mat); //initializing

    /*factors in d_matrix calculated*/
    //calculating n_parts and k_parts in
    for (n=0; n<=num; ++n) {
        n_part[n] = pow(-1.0,n)*sqrtfac_mat[n]*sqrtfac_mat[num-n]*cosN*pow(tan_beta, n);
    }
    for (k=0; k<=num; ++k) {
        k_part[k] = 1/pow(-1.0,k)*sqrtfac_mat[k]*sqrtfac_mat[num-k]*pow(1/tan_beta, k);
    }

    for (n=0; n<=num; ++n) { //we have N+1 state b/c we include n=0 and n=N.
        temp=0;
        input_state[n].real()=0;
        input_state[n].imag()=0;
        for (k=0; k<=num; ++k) {
            s_part = cal_spart(n,k,num);
            temp = s_part*k_part[k]*sin((k+1)*M_PI/(num+2));
            input_state[n].real()=input_state[n].real()+temp*cos(M_PI/2.0*(k-n));
            input_state[n].imag()=input_state[n].imag()+temp*sin(M_PI/2.0*(k-n));
        }//end k
        // FIXED:This correction corrects the probability overshoot.
        input_state[n].real()=input_state[n].real()*n_part[n]/sqrt(num/2.0+1);
        input_state[n].imag()=input_state[n].imag()*n_part[n]/sqrt(num/2.0+1);
    }//end n
}

/*########### Measurement Functions ###########*///Measurement function
inline bool Phase::outcome(const double phi, const double PHI, const int N) {
    //N is the number of photons currently available, not equal to 'num'
    const double theta = (phi-PHI)/2;
    const double cos_theta=cos(theta)/sqrt_cache[N];
    const double sin_theta=sin(theta)/sqrt_cache[N];
    int n;
    double prob0 = 0.0, prob1 = 0.0;
    for (n=0; n<N; ++n) {
        //if C_0 is measured
        update0[n] = state[n+1]*sin_theta*sqrt_cache[n+1]+state[n]*cos_theta*sqrt_cache[N-n];
        prob0 += abs(update0[n]*conj(update0[n]));
        //if C_1 is measured
        //This is cache-optimized: we update state[n] in-place
        //state[n] = state[n+1]*cos_theta*sqrt_cache[n+1]-state[n]*sin_theta*sqrt_cache[N-n];
        //prob1 += abs(state[n]*conj(state[n]));
    }
    //FIXME: Why don't prob0+prob1 always sum to 1?
    //The condition is set this way so that when the rounding error starts to kick in at num=82
    //(the state is valid but the numerical value started to go off slightly)
    //it does not effect the random selection of output path.
    if ((double(rand())/RAND_MAX) <= prob0) {
        //measurement outcome is 0
        state[N] = 0;
        prob0 = 1.0/sqrt(prob0);
        for(n=0; n<N; ++n) {
            state[n] = update0[n] * prob0;
        }
        return 0;
    } else {
        //measurement outcome is 1
        for(n=0; n<N; ++n) {
            //FIXME: Why don't we multiply the RHS with (1.0-prob0) and save the
            //second for loop here? Can't we get rid of prob1 completely?
            //If not, rename prob0 to prob, and reuse it here. It saves a double
            //in a critical region of the code.
            state[n] = state[n+1]*cos_theta*sqrt_cache[n+1]-state[n]*sin_theta*sqrt_cache[N-n];
            prob1 += abs(state[n]*conj(state[n]));
        }
        state[N] = 0;
        prob1 = 1.0/sqrt(prob1);
        for(n=0; n<N; ++n) {
            state[n] *= prob1;
        }
        return 1;
    }
}

inline bool Phase::noise_outcome(const double phi, const double PHI, const int N) {
    //N is the number of photons currently available, not equal to 'num'
    const double theta = (phi-PHI)/2.0;
    const double cos_theta = cos(theta)/sqrt_cache[N];
    const double sin_theta = sin(theta)/sqrt_cache[N];
    //noise in operator: currently not in use
    const double oper_n0 = next_grand(0.0, DEV_N);//n_x
    const double oper_n2 = next_grand(0.0, DEV_N);//n_z
    const double oper_n1 = sqrt(1.0-(oper_n0*oper_n0+oper_n2*oper_n2));
    const dcmplx U00(sin_theta*oper_n1, -oper_n0*sin_theta);
    const dcmplx U01(cos_theta, oper_n2*sin_theta);
    const dcmplx U10(cos_theta, -oper_n2*sin_theta);
    const dcmplx U11(sin_theta*oper_n1, sin_theta*oper_n0);
    int n;
    double prob0 = 0.0, prob1 = 0.0;
    for (n=0; n<N; ++n) {
        //if C_0 is measured
        update0[n] = state[n+1]*U00*sqrt_cache[n+1]+state[n]*U01*sqrt_cache[N-n];
        prob0 += abs(update0[n]*conj(update0[n]));
        //if C_1 is measured
        //This is cache-optimized: we update state[n] in-place
        //state[n] = state[n+1]*U10*sqrt_cache[n+1]-state[n]*U11*sqrt_cache[N-n];
        //prob1 += abs(state[n]*conj(state[n]));
    }

    if (next_urand() <= prob0) {
        //measurement outcome is 0
        state[N] = 0;
        prob0 = 1.0/sqrt(prob0);
        for(n=0; n<N; ++n) {
            state[n] = update0[n] * prob0;
        }
        return 0;
    } else {
        //measurement outcome is 1
        for(n=0; n<N; ++n) {
            state[n] = state[n+1]*U10*sqrt_cache[n+1]-state[n]*U11*sqrt_cache[N-n];
            prob1 += abs(state[n]*conj(state[n]));
        }
        state[N] = 0;
        prob1 = 1.0/sqrt(prob1);
        for(n=0; n<N; ++n) {
            state[n] *= prob1;
        }
        return 1;
    }
}

inline void Phase::state_loss(const int N) {
    //state update if the photon is loss
    double factor=1/sqrt(2*N);
    for(int n=0; n<N; ++n) {
        state[n] = (state[n]*sqrt_cache[N-n]+state[n+1]*sqrt_cache[n+1])*factor;
    }
    state[N] = 0;
}

inline double Phase::mod_2PI(double PHI) {
    while(PHI>=2*M_PI) {
        PHI=PHI-2*M_PI;
    }
    while (PHI<0) {
        PHI=PHI+2*M_PI;
    }
    return PHI;
}

/*########### RNG Functions ###########*/
inline double Phase::rand_Gaussian(const double mean, /*the average theta*/
                                   const double dev /*deviation for distribution*/
                                  ) {
    /*creating random number using Box-Muller Method/Transformation*/
    double U1,U2; /*uniformly distributed random number input*/
    double r;

    /*create input between [-1,1]*/
    do {
        U1=2.0*double(rand())/RAND_MAX-1.0;
        U2=2.0*double(rand())/RAND_MAX-1.0;
        r=U1*U1+U2*U2;
    } while(r==0.0||r>=1.0);
    /*using Box-Muller Transformation*/

    return U1*sqrt(-2*log(r)/r)*dev+mean;
}/*end of rand_Gaussian*/

void Phase::init_urandom_number_cache(const int n) {
#ifndef CUDA
    for (int i=0; i<n; ++i) {
        // Comment out this line if not using vectorized version
        // urandom_numbers[i] = double(rand())/RAND_MAX;
    }
#endif
    index_urandom_numbers = 0;
}

inline double Phase::next_urand() {
    // Warning: for speed, there is no bound checking! See
    // nadeausoftware.com/articles/2013/12/
    // c_tip_considering_use_exceptions_asserts_and_bounds_checking

#ifndef CUDA
    // Flip these lines to change between vectorized and non-vectorized variants
    return double(rand())/RAND_MAX;
    //return urandom_numbers[index_urandom_numbers++];
#else
    return urandom_numbers[index_urandom_numbers++];
#endif
}

void Phase::init_grandom_number_cache(const int n) {
    for (int i=0; i<n; ++i) {
        // Comment out this line if not using vectorized version
        // grandom_numbers[i] = rand_Gaussian(0.0, 1.0);
    }
    index_grandom_numbers = 0;
}

inline double Phase::next_grand(const double mean, const double dev) {
    // Warning: for speed, there is no bound checking! See
    // nadeausoftware.com/articles/2013/12/
    // c_tip_considering_use_exceptions_asserts_and_bounds_checking

#ifndef CUDA
    // Flip these lines to change between vectorized and non-vectorized variants
    return rand_Gaussian(mean, dev);
    //return grandom_numbers[index_grandom_numbers++]*dev+mean;
#else
    return grandom_numbers[index_grandom_numbers++]*dev+mean;
#endif
}

void Phase::status_rand() {
    cout << index_grandom_numbers << "/" << n_grandom_numbers << endl;
}
