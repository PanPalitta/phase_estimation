/*NOTES on using fitness functions for phase estimation problem with loss*/
/* Both avg_fitness() and fitness() contain the same code.
 * The policies are learned with out loss (loss in avg_fitness() set to zero).
 * Then the policies are selected based on its average fitness value over
 * the lossy case (loss in fitness() set to other than zero) which is called
 * through avg_Final_select() in OptAlg class.
 * */

#include "phase_loss_opt.h"

Phase::Phase(int numvar) {
    int i;
    lower=0;
    upper=2*M_PI;
    num=numvar;
    num_repeat=10*num*num;
    lower_bound=new double[num];
    upper_bound=new double[num];
    for(i=0; i<num; i++) {
        lower_bound[i]=lower;
        upper_bound[i]=upper;
    }

    state=new dcmplx[num+1];
    update0=new dcmplx[num+1];
    update1=new dcmplx[num+1];
}

Phase::~Phase() {
    delete state;
    delete update0;
    delete update1;
    delete input_state;
}

double Phase::fitness(double *soln) {
    dcmplx sharp=0;
    bool dect_result;
    double PHI;
    double phi;
    int m,k,d;
    int K=10*num*num;

    double coin;
    double PHI_in;
    double loss=0.2;//loss level
    double theta_dev=0.0;//M_PI;//phase noise level

    WK_state();

    for(k=0; k<K; k++) {
        phi=double(rand())/RAND_MAX*(upper-lower)+lower;
        PHI=0;
        /*copy input state*/
        for (m=0; m<=num; m++) {
            state[m]=input_state[m];
        }
        /*measurement*/
        d=0;
        for (m=0; m<num; m++) {
            //randomly decide whether loss occurs
            coin=double(rand())/RAND_MAX;
            if(coin<=loss) {
                state_loss(num-m);//update only the state using loss function
            }
            else {
                PHI_in=rand_Gaussian(PHI,theta_dev);
                PHI_in=mod_2PI(PHI_in);//noisy PHI
                dect_result=noise_outcome(phi,PHI_in,num-m);
                //dect_result=outcome(phi,PHI_in,num-m);
                if(dect_result==0) {
                    PHI=PHI-soln[d];
                }
                else {
                    PHI=PHI+soln[d];
                }
                PHI=mod_2PI(PHI);
                d++;
            }
        }
        sharp.real()+=cos(phi-PHI);
        sharp.imag()+=sin(phi-PHI);
    }
    return abs(sharp)/double(K);
}

double Phase::avg_fitness(double *soln, int K) {
    dcmplx sharp=0;
    bool dect_result;
    double PHI;
    double phi;
    int m,k,d;

    double coin;
    double PHI_in;
    double loss=0.0;//loss level
    double theta_dev=0.0;//M_PI;//phase noise level


    WK_state();

    for(k=0; k<K; k++) {
        phi=double(rand())/RAND_MAX*(upper-lower)+lower;
        PHI=0;
        /*copy input state*/
        for (m=0; m<=num; m++) {
            state[m]=input_state[m];
        }
        /*measurement*/
        d=0;
        for (m=0; m<num; m++) {
            //randomly decide whether loss occurs
            coin=double(rand())/RAND_MAX;
            if(coin<=loss) {
                state_loss(num-m);//update only the state using loss function
            }
            else {
                PHI_in=rand_Gaussian(PHI,theta_dev);
                PHI_in=mod_2PI(PHI_in);//noisy PHI
                dect_result=noise_outcome(phi,PHI_in,num-m);
                //dect_result=outcome(phi,PHI_in,num-m);
                if(dect_result==0) {
                    PHI=PHI-soln[d];
                }
                else {
                    PHI=PHI+soln[d];
                }
                PHI=mod_2PI(PHI);
                d++;
            }
        }
        sharp.real()+=cos(phi-PHI);
        sharp.imag()+=sin(phi-PHI);
    }
    return abs(sharp)/double(K);
}

/*private functions*/
/*########### Generating Input State ###########*/
/*Aux functions*/
void Phase::one_over_fac(double *over_mat) { //calculate one over factorial matrix
    int i;
    for(i=0; i<=num; i++) {
        if(i==0) {
            over_mat[i]=1;
        }
        else {
            over_mat[i]=over_mat[i-1]/i;
        }
    }//end i
}

void Phase::sqrtfac(double *fac_mat) { //calculate sqrt of factorial matrix
    int i;
    for(i=0; i<=num; i++) {
        if(i==0) {
            fac_mat[i]=1;
        }
        else {
            fac_mat[i]=fac_mat[i-1]*sqrt(i);
        }
    }//end i
}

double Phase::cal_spart(int n, int k, int N) { //calculating the Wigner d-matrix element
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
    for (s=s_min; s<=s_max; s++) {
        s_part=s_part+pow(-1.0,s)*overfac_mat[k-s]*overfac_mat[s]*overfac_mat[n-k+s]*overfac_mat[N-n-s]*pow(tan_beta,2*s);
    }
    return s_part;
}


/*Generation function*/
void Phase::WK_state() {
    int n,k;
    double beta=M_PI/2;
    double n_part[num+1];
    double k_part[num+1];
    double s_part;
    double temp;
    //double total=0;

    input_state=new dcmplx[num+1];
    sqrtfac_mat=new double[num+1];
    overfac_mat=new double[num+1];

    cosN=pow(cos(beta/2),num); //calculating cosine to the power of N.
    tan_beta=tan(beta/2); //calculating tangent of beta/2

    sqrtfac(sqrtfac_mat); //initializing
    one_over_fac(overfac_mat); //initializing

    /*factors in d_matrix calculated*/
    //calculating n_parts and k_parts in
    for(n=0; n<=num; n++) {
        n_part[n]=pow(-1.0,n)*sqrtfac_mat[n]*sqrtfac_mat[num-n]*cosN*pow(tan_beta,n);
    }
    for(k=0; k<=num; k++) {
        k_part[k]=1/pow(-1.0,k)*sqrtfac_mat[k]*sqrtfac_mat[num-k]*pow(1/tan_beta,k);
    }

    for(n=0; n<=num; n++) { //we have N+1 state b/c we include n=0 and n=N.
        temp=0;
        input_state[n].real()=0;
        input_state[n].imag()=0;
        for(k=0; k<=num; k++) {
            s_part=cal_spart(n,k,num);
            temp=s_part*k_part[k]*sin((k+1)*M_PI/(num+2));
            input_state[n].real()+=temp*cos(M_PI/2*(k-n));
            input_state[n].imag()+=temp*sin(M_PI/2*(k-n));
        }//end k
        input_state[n].real()=input_state[n].real()*n_part[n]/sqrt(num/2+1);
        input_state[n].imag()=input_state[n].imag()*n_part[n]/sqrt(num/2+1);
    }//end n
    //releasing memories
    delete sqrtfac_mat;
    delete overfac_mat;
}
/*########### Measurement Functions ###########*///Measurement function
bool Phase::outcome(double phi, double PHI, int N) {
//N is the number of photons currently available, not equal to 'num'
    int n;
    bool result=0;
    double theta=(phi-PHI)/2;
    double prob[2]= {0,0};

    double cos_theta=cos(theta)/sqrt(N);
    double sin_theta=sin(theta)/sqrt(N);

    for (n=0; n<N; n++) {
        //if C_0 is measured
        update0[n]=state[n+1]*sin_theta*sqrt(n+1)+state[n]*cos_theta*sqrt(N-n);
        prob[0]+=abs(update0[n]*conj(update0[n]));
        //if C_1 is measured
        update1[n]=state[n+1]*cos_theta*sqrt(n+1)-state[n]*sin_theta*sqrt(N-n);
        prob[1]+=abs(update1[n]*conj(update1[n]));
    }


    if ((double(rand())/RAND_MAX)<=(prob[0]/(prob[0]+prob[1]))) {
        //measurement outcome is 0
        for(n=0; n<N; n++) {
            state[n]=update0[n]/sqrt(prob[0]);
        }
        state[N]=0;
        result=0;
    }
    else { //measurement outcome is 1
        for(n=0; n<N; n++) {
            state[n]=update1[n]/sqrt(prob[1]);
        }
        state[N]=0;
        result=1;
    }
    return result;
}

bool Phase::noise_outcome(double phi, double PHI, int N) {
//N is the number of photons currently available, not equal to 'num'
    int n;
    bool result=0;
    double theta=(phi-PHI)/2;
    double prob[2]= {0,0};
    double oper_n[3];
    dcmplx U[2][2];

    double dev_n=0.0;//level of noise in operator

    double cos_theta=cos(theta)/sqrt(N);
    double sin_theta=sin(theta)/sqrt(N);
    //noise in operator: currently not in use
    oper_n[0]=rand_Gaussian(0,dev_n);//n_x
    oper_n[2]=rand_Gaussian(0,dev_n);//n_z
    oper_n[1]=sqrt(1-(oper_n[0]*oper_n[0]+oper_n[2]*oper_n[2]));

    U[0][0].real()=sin_theta*oper_n[1];
    U[0][0].imag()=-oper_n[0]*sin_theta;
    U[0][1].real()=cos_theta;
    U[0][1].imag()=oper_n[2]*sin_theta;
    U[1][0].real()=cos_theta;
    U[1][0].imag()=-oper_n[2]*sin_theta;
    U[1][1].real()=sin_theta*oper_n[1];
    U[1][1].imag()=sin_theta*oper_n[0];

    for (n=0; n<N; n++) {
        //if C_0 is measured
        update0[n]=state[n+1]*U[0][0]*sqrt(n+1)+state[n]*U[0][1]*sqrt(N-n);
        prob[0]+=abs(update0[n]*conj(update0[n]));
        //if C_1 is measured
        update1[n]=state[n+1]*U[1][0]*sqrt(n+1)-state[n]*U[1][1]*sqrt(N-n);
        prob[1]+=abs(update1[n]*conj(update1[n]));
    }


    if ((double(rand())/RAND_MAX)<=(prob[0]/(prob[0]+prob[1]))) {
        //measurement outcome is 0
        for(n=0; n<N; n++) {
            state[n]=update0[n]/sqrt(prob[0]);
        }
        state[N]=0;
        result=0;
    }
    else { //measurement outcome is 1
        for(n=0; n<N; n++) {
            state[n]=update1[n]/sqrt(prob[1]);
        }
        state[N]=0;
        result=1;
    }
    return result;
}

void Phase::state_loss(int N) {
    //state update if the photon is loss
    int n;
    double factor=1/sqrt(2*N);

    for(n=0; n<N; n++) {
        update0[n]=state[n]*sqrt(N-n)+state[n+1]*sqrt(n+1);
    }
    for(n=0; n<N; n++) {
        state[n]=update0[n]*factor;
    }
    state[N]=0;
}

double Phase::rand_Gaussian(double mean, /*the average theta*/
                            double dev /*deviation for distribution*/
                           ) {
    /*creating random number using Box-Muller Method/Transformation*/
    double Z0;//,Z1;
    double U1,U2; /*uniformly distributed random number input*/
    double r;

    /*create input between [-1,1]*/
    do {
        U1=2.0*double(rand())/RAND_MAX-1.0;
        U2=2.0*double(rand())/RAND_MAX-1.0;
        r=U1*U1+U2*U2;
    } while(r==0.0||r>=1.0);
    /*using Box-Muller Transformation*/
    Z0=U1*sqrt(-2*log(r)/r);

    return Z0*dev+mean;
}/*end of rand_Gaussian*/

double Phase::mod_2PI(double PHI) {
    while(PHI>=2*M_PI) {
        PHI=PHI-2*M_PI;
    }
    while (PHI<0) {
        PHI=PHI+2*M_PI;
    }
    return PHI;
}
