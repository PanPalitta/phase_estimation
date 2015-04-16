#include <cstdlib>
#include <iostream>
#include <typeinfo>
#include <mpi.h>

#include "phase_loss_opt.h"
#include "mpi_de.h"
#include "subfunc.h"
#include "mpi_pso.h"
#ifdef CUDA
#include "rng_gpu.h"
#endif

using namespace std;

typedef complex<double> dcmplx;

int main(int argc, char **argv) {

    /*mpi handlers*/
    int my_rank;
    int nb_proc;
    time_t start_time;
    MPI_Status status;
    MPI_Datatype MPI_TYPE;
    int tag=1;

    /*variables*/
    int numvar;
    double *solution;//the type of this array must correspond to that of the solution of the problem.
    int p,t;
    double final_fit;
    double *soln_fit;
    int *can_per_proc;
    double *x;
    double *y;

    /*parameter settings*/
    int pop_size=12;
    int N_begin=4;
    int N_cut=7;
    int N_end=15;
    int iter=100;
    int iter_begin=300;
    int T_cut_off=N_cut;
    double prev_dev=0.01*M_PI;
    double new_dev=0.25*M_PI;
    int repeat=10;
	
	int data_start=N_begin;
    int data_end=94;
    double t_goal=0.98; //probability for calculating quantile
    int data_size=data_end-data_start;
    double slope=0.0, intercept=0.0;
	double mean_x=0.0, SSres=0.0;
	double TSSres,Tmean_x;
	double fit_goal,error;

    if(N_cut<N_begin) {
        cout<<"please select new N_cut>"<<N_begin<<":";
        cin >> N_cut;
    }
    else {}

    /*start mpi*/
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nb_proc);
#ifdef CUDA
    setDevice(my_rank, nb_proc);
#endif
    srand(time(NULL)+my_rank);

    soln_fit=new double[pop_size];//create an array to store global fitness from each candidate.
    solution=new double[N_begin];

    x=new double[N_end-data_start];
    y=new double[N_end-data_start];
    //calculating number of candidate per processor and stores the number in an array.
    can_per_proc=new int[nb_proc];
    for(p=0; p<nb_proc; ++p) {
        can_per_proc[p]=0;   //make sure the array started with zeroes.
    }
    for(p=0; p<pop_size; ++p) {
        can_per_proc[p%nb_proc]+=1;
    }


    if(my_rank==0) {
        output_header();   //write header of result files
    }
    else {}

    //solution=new double[N_end];

    for(numvar=N_begin; numvar<=N_end; ++numvar) {
        //cout<<numvar<<":";
        t=0;
		x[numvar-data_start]=log10(numvar);//collect x data
		
        Problem<double>* prob_ptr=new Phase(numvar);
        OptAlg<double>* opt= new DE<double>(prob_ptr);

        start_time=time(NULL);

        if(numvar<N_cut) {
            opt->Init_population(can_per_proc[my_rank]);
        }
        //create candidates according to how many should be on that particular processor.
        else {
            //previous solution must be send from root processor.
            //set MPI_type ***NEED A GENERAL SOLUTION FOR THIS TASK***
            if(typeid(solution[0])==typeid(double)) {
                MPI_TYPE=MPI_DOUBLE;
            }
            else if(typeid(solution[0])==typeid(dcmplx)) {
                MPI_TYPE=MPI_DOUBLE_COMPLEX;
            }
            else {}

            if(my_rank==0) {
                for(p=1; p<nb_proc; ++p) {
                    MPI_Send(&solution[0],numvar,MPI_TYPE,p,tag,MPI_COMM_WORLD);
                }
            }
            else {
                MPI_Recv(&solution[0],numvar,MPI_TYPE,0,tag,MPI_COMM_WORLD,&status);
            }

            opt->Init_previous(prev_dev,new_dev,can_per_proc[my_rank],solution);
            //each processor initialize the candidates.
        }

        opt->put_to_best(my_rank,pop_size,nb_proc);

        delete [] solution;
        solution=new double[numvar+1];//This initialization must happen after the initialization of candidates.

        //setting the success criterion
        if(numvar<T_cut_off) {
            opt->set_success(iter_begin,0);
        }
        else if(numvar>=data_end) {
            opt->set_success(1000,1);   //stop after perform 1000 iteration or exceeding the goal.
        }
        else {
            opt->set_success(iter,0);
        }

        do {
            opt->update_popfit();
            opt->combination(my_rank,pop_size,nb_proc);//this is where a lot of comm goes on between processors
            opt->selection(my_rank,pop_size,nb_proc);
            ++t;

            //root check for success
            final_fit=opt->Final_select(my_rank,pop_size,nb_proc,soln_fit,solution);//again, communicate to find the best solution that exist so far

			if(numvar>=data_end){
			    y[numvar-data_start]=log10(pow(final_fit,-2)-1);
				TSSres=SSres;
				Tmean_x=mean_x;
				//error=opt->error_update(data_size,&TSSres,&Tmean_x,slope,intercept,y,x);
				//cout<<numvar<<":error="<<error<<","<<abs(y[numvar-data_start]-slope*x[numvar-data_start]-intercept)<<endl;
				fit_goal=1/sqrt(pow(10.0,intercept)*pow(double(numvar),slope)+1);				
			}	


            opt->success = opt->check_success(t, numvar, final_fit, slope,
                                              intercept);

            //cout<<t<<":"<<final_fit<<endl;

        } while(opt->success==0);

        final_fit=opt->avg_Final_select(solution,repeat,my_rank,pop_size,nb_proc,soln_fit);

        if(my_rank==0) {
            output_result(numvar,final_fit,solution,start_time);
        }
        else {}
        //if(my_rank==0){cout<<t<<":"<<final_fit<<endl;}else{}
		
		//collect data for linear fit
        if(numvar>=data_start&&numvar<data_end) {
            y[numvar-data_start]=log10(pow(final_fit,-2)-1);
        }
        else {}

        if(numvar==data_end-1) {
            opt->linear_fit(data_size,x,y,&slope,&intercept,&mean_x);
			t_goal=(t_goal+1)/2;
			t_goal=opt->quantile(t_goal);
			error=opt->error_interval(x,y,mean_x,data_size,&SSres,slope,intercept);
			error=error*t_goal;
			cout<<slope<<","<<intercept<<endl;
			//cout<<numvar<<":error="<<error<<endl;
        }
        else if (numvar>=data_end){
			SSres=TSSres;
			mean_x=Tmean_x;
			++data_size;
		}
		else {}

	//testing loss
	

	if(my_rank==0){
		final_fit=prob_ptr->fitness(solution);
		cout<<numvar<<"\t"<<final_fit<<endl;
	}

        opt->~OptAlg();
        prob_ptr->~Problem();
    }

    delete [] solution;
    delete [] soln_fit;
    delete [] x;
    delete [] y;
    delete [] can_per_proc;

    MPI_Finalize();
    cout<<"done"<<endl;

    return 0;
}

