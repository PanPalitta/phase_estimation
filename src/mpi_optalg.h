#ifndef MPIOPTALG_H
#define MPIOPTALG_H

#include<complex>
#include<ctime>
#include<cstdlib>
#include<typeinfo>
#include "mpi.h"
#include "problem.h"
#include "candidate.h"

using namespace std;
using namespace MPI;
typedef complex<double> dcmplx;

template<typename typeT>
class OptAlg{
public:
	OptAlg(){};
	OptAlg(Problem<typeT> *problem_ptr){prob=problem_ptr;num=prob->num;}
	~OptAlg(){};

	virtual void put_to_best(int my_rank, int total_pop, int nb_proc){};
	virtual void combination(int my_rank, int total_pop, int nb_proc){};
	virtual void selection(int my_rank, int total_pop, int nb_proc){};
	virtual void write_param(double *param_array){};
	virtual void fit_to_global(){};
	virtual void find_global(int my_rank, int total_pop, int nb_proc){};

	Problem<typeT>* prob;

	//functions and variables that are not algorithm specific
	void Init_population(int psize);
	void Init_previous(double prev_dev, double new_dev,int psize, typeT *prev_soln);
	void Cont_fitness(int p);
	void Best_fitness(int p);
	void update_popfit();
	double Final_select(int my_rank, int total_pop,int nb_proc, double *fit,typeT *solution);
	double avg_Final_select(typeT* solution,int repeat, int my_rank, int total_pop, int nb_proc, double *soln_fit);

	void set_success(int iter,bool goal);
	bool check_success(int t, int D, double fit,double slope, double intercept);
	void linear_fit(int data_size,double *x, double *y, double *slope, double *intercept);

	double rand_Gaussian(double mean, double dev);
	void dev_gen(double *dev_array,double prev_dev, double new_dev, int cut_off);

    bool success;

protected:
	int num,pop_size,T,t;
	Candidate<typeT> *pop;

	bool goal;
};


template<typename typeT>
void OptAlg<typeT>::Init_population(int psize){
	int i,p;
	int sum_im=0;
	dcmplx temp[num];
	typeT input[num];
	//store the variables
	pop_size=psize;
	pop=new Candidate<typeT>[pop_size];

    srand(time(NULL));

	for(p=0;p<pop_size;p++){
		//generate candidate
		for(i=0;i<num;i++){
			input[i]=double(rand())/RAND_MAX*(prob->upper_bound[i]-prob->lower_bound[i])+prob->lower_bound[i];
			}
		//store it in candidate object
		pop[p].init_can(num);
		pop[p].update_cont(input);
		Cont_fitness(p);
		}

	}

template<typename typeT>
void OptAlg<typeT>::Init_previous(double prev_dev, double new_dev,int psize, typeT *prev_soln){
	int i,p;
	int sum_im=0;
	int dev_cut_off=num-1;
	dcmplx temp[num];
	typeT input[num];
	double dev[num];
	//store the variables
	pop_size=psize;
	pop=new Candidate<typeT>[pop_size];

	prev_soln[num-1]=prev_soln[num-2];
	dev_gen(dev,prev_dev,new_dev,dev_cut_off);

    srand(time(NULL));

	for(p=0;p<pop_size;p++){
		//generate candidate

		for(i=0;i<num;i++){
			temp[i].real()=rand_Gaussian(dcmplx(prev_soln[i]).real(),dev[i]);
			if(dcmplx(prev_soln[i]).imag()==0){temp[i].imag()=0;}
			else{temp[i].imag()=rand_Gaussian(dcmplx(prev_soln[i]).imag(),dev[i]);}

			if(temp[i].imag()==0){}
			else{sum_im+=1;}
			}


		if(sum_im==0){
			for(i=0;i<num;i++){input[i]=abs(temp[i]);}
			}
		else if(sum_im!=0){//TEST THIS FOR COMPLEX VAR LATER!!!!!!
			for(i=0;i<num;i++){
         //           input[i]=temp[i];
				static_cast<dcmplx>(input[i])=temp[i];//need this line if the typeT is not complex
				//cout<<i<<"="<<input[i]<<endl;
				}
			}
		else{}
		//prob->normalize(input);
		prob->modulo(input);
		//store it in candidate object
		pop[p].init_can(num);
		pop[p].update_cont(input);
		Cont_fitness(p);
		}

	}

template<typename typeT>
void OptAlg<typeT>::Cont_fitness(int p){
	double fit;
	typeT soln[num];

		pop[p].read_cont(soln);
		fit=prob->avg_fitness(soln,prob->num_repeat);
		fit+=prob->avg_fitness(soln,prob->num_repeat);
		pop[p].write_contfit(fit,2);

	}

template<typename typeT>
void OptAlg<typeT>::Best_fitness(int p){
	double fit;
	typeT soln[num];

		pop[p].read_best(soln);
		fit=prob->avg_fitness(soln,prob->num_repeat);
		pop[p].write_bestfit(fit);
	}

template<typename typeT>
void OptAlg<typeT>::update_popfit(){
	int p;
	for(p=0;p<pop_size;p++){
		Best_fitness(p);
		}
}
/*##############################Success Criteria#################################*/

template<typename typeT>
void OptAlg<typeT>::set_success(int iter,bool goal_in){
    success=0;
	T=iter;
    goal=goal_in;
	}

template<typename typeT>
bool OptAlg<typeT>::check_success(int t, int D, double fit,double slope, double intercept){
    double fit_goal;
    double fit_del=0;

    if(goal==0){
        if(t>=T){return 1;}
        else{return 0;}
        }
    else{
        fit_goal=1/sqrt(pow(10.0,intercept)*pow(double(D),slope+fit_del/2)+1);
        if(t>=T){return 1;}
        else{
            if(fit>=fit_goal){return 1;}
            else{return 0;}
            }
        }
	}

template<typename typeT>
void OptAlg<typeT>::linear_fit(int data_size,double *x, double *y, double *slope, double *intercept){
    int v;
    double sum_x=0;
    double sum_xx=0;
    double sum_y=0;
    double sum_xy=0;
    double mean_x=0;

    double res_y,res_x;
    double fit_del;

            for(v=0;v<data_size;v++){
                sum_x=sum_x+x[v];
                sum_xx=sum_xx+x[v]*x[v];
                sum_y=sum_y+y[v];
                sum_xy=sum_xy+x[v]*y[v];
                mean_x=mean_x+x[v];
                }

            mean_x=mean_x/double(data_size);

            *slope=(sum_xy-sum_y*sum_x/double(data_size))/(sum_xx-sum_x*sum_x/double(data_size));
            *intercept=sum_y/double(data_size)-*slope*sum_x/double(data_size);

                /*error in slope*/
                /*res_y=0;
                res_x=0;
                for(v=0;v<data_size;v++){
                res_y=res_y+pow(y[v]-intercept-slope*x[v],2);
                res_x=res_x+pow((x[v]-mean_x),2);
                }
                fit_del=sqrt(res_y/res_x/double(data_size-2));

                //cout<<"linear fit: m="<<slope<<",b="<<intercept<<endl;
                /*extrapolation from two data points*/
                /*m=(y[numvar-1]-y[numvar])/(x[numvar-1]-x[numvar]);
                b=y[numvar]-m*x[numvar];*/
                //printf("m=%lf, fit_del=%lf, b=%lf\n",m,fit_del,b);

}



/*##############################Final Selections#################################*/
template<typename typeT>
double OptAlg<typeT>::Final_select(int my_rank, int total_pop,int nb_proc, double *fit,typeT *solution){
	double *soln;
	int p;
	int indx;
	MPI_Status status;
	MPI_Datatype MPI_TYPE;
	int tag=1;
	double global_fit;

	if(typeid(solution[0])==typeid(double)){MPI_TYPE=MPI_DOUBLE;}
	else if(typeid(solution[0])==typeid(dcmplx)){MPI_TYPE=MPI_DOUBLE_COMPLEX;}
	else{}

	fit_to_global();//ensuring that global_best contains the solutions

	//roots needed all the global fitness
	for(p=1;p<total_pop;p++){
		if(my_rank==p%nb_proc){
			global_fit=pop[int(p/nb_proc)].read_globalfit();
			MPI_Send(&global_fit,1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
        }
		else if(my_rank==0){MPI_Recv(&fit[p],1,MPI_DOUBLE,p%nb_proc,tag,MPI_COMM_WORLD,&status);}
		else{}
		}

	MPI_Barrier(MPI_COMM_WORLD);

	//find the candidate that is the solution and send the index to all processor.
	if(my_rank==0){
		soln=&fit[0];
		indx=0;
		for(p=1;p<total_pop;p++){
			if(*soln<fit[p]){soln=&fit[p];indx=p;}
			else{}
			}
		for(p=1;p<nb_proc;p++){MPI_Send(&indx,1,MPI_INT,p%nb_proc,tag,MPI_COMM_WORLD);}
	}
	else if(my_rank!=0){MPI_Recv(&indx,1,MPI_INT,0,tag,MPI_COMM_WORLD,&status);}
	else{}

	MPI_Barrier(MPI_COMM_WORLD);

	//send the best fitness function to all processor
	if(my_rank==0){
		global_fit=*soln;
		for(p=1;p<nb_proc;p++){
		MPI_Send(&global_fit,1,MPI_DOUBLE,p,tag,MPI_COMM_WORLD);
		}
		}
	else if(my_rank!=0){
	MPI_Recv(&global_fit,1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD,&status);
	}
	else{}

	MPI_Barrier(MPI_COMM_WORLD);

	//sending the solution back to root //need to check data type
	if(my_rank==indx%nb_proc){
		pop[int(indx/nb_proc)].read_global(solution);
		MPI_Send(&solution[0],num,MPI_TYPE,0,tag,MPI_COMM_WORLD);
		}
	else if(my_rank==0){MPI_Recv(&solution[0],num,MPI_TYPE,indx%nb_proc,tag,MPI_COMM_WORLD,&status);}
	else{}

	MPI_Barrier(MPI_COMM_WORLD);

	return global_fit;
	}

template<typename typeT>
double OptAlg<typeT>::avg_Final_select(typeT* solution,int repeat, int my_rank, int total_pop, int nb_proc, double *soln_fit){
	MPI_Status status;
	MPI_Datatype MPI_TYPE;
	int tag=1;
	double *soln;
	double final_fit;
	int p,i,indx;
	typeT array[num];
	double fit[pop_size];

	if(typeid(array[0])==typeid(double)){MPI_TYPE=MPI_DOUBLE;}
	else if(typeid(array[0])==typeid(dcmplx)){MPI_TYPE=MPI_DOUBLE_COMPLEX;}
	else{}

	fit_to_global();//move solution to global_best array in case we're using DE.

	//Need to calculate fitness again for 'repeat' times, independently on each
	for(p=0;p<pop_size;p++){
		pop[p].read_global(array);
		//cout<<"previous global_fit ="<<pop[p].read_globalfit()<<",";
		fit[p]=0;
		for(i=0;i<repeat;i++){fit[p]+=prob->avg_fitness(array,prob->num_repeat);}
		fit[p]=fit[p]/repeat;
		//cout<<"new_fit="<<fit[p]<<endl;
		//cout<<my_rank<<","<<p<<":"<<fit[p]<<endl;
		}

	//filling the fitness table in root
	for(p=0;p<total_pop;p++){
        if(p%nb_proc==0){
            soln_fit[p]=fit[p/nb_proc];//no need for transmitting data
        }
        else{
            if(my_rank==p%nb_proc){
                //cout<<"my_rank="<<my_rank<<",p="<<p<<":"<<fit[p/nb_proc]<<endl;
                MPI_Send(&fit[p/nb_proc],1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
                }
            else if(my_rank==0){
                MPI_Recv(&soln_fit[p],1,MPI_DOUBLE,p%nb_proc,tag,MPI_COMM_WORLD,&status);
                //cout<<"my_rank="<<my_rank<<":"<<soln_fit[p]<<endl;
                }
            else{}
            }

		}

    MPI_Barrier(MPI_COMM_WORLD);

	if(my_rank==0){
		soln=&soln_fit[0];
		indx=0;
		//cout<<"0:"<<soln_fit[0]<<endl;
		for(p=1;p<total_pop;p++){
            //cout<<p<<":"<<soln_fit[p]<<endl;
			if(*soln<soln_fit[p]){soln=&soln_fit[p];indx=p;}
			else{}
			}
        final_fit=*soln;
        for(p=1;p<nb_proc;p++){MPI_Send(&final_fit,1,MPI_DOUBLE,p,tag,MPI_COMM_WORLD);}
	}
	else{MPI_Recv(&final_fit,1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD,&status);}

	MPI_Barrier(MPI_COMM_WORLD);

    if(my_rank==0){
        for(p=1;p<nb_proc;p++){MPI_Send(&indx,1,MPI_INT,p,tag,MPI_COMM_WORLD);}
        }
    else{MPI_Recv(&indx,1,MPI_INT,0,tag,MPI_COMM_WORLD,&status);}

    MPI_Barrier(MPI_COMM_WORLD);

	//get solution from the processor
	if(my_rank==indx%nb_proc){
		pop[indx/nb_proc].read_global(array);
		MPI_Send(&array[0],num,MPI_TYPE,0,tag,MPI_COMM_WORLD);
		}
	else if(my_rank==0){
		MPI_Recv(&array[0],num,MPI_TYPE,indx%nb_proc,tag,MPI_COMM_WORLD,&status);
		}
	else{}

	MPI_Barrier(MPI_COMM_WORLD);

    if(my_rank==0){for(i=0;i<num;i++){solution[i]=array[i];}}

	return final_fit;
	}

/*####################Auxilliary functions####################*/

template<typename typeT>
double OptAlg<typeT>::rand_Gaussian(double mean, /*the average theta*/
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

template<typename typeT>
void OptAlg<typeT>::dev_gen(double *dev_array,double prev_dev, double new_dev, int cut_off){
	int i;
	for(i=0;i<num;i++){
		if(i<cut_off){dev_array[i]=prev_dev;}
		else{dev_array[i]=new_dev;}
		}
	}
#endif // OPTALG_H
