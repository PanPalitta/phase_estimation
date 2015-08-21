#ifndef DE_H
#define DE_H

#include "mpi_optalg.h"

template <typename typeT>
class DE : public OptAlg<typeT>
{
public:
	DE() {};
	DE(Problem<typeT> *problem_ptr):F(0.1),Cr(0.6) {
		this->prob=problem_ptr;
		this->num=this->prob->num;
	}
	~DE() {};

	void put_to_best(int my_rank, int total_pop, int nb_proc);
	void combination(int my_rank, int total_pop, int nb_proc);
	void selection(int my_rank, int total_pop, int nb_proc);
	void write_param(double *param_array);
	void fit_to_global();
	void find_global(int my_rank, int total_pop, int nb_proc) {};

private:
	double F, Cr;
};

template <typename typeT>
void DE<typeT>::put_to_best(int my_rank, int total_pop, int nb_proc)
{
	int p;
	for(p=0; p<this->pop_size; ++p) { //pop_size here is the number of candidates assigned for a processor from the initialization.
		this->pop[p].update_best();
	}
}

template <typename typeT>
void DE<typeT>::write_param(double *param_array)
{
	F=param_array[1];
	Cr=param_array[2];
}

template <typename typeT>
void DE<typeT>::combination(int my_rank, int total_pop, int nb_proc)
{
	MPI_Status status;
	MPI_Datatype MPI_TYPE;
	int tag=1;
	int i,f,p;
	double coin;
	int fam_size=3;
	typeT all_soln[total_pop][this->num];
	int fam[fam_size];
	typeT input[this->num];
	typeT can[this->num];

	srand(0+my_rank);

	if(typeid(input[0])==typeid(double)) {
		MPI_TYPE=MPI_DOUBLE;
	} else if(typeid(input[0])==typeid(dcmplx)) {
		MPI_TYPE=MPI_DOUBLE_COMPLEX;
	} else {}
	//get the solution from all processor to root ***POTENTIAL BOTTLENECK***
	for(p=0; p<total_pop; ++p) {
		if(p%nb_proc!=0) { //if candidate is not in root, send the solution to root.
			if(my_rank==p%nb_proc) {
				this->pop[int(p/nb_proc)].read_best(input);
				MPI_Send(&input,this->num,MPI_TYPE,0,tag,MPI_COMM_WORLD);
			} else if(my_rank==0) {
				MPI_Recv(&input,this->num,MPI_TYPE,p%nb_proc,tag,MPI_COMM_WORLD,&status);
				for(i=0; i<this->num; ++i) {
					all_soln[p][i]=input[i];
				}
			}
		} else if(p%nb_proc==0) { //if candidate is in root, read the solution into the memory
			if(my_rank==0) {
				this->pop[int(p/nb_proc)].read_best(input);
				for(i=0; i<this->num; ++i) {
					all_soln[p][i]=input[i];
				}
			}
		}
	}// p loop

	MPI_Barrier(MPI_COMM_WORLD);

	//select family members for the candidate on the processor other that root
	for(p=0; p<total_pop; ++p) {
		if(my_rank==0) {
			//srand(time(NULL)+p);
			//temp=rand();//flush out the first sampling which is not random
			if(total_pop<=fam_size+1) {
				for(f=0; f<fam_size; ++f) {
					fam[f]=rand()%total_pop;
				}
			} else {
				do fam[0]=rand()%total_pop;
				while(fam[0]==p);
				do fam[1]=rand()%total_pop;
				while(fam[1]==p||fam[1]==fam[0]);
				do fam[2]=rand()%total_pop;
				while(fam[2]==p||fam[2]==fam[0]||fam[2]==fam[1]);
			}

			//cout<<"family="<<fam[0]<<","<<fam[1]<<","<<fam[2]<<endl;

			//create donor
			for(i=0; i<this->num; ++i) {
				//create donor
				input[i]=all_soln[fam[0]][i]+F*(all_soln[fam[1]][i]-all_soln[fam[2]][i]);
				//cross-over
				coin=double(rand())/RAND_MAX;
				if(coin<=Cr||i==rand()%this->num) {
					can[i]=input[i];
				}
				//if(coin<=Cr) {
				//    can[i]=input[i];
				//}
				else {
					can[i]=all_soln[p][i];
				}
			}
			//send the new candidate back
			MPI_Send(&can,this->num,MPI_TYPE,p%nb_proc,tag,MPI_COMM_WORLD);
		}//my_rank==0

		if(my_rank==p%nb_proc) { //receive and store the new candidates in contender
			MPI_Recv(&can,this->num,MPI_TYPE,0,tag,MPI_COMM_WORLD,&status);
			//this->prob->normalize(input);
			this->prob->modulo(can);
			this->pop[int(p/nb_proc)].update_cont(can);
		}

	}//p loop

	MPI_Barrier(MPI_COMM_WORLD);

	for(p=0; p<this->pop_size; ++p) {
		this->Cont_fitness(p);   //compute the fitness
	}

	MPI_Barrier(MPI_COMM_WORLD);
}

template <typename typeT>
void DE<typeT>::selection(int my_rank, int total_pop, int nb_proc)
{
	int p;
	for(p=0; p<this->pop_size; ++p) {
		if(this->pop[p].read_bestfit(0)<this->pop[p].read_contfit(0)) {
			this->pop[p].update_best();
			//cout<<"can select:"<<this->pop[p].read_bestfit()<<endl;
		} else {
			/*do nothing cout<<endl;*/
		}
	}
};

template <typename typeT>
void DE<typeT>::fit_to_global()
{
	int p;
	for(p=0; p<this->pop_size; ++p) {
		this->pop[p].put_to_global();
		//cout<<p<<":"<<this->pop[p].read_globalfit(0)<<","<<this->pop[p].read_globalfit(1)<<endl;
	}
}


#endif // DE_H
