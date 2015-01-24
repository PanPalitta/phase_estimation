#ifndef CANDIDATE_H
#define CANDIDATE_H

#include <cstring>

typedef complex<double> dcmplx;

template<typename typeT>
class Candidate
{
public:
    Candidate() {};
    ~Candidate();

    void init_can(int numvar);
    void init_velocity();

    void update_cont(typeT *input);
    void update_vel(typeT *input);
    void update_best();
    void update_global(typeT *input);
    void put_to_global();

    void read_cont(typeT *output);
    void read_vel(typeT *output);
    void read_best(typeT *output);
    void read_global(typeT *output);

    double read_contfit() {
        return cont_fit;
    }
    double read_bestfit() {
        return best_fit;
    }
    double read_globalfit() {
        return global_fit;
    }
    int read_bestt() {
        return best_times;
    }

    void write_contfit(double fit,int tt) {
        cont_fit=fit/tt;
        times=tt;
    }
    void write_bestfit(double fit) {
        best_fit=best_fit*best_times+fit;
        best_times+=1;
        best_fit=best_fit/best_times;
    }
    void write_globalfit(double fit) {
        global_fit=fit;
    }

private:
    int num;
    typeT *can_best;
    typeT *contender;
    typeT *velocity;
    typeT *global_best;
    double best_fit,cont_fit,global_fit; //should this be double or another type?
    int times, best_times, global_times; //number of samples used to calculate average best_fit

};

template<typename typeT>
Candidate<typeT>::~Candidate() {
    delete[] can_best;
    delete[] contender;
    delete[] velocity;
    delete[] global_best;
}

template<typename typeT>
void Candidate<typeT>::init_can(int numvar) {
    num=numvar;
    can_best=new typeT[num];
    contender=new typeT[num];
    global_best=new typeT[num];
}

template<typename typeT>
void Candidate<typeT>::init_velocity() {
    velocity=new typeT[num];
}


template<typename typeT>
void Candidate<typeT>::update_cont(typeT *input) {
    memcpy(contender, input, num*sizeof(typeT));
}

template<typename typeT>
void Candidate<typeT>::update_vel(typeT *input) {
    int i;
    for (i=0; i<num; i++) {
        velocity[i]=input[i];
    }
}

template<typename typeT>
void Candidate<typeT>::update_best() {
    int i;
    for (i=0; i<num; i++) {
        can_best[i]=contender[i];
    }
    best_fit=cont_fit;
    best_times=times;
}

template<typename typeT>
void Candidate<typeT>::update_global(typeT *input) {
    int i;
    for (i=0; i<num; i++) {
        global_best[i]=input[i];
    }
}

template<typename typeT>
void Candidate<typeT>::put_to_global() {
    int i;
    for (i=0; i<num; i++) {
        global_best[i]=can_best[i];
    }
    global_fit=best_fit;
    global_times=best_times;
}

template<typename typeT>
void Candidate<typeT>::read_cont(typeT *output) {
    int i;
    for (i=0; i<num; i++) {
        output[i]=contender[i];
    }
}

template<typename typeT>
void Candidate<typeT>::read_vel(typeT *output) {
    int i;
    for (i=0; i<num; i++) {
        output[i]=velocity[i];
    }
}

template<typename typeT>
void Candidate<typeT>::read_best(typeT *output) {
    int i;
    for (i=0; i<num; i++) {
        output[i]=can_best[i];
    }
}

template<typename typeT>
void Candidate<typeT>::read_global(typeT *output) {
    int i;
    for (i=0; i<num; i++) {
        output[i]=global_best[i];
    }
}


#endif // CANDIDATE_H
