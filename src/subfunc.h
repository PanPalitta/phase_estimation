/*supporting functions for main function*/
#include<fstream>

void output_header() {
    ofstream output_file;
    ofstream time_file;

    output_file.open("output.dat",ios::app);
    output_file<<"#N \t Sharpness \t Policy"<<endl;
    output_file.close();
    time_file.open("time.dat",ios::app);
    time_file<<"#N \t Time"<<endl;
    time_file.close();
}

template<typename typeT>
void output_result(int num, double final_fit, typeT *solution, time_t start_time) {
    int i;
    ofstream output_file;
    ofstream time_file;

    output_file.open("output.dat",ios::app);
    output_file<<num<<"\t"<<final_fit<<"\t";
    for(i=0; i<num; ++i) {
        output_file<<solution[i]<<"\t";
    }
    output_file<<endl;
    output_file.close();

    time_file.open("time.dat",ios::app);
    time_file<<num<<"\t"<<difftime(time(NULL),start_time)<<endl;
    time_file.close();

}

