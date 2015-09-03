#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <map>

using namespace std;

void output_header(char const *output_filename, char const *time_filename) {
    ofstream output_file;
    ofstream time_file;

    output_file.open(output_filename, ios::app);
    output_file << "#N \t Sharpness \t Policy" << endl;
    output_file.close();
    time_file.open(time_filename, ios::app);
    time_file << "#N \t Time" << endl;
    time_file.close();
    }

template<typename typeT>
void output_result(int num, double final_fit, typeT *solution, 
                   time_t start_time, char const *output_filename,
                   char const *time_filename) {
    int i;
    ofstream output_file;
    ofstream time_file;

    output_file.open(output_filename, ios::app);
    output_file << num << "\t" << final_fit << "\t";
    for(i = 0; i < num; ++i) {
        output_file << solution[i] << "\t";
        }
    output_file << endl;
    output_file.close();

    time_file.open(time_filename, ios::app);
    time_file << num << "\t" << difftime(time(NULL), start_time) << endl;
    time_file.close();

    }

template void output_result<double>(int num, double final_fit, double *solution, 
                                    time_t start_time, 
                                    char const *output_filename,
                                    char const *time_filename);

map<string, string> parse(ifstream & cfgfile)
{
    map<string, string> options;
    string id, eq, val;

    while(cfgfile >> id >> eq >> val)
    {
      if (id[0] == '#') continue;  // skip comments
      //cout << id << " " << eq << " " << val << endl;
      if (eq != "=") throw runtime_error("Parse error");

      options[id] = val;
    }
    return options;
}

void read_config_file(char const *filename, int *pop_size, int *N_begin, 
                      int *N_cut, int *N_end, int *iter, int *iter_begin, 
                      int *repeat, int *seed, string *output_filename, 
                      string *time_filename) {
    // Setting defaults
    *pop_size = 20;
    *N_begin = 4;
    *N_cut = 5;
    *N_end = 10;
    *iter = 100;
    *iter_begin = 300;
    *repeat = 10;
    *seed = time(NULL);
    *output_filename = "output.dat";
    *time_filename = "time.dat";
    // Parsing config file if there is one
    if (filename == NULL) return;
    ifstream config_file(filename);
    if(!config_file.is_open()) throw runtime_error("Config file cannot be opened!");
    map<string, string> options = parse(config_file);
    config_file.close();
    map<string, string>::iterator it;
    for (it = options.begin(); it != options.end(); ++it) {
       if (it->first == "pop_size") {
           istringstream (it->second) >> *pop_size;
       } else if (it->first == "N_begin") {
           istringstream (it->second) >> *N_begin;
       } else if (it->first == "N_cut") {
           istringstream (it->second) >> *N_cut;
       } else if (it->first == "N_end") {
           istringstream (it->second) >> *N_end;
       } else if (it->first == "iter") {
           istringstream (it->second) >> *iter;
       } else if (it->first == "iter_begin") {
           istringstream (it->second) >> *iter_begin;
       } else if (it->first == "repeat") {
           istringstream (it->second) >> *repeat;
       } else if (it->first == "output_filename") {
           *output_filename = it->second;
       } else if (it->first == "time_filename") {
           *time_filename = it->second;
       } else if (it->first == "random_seed") {
           istringstream (it->second) >> *seed;
       } else {
           throw runtime_error("Unknown configuration option");
       }
    }
    if (*N_cut < *N_begin) {
        throw runtime_error("Select new N_cut");
    }
}
