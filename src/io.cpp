#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <map>

using namespace std;

void output_header( char const *output_filename, /*!<Pointer to output file where the sharpness will be printed.*/
                    char const *time_filename /*!<Pointer to output file where the clock time for optimizing a policy.*/
                  ) {
    ofstream output_file;
    ofstream time_file;

    output_file.open(output_filename, ios::app);
    output_file << "#N \t Sharpness \t Mean \t Policy" << endl;
    output_file.close();
    time_file.open(time_filename, ios::app);
    time_file << "#N \t Time" << endl;
    time_file.close();
    }

void output_result( int num, /*!< Number of variables in a policy.*/
		    int num_fit, /*!<number of fitness values*/
                    double *final_fit, /*!<fitness values.*/
                    double *solution, /*!< Policy*/
                    time_t start_time, /*!< Clock time for finding the policy.*/
                    char const *output_filename, /*!< Pointer to file that store the fitness value and the policy.*/
                    char const *time_filename /*!< Pointer to file that store time.s*/
                  ) {
    int i;
    ofstream output_file;
    ofstream time_file;

    output_file.open(output_filename, ios::app);
    output_file << num << "\t"; 
    for(i = 0; i< num_fit; i++){
	output_file << final_fit[i] << "\t";
    }
    for(i = 0; i < num; ++i) {
        output_file << solution[i] << "\t";
        }
    output_file << endl;
    output_file.close();

    time_file.open(time_filename, ios::app);
    time_file << num << "\t" << difftime(time(NULL), start_time) << endl;
    time_file.close();

    }

void output_result(int num, double final_fit, double *solution,
                   time_t start_time,
                   char const *output_filename,
                   char const *time_filename);

map<string, string> parse(ifstream & cfgfile) {
    map<string, string> options;
    string id, eq, val;

    while(cfgfile >> id >> eq >> val) {
        if (id[0] == '#') continue;  // skip comments
        //cout << id << " " << eq << " " << val << endl;
        if (eq != "=") throw runtime_error("Parse error");

        options[id] = val;
        }
    return options;
    }

void read_config_file(  char const *filename, /*!< File that contains program parameters.*/
                        int *pop_size, int *N_begin,
                        int *N_cut, int *N_end, int *iter, int *iter_begin,
                        int *repeat, int *seed, string *output_filename,
                        string *time_filename, string *optimization,
                        double *prev_dev, double *new_dev, double *t_goal, int *data_end) {
    /*! This fucntion reads the user specified configuration file or set the parameters to defaults.
    Parameters that can be set by users are
    population size: *pop_size
    smallest number of variables to be searched: *N_begin
    number of variables to switch from uniform initialization: *N_cut
    largest number of variables: *N_end
    number of iteration when initialization is uniformed: *iter_begin
    number of iteration before accepting/checking a policy: *iter
    number of times the fitness value is compute in the final selection process: *repeat
    seed for rng: *seed
    the name of the file that the fitness value and policy will be stored: *output_filename
    the name of the file that the clock time will be stored: *time_filename = "time.dat"
    name of the optimization algorithm: *optimization = "de"
    parameter for initialization using previous result -- for variables with previous data: *prev_dev = 0.01;
    parameter for initialization using previous result -- for new variables: *new_dev = 0.25;
    goal to be used in accept/reject criterion *t_goal
    the number of variable to begin using accept/reject criterion: *data_end
    */
    // Setting defaults
    *pop_size = 12;
    *N_begin = 4;
    *N_cut = 5;
    *N_end = 10;
    *iter = 100;
    *iter_begin = 300;
    *repeat = 10;
    *seed = time(NULL);
    *output_filename = "output.dat";
    *time_filename = "time.dat";
    *optimization = "de";

    *data_end = 8;
    *prev_dev = 0.01;
    *new_dev = 0.25;
    *t_goal = 0.98; //probability for calculating quantile

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
            }
        else if (it->first == "N_begin") {
            istringstream (it->second) >> *N_begin;
            }
        else if (it->first == "N_cut") {
            istringstream (it->second) >> *N_cut;
            }
        else if (it->first == "N_end") {
            istringstream (it->second) >> *N_end;
            }
        else if (it->first == "iter") {
            istringstream (it->second) >> *iter;
            }
        else if (it->first == "iter_begin") {
            istringstream (it->second) >> *iter_begin;
            }
        else if (it->first == "repeat") {
            istringstream (it->second) >> *repeat;
            }
        else if (it->first == "output_filename") {
            *output_filename = it->second;
            }
        else if (it->first == "time_filename") {
            *time_filename = it->second;
            }
        else if (it->first == "random_seed") {
            istringstream (it->second) >> *seed;
            }
        else if (it->first == "optimization") {
            *optimization = it->second;
            if (*optimization != "de" && *optimization != "pso") {
                throw runtime_error("Unknown optimization algorithm");
                }
            }
        else if (it->first == "data_end") {
            istringstream (it->second) >> *data_end;
            }
        else if (it->first == "prev_dev") {
            istringstream (it->second) >> *prev_dev;
            }
        else if (it->first == "new_dev") {
            istringstream (it->second) >> *new_dev;
            }
        else if (it->first == "t_goal") {
            istringstream (it->second) >> *t_goal;
            }
        else {
            throw runtime_error("Unknown configuration option");
            }
        }
    if (*N_cut < *N_begin) {
        throw runtime_error("Select new N_cut");
        }
    }
