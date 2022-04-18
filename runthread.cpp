#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <pthread.h>
#include <signal.h>
#include "md5.h"

using namespace std;

int number_of_active_threads = 0;
pthread_mutex_t var_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t io_mutex = PTHREAD_MUTEX_INITIALIZER;
bool listing_ = false;

struct run_parameters
{
    string fname;
    string command;
};

bool file_exists(string &file_name)
{
    struct stat file_info;

    if(stat(file_name.c_str(),&file_info) == 0)
        return true;
    else
        return false;
}



void* run(void *arg)
{
    run_parameters *p = (run_parameters*) arg;
    string ip = p->fname + md5(p->command) + ".ip";
    string fn = p->fname + md5(p->command) + ".fn";

    pthread_mutex_lock(&io_mutex);
    ofstream file1(ip.c_str());
    file1 << 0;
    file1.flush();
    file1.close();
    pthread_mutex_unlock(&io_mutex);

    cout << p->command.c_str() << endl;
	if(!listing_){
		if (system(p->command.c_str()) == 0)
		{
			pthread_mutex_lock(&io_mutex);
			ofstream file2(fn.c_str());
			file2 << 0;
			file2.close();
			unlink(ip.c_str());
			pthread_mutex_unlock(&io_mutex);
		}
	}
	else{
		pthread_mutex_lock(&io_mutex);
		ofstream file2(fn.c_str());
		file2 << 0;
		file2.close();
		unlink(ip.c_str());
		pthread_mutex_unlock(&io_mutex);
	}
    pthread_mutex_lock(&var_mutex);
    number_of_active_threads--;
    pthread_mutex_unlock(&var_mutex);
    pthread_exit(NULL);
    return NULL;
}

void my_handler(int s){
	printf("MULTI_THREAD_RUNNER Caught signal %d\n",s);
	int cinin = 0;
	cin >> cinin;

	switch (cinin){
		case 1:
			return;
			break;
		case 2:
			exit(1);
			break;
		case 3:
			abort();
			break;
		default :
			return;
	}

	{


		exit(1);
	}
	return;
}
int main(int argc, char *argv[])
{


bool wait_ = false;


struct sigaction sigIntHandler;


double iteration_time_limit = -1;
double global_time_limit = -1;

sigIntHandler.sa_handler = my_handler;
sigemptyset(&sigIntHandler.sa_mask);
sigIntHandler.sa_flags = 0;

sigaction(SIGINT, &sigIntHandler, NULL);

    ifstream command_list_file(argv[1]);

	int continue_from = 0;
    vector<string> commands;
    string s;
    int n_threads = 1;

    int opt;
    bool fForce = false, fClean = false, fRepeat = false, run_backwards = false;

    int sleeptime = 3000;

    while ((opt = getopt(argc, argv, "rcfm:s:d:wl")) != -1)
    {
        switch (opt)
        {
		case 'r':
			 fRepeat = true;
			 break;
        case 'c':
            fClean = true;
            break;
        case 'f':
            fForce = true;
            break;
        case 'm':
            n_threads = atoi(optarg);
            break;
        case 's':
            sleeptime = atoi(optarg);
            break;
		case 'd':
			continue_from = atoi(optarg);
			break;
		case 'w':
			wait_ = true;
			break;        
		case 'l':
			listing_ = true;
			break;
// 		case 'b':
// 			run_backwards = true;
// 			break;

        default:
            cerr << "Usage: " << argv[0] << "  <file_name> [-c -f -m <# of threads> -s <sleep time> -d <continue from>]" << endl;
            cout << "Continue anyway [Y/N] ? ";
            cout.flush();
            char c;
            cin >> c;
            if (c != 'Y' || c != 'y')
            {
                return 1;
            }
        }
    }

    while (getline(command_list_file, s))
        commands.push_back(s);

	if (continue_from <0) continue_from = 0;
	if(continue_from>=commands.size()) continue_from = commands.size() -1;



    string folder_name(".runnertemp/");
    mkdir(folder_name.c_str(), S_IWOTH || S_IROTH);

    if (fClean)
    {
        for (unsigned i = continue_from; i < commands.size(); ++i)
        {
            string ip = folder_name + md5(commands[i]) + ".ip";
            string fn = folder_name + md5(commands[i]) + ".fn";
            if (file_exists(fn)) unlink(fn.c_str());
            if (file_exists(ip)) unlink(ip.c_str());
        }
        return 0;
    }

    int temp;
	unsigned command_id = continue_from;
// 	if (run_backwards) command_id = commands.size() -1;
    int dummy;
	int n_continue = 1;
	int n_pass = 0;
    while (true)
    {
        pthread_mutex_lock(&var_mutex);
        temp = number_of_active_threads;
        pthread_mutex_unlock(&var_mutex);
        if (temp < n_threads)
        {
            if (command_id < commands.size())
            {
                string ip = folder_name + md5(commands[command_id]) + ".ip";
                string fn = folder_name + md5(commands[command_id]) + ".fn";

                pthread_mutex_lock(&io_mutex);
                bool flag = file_exists(fn);
                pthread_mutex_unlock(&io_mutex);
                if (flag && !fRepeat)
                {
                    command_id++;
                    continue;
                }
                pthread_mutex_lock(&io_mutex);
                flag = file_exists(ip);
                pthread_mutex_unlock(&io_mutex);
                if (flag && !fForce)
                {
                    command_id++;
                    continue;
                }
                run_parameters arg;
                arg.fname = folder_name;
                arg.command = commands[command_id];
                pthread_t th;
                if ((wait_)&&(command_id>continue_from)&&(n_continue<=0)){
                    cout << " Enter number of instances to continue" << endl;
					cin >> n_continue;
					if (n_continue == 0){
						cout << " Enter number of instances to pass" << endl;
						
						cin >> n_pass ; 
						command_id+= n_pass;
					}
//                     cin >> dummy;
                }
                n_continue--;
				cout << "________________________________________________________" << endl;
				cout << "MULTI_THREAD_RUNNER: RUNNING COMMAND " << command_id +1  << " OF " << commands.size() << " OF FILE "<<  argv[1] << endl;
				cout << "________________________________________________________" << endl;
				pthread_create(&th, NULL, run, &arg);
                pthread_mutex_lock(&var_mutex);
                number_of_active_threads++;
                pthread_mutex_unlock(&var_mutex);
                command_id++;
            }
        }
        pthread_mutex_lock(&var_mutex);
        temp = number_of_active_threads;
        pthread_mutex_unlock(&var_mutex);
        if ( (temp == 0) && (command_id >= commands.size()) )
        {
			if (fRepeat) {
			  command_id = 0;
			  continue;
			}
            else break;
        }
        usleep(sleeptime);
    }
    pthread_exit(NULL);
    return 0;
}
