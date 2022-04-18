#include "../message_list.h"

#include "../master6.h"

#include "parallel_cp_stuff.h"
#include "../input_list_creator.h"

#include <getopt.h>
#include <signal.h>




#define USING_IOPTIMIZE

#ifdef USING_IOPTIMIZE
	#include "start_point_generator6.h"
#endif

using namespace std;

int N_TASKS = 8;
int N_NODES = 1;
int N_PROCS = 20;
// int MAX_BUFFER_SIZE = 1024;
int N_SLAVES = 8;
int n_cores_per_slave = 1;
int parallellizzation_level = 8;
int change_to_random_obj = 0;

double IMP = 0.01;

bool USE_PARALLEL_CP =true;
// int basis_eval_function_id;
// int objective_function_id;
// int sample_function_id;
bool use_first_output = true;
double variable_perturbation = 1;

bool using_analytic_center = false;

int printlevel = 1;

bool bwait = false;
bool solve_original_centering_problem = false;

static constexpr unsigned default_seed = 5489u;

class Seeder2{
	unsigned seed;
	minstd_rand0 my_seed_seeder;
public:
	Seeder2(unsigned arg_seed = default_seed):seed(arg_seed){};

	unsigned nextSeed(){
		return ++seed;
	};
	unsigned nextSeed_lcm(){
		unsigned retval = seed;
		seed = my_seed_seeder();
		return retval;
	};

};

Seeder2 my_seed_generator(2);

Master<Pcp_input, Pcp_executable, Pcp_output> my_master(N_SLAVES);


bool master_already_spawn = false;

double starting_objective(bool set, double s){
	static double starting_objective_value = 1e30;

	if (set) {
		starting_objective_value = s;
	}

	return starting_objective_value;

}


void my_handler(int s){
	printf("Caught signal %d\n",s);


	if(s == 2){
		if (master_already_spawn){
			my_master.send_breaking_command_to_all_slaves(use_breaking_command_kill_all);
		}

		exit(1);
	}
	return;
}

input_list_creator list_creator;

#ifdef USING_IOPTIMIZE
int set_new_random_objective(Pcp_input &in, ioptimize_wrapper &iop, std::ostream &out =cout){
	vector<double> random_point;
// 	vector<double> random_point;
//  	iop.presetting_for_sampler(1,out);
	out << "center: ";
	for (int i = 0;i <  iop.get_current_center_point().size(); ++i){
		out << iop.get_current_center_point()[i] << " ";

	}
	out << endl;
	out << "randoms: ";
// 	iop.initialize

for(int i =0; i < 3; ++i ){
		cout << iop.get_next_random_point(random_point,out) << endl;;

	if (iop.get_current_center_point().size() != random_point.size()){
		out << "iop.get_current_center_point().size() != random_point.size() " << endl;
		out << iop.get_current_center_point().size() <<" != " << random_point.size() << endl;
		out <<"__LINE__: "<<__LINE__ << endl;
		return -1;
	}
	out << "random_point2: " ;
	for (int i = 0;i < random_point.size(); ++i){
// 		random_point[i] -= iop.get_current_center_point()[i];
		out << random_point[i] << " ";
	}
	out << endl;
}
// 	random_point[0] = 1;
	in.set_objective_function(random_point);

	return 0;

}
#endif

Pcp_output run(string filename, int failed_to_find_limit = 1, double pert_ = 0, double global_time_limit = -1, double iteration_time_limit = -1, std::ostream &out = std::cout){
	bool fcontinue = true;
	int failed_to_find_counter = 0;

	int use_input_type = use_input_type1;

	double gap;
	double relaxation_objective;

	int iteration_counter = 0;
	int counter = 0;

	time_t timestart1;
	time_t timeend1;
	time_t global_time_begin;
	time_t global_time_end;
	global_time_begin = time(NULL);

	vector<mySerializableRowCut> allrowcuts;
	vector<bool> cutapplied;
	vector<mySerializableRowCut> currentrowcuts;
	vector<mySerializableSolution> allsolutionlist;
	vector<mySerializableSolution> currentsolutionlist;

	mySerializableSolution best_solution;
	Pcp_output best_output(false);

	double best_integer_objective = starting_objective(false, 1);

	double best_noninteger_objective;
	bool best_incumbent_changed = false;

	bool original_objective_sense;

	double GAP;

	bool problem_infeasible = false;

	bool print = true;


	Cut_type t(2);

	#ifdef USING_IOPTIMIZE
	static int analytic_center_version = 0;
	ioptimize_wrapper iop_caller;
// 	vector<double> 	current_analytic_center;
	vector<double> current_lp_optimum;
	#endif

	#ifdef USING_IOPTIMIZE
	if(using_analytic_center){

		iop_caller.set_seed(1);
		iop_caller.setCenProbType(solve_original_centering_problem);
		cout << "solve_original_centering_problem "<< solve_original_centering_problem << endl;
		iop_caller.initialize(filename);


		bool hit_n_run= false, short_dikin = false, long_dikin =false;
		int walk_set = 0;
		int n_set = 0;

		if (list_creator.get_n_hit_and_run() > 0){
			hit_n_run = true;
			walk_set = 3;
			n_set ++;
		}
		if (list_creator.get_n_long_dikin() > 0){
			long_dikin = true;
			walk_set = 2;
			n_set ++;

		}
		if (list_creator.get_n_short_dikin() > 0){
			short_dikin = true;
			n_set ++;
			walk_set = 1;
		}
		iop_caller.set_walk_type(walk_set);

		if (n_set != 1 ){
			cout << "number of walk types needed is not 1 but " << n_set << endl;
			cout << "hit_n_run  : " << hit_n_run << endl;
			cout << "long_dikin : " << long_dikin<< endl;
			cout << "short_dikin: " << short_dikin<< endl;

		}

		iop_caller.update_sides(-1e20,1e20);
		out<< "calculating analytic center "<< endl ;
		cout<< "calculating analytic center "<< endl ;

		vector<double> sx;
		vector<double> sx2;

// 		iop_caller.setCenProbType(false);
// 		iop_caller.use_paramOpt();

 		iop_caller.calculate_analytic_center_self();
//  		out << __LINE__ << endl;
// 		out << " center point " << endl;
//
// 		vector_line_print<double>(iop_caller.get_current_center_point(),out);
// 		iop_caller.setCenProbType(true);
// 		iop_caller.use_paramOpt();
// 		iop_caller.calculate_analytic_center_self();
// 		out << __LINE__ << endl;
// 		out << "apprx center point " << endl;

// 		vector_print<double>(iop_caller.get_current_center_point(),out);

// 		iop_caller.setCenProbType();

		// 		iop_caller.init_sampler();
// 		out << __LINE__ << endl;
// 		iop_caller.get_next_random_point(sx);
// 		out << __LINE__ << endl;
//
// 		out << __LINE__ << endl;
		iop_caller.presetting_for_sampler();
// 		out << __LINE__ << endl;

// 		iop_caller.get_next_random_point(sx);
// 		out << "random point: " << endl;
// 		vector_line_print<double>(sx,out);
// 		iop_caller.get_next_random_point(sx);
// // 		out << "random point: " << endl;
// 		vector_line_print<double>(sx,out);
// 		out << __LINE__ << endl;

// 		out << "TRYING TO READ VIA ANOTHER IOP INSTANCE" << endl;
// 		out << "call calculate_analytic_center_to file " << endl;

// 		iop_caller.calculate_analytic_center_to_file("analytic_center_file.txt");

// 		out << __LINE__ << endl;
//
// 		ioptimize_wrapper reader;
// 		out << __LINE__ << endl;
//
// 		reader.set_seed(5489u);
// 		out << __LINE__ << endl;
// 		reader.setCenProbType(solve_original_centering_problem);
// 		out << __LINE__ << endl;
// 		reader.initialize(filename);
// 		out << __LINE__ << endl;
// // 		reader.use_paramOpt();
// 		out << __LINE__ << endl;
// 		reader.set_walk_type(walk_set);
// 		out << __LINE__ << endl;
// 		reader.update_sides(-1e20,1e20);
// 		out << __LINE__ << endl;
//
// 		reader.set_analytic_center_from_file("analytic_center_file.txt", out);
// 		out << __LINE__ << endl;
// 		out << "read center: " << endl;
// 		vector_line_print<double>(reader.get_current_center_point(),out);
//
//
//
// 		out << __LINE__ << endl;
// 		out << "reader.get_next_random_point(sx2,out): " << reader.get_next_random_point(sx2,out) << endl;;
// 		out << __LINE__ << endl;
// 		vector_line_print<double>(sx2,out);
//
// 		out << __LINE__ << endl;
// 		out << "reader.get_next_random_point(sx2,out): " << reader.get_next_random_point(sx2,out) << endl;;
// 		out << __LINE__ << endl;
// 		vector_line_print<double>(sx2,out);
//
// 		out << __LINE__ << endl;
//
// 		out << __LINE__ << endl;


	}

	#endif

// 	return best_output;


	//int counter = 0;
	while(fcontinue){
		if (global_time_limit > 0){
			global_time_end = time(NULL);
			if (global_time_limit < (global_time_end - global_time_begin)){
				cout << "time limit reached" << endl;
				break;
			}
		}
		if(USE_PARALLEL_CP){
			if(!master_already_spawn){

				N_SLAVES = list_creator.get_n_all();
				if (N_SLAVES == 0){
					list_creator.set_n_original_objective(8);
					N_SLAVES = 8;
				}
				parallellizzation_level = N_SLAVES;

				my_master = Master<Pcp_input, Pcp_executable, Pcp_output>(N_SLAVES);


				if(print) out << "(!master_already_spawn)" << endl;
				int n_spawns = my_master.default_spawn_via_pvm("parallel_cp_slave", NULL, 0, N_SLAVES);
				if(print) out << "(!master_already_spawn2) nspawns: " << n_spawns << endl;

				Options op(n_cores_per_slave,n_cores_per_slave,n_cores_per_slave,0,0,true,1,0,0,1);
				if(print) out << "(!master_already_spawn3)" << endl;

				my_master.send_single_option_to_slaves(op);
				if(print) out << "(!master_already_spawn4)" << endl;

 				vector<Pcp_input> inp_list = list_creator.create_input_list_CP<Pcp_input>();

// 				Pcp_input pcp_in(1);
// 				pcp_in.set_change_obj(change_to_random_obj);
// 				Cut_type t(2);
// 				pcp_in.set_cut_type(t);
// 				pcp_in.set_perturbation_level(pert_);
// 				pcp_in.set_serialization_size(0);

				if(print) out << "(!master_already_spawn5) parallellizzation_level: " << parallellizzation_level<< endl;

				for (int i = 0; i < parallellizzation_level; ++i){
					inp_list[i].OVPR_to_OBJECTIVE();
					Pcp_executable pcp_ex(filename,my_seed_generator.nextSeed_lcm());
					my_master.add_executable1(pcp_ex);

					inp_list[i].set_input_type(1);
					inp_list[i].set_cut_type(t);
					inp_list[i].set_serialization_size(0);
// 					inp_list[i].set_objective_function_size(0);
					inp_list[i].set_new_rhs_for_problem(best_integer_objective - IMP);

					#ifdef USING_IOPTIMIZE
					if(inp_list[i].get_obj_type() == OBJECTIVE_TYPE_RANDOM_FROM_HIT_RUN){
						iop_caller.set_walk_type(3);
						out << "inp_list["<<i<<"]: OBJECTIVE_TYPE_RANDOM_FROM_HIT_RUN "  << endl;
						set_new_random_objective(inp_list[i],iop_caller,out);

					}
					if(inp_list[i].get_obj_type() == OBJECTIVE_TYPE_RANDOM_FROM_LONG_DIKIN){
						iop_caller.set_walk_type(2);
						out << "inp_list["<<i<<"]: OBJECTIVE_TYPE_RANDOM_FROM_LONG_DIKIN "  << endl;

						set_new_random_objective(inp_list[i],iop_caller,out);

					}
					if(inp_list[i].get_obj_type() == OBJECTIVE_TYPE_RANDOM_FROM_SHORT_DIKIN){
						iop_caller.set_walk_type(1);
						out << "inp_list["<<i<<"]: OBJECTIVE_TYPE_RANDOM_FROM_SHORT_DIKIN"  << endl;
						set_new_random_objective(inp_list[i],iop_caller,out);

					}
					#endif

					my_master.add_input1(inp_list[i]);

				}

				if(print) out << "(!master_already_spawn6)" << endl;

				int code = my_master.distribute_input_and_executable_list(cout);
				if(print) out << "(!master_already_spawn7)" << endl;

				if(code>0){
					out << "problem in distributing distribute_input_and_executable_list code: " << code;
				}
				master_already_spawn = true;
				if(print) out << "(!master_already_spawn end)" << endl;

			}
			else{
				if(print) out << "ELSE (!master_already_spawn)" << endl;

				vector<Pcp_input> inp_list = list_creator.create_input_list_CP<Pcp_input>();


// 				Pcp_input new_input(1);
// 				new_input.set_change_obj(change_to_random_obj);
// 				new_input.set_cut_type(t);
// 				new_input.set_perturbation_level(pert_);
// 				new_input.add(currentrowcuts);
// 				new_input.set_input_type(1);
// 				new_input.set_apply_cuts(true);
				my_master.clear_input1_list();
				vector<mySerializableRowCut> cuts_to_be_added; // add all unused cuts

				if(print)  	out << "cutapplied.size() "<< cutapplied.size() << endl;
				int tttl = 0;
				for(int i =0; i < cutapplied.size();++i){
					if(!cutapplied[i]){
						tttl++;
						cuts_to_be_added.push_back(allrowcuts[i]);
						cutapplied[i] = true;
					}

				}
// 				if(print) out << "tttl "<< tttl << endl;
				if(print) out << "cuts_to_be_added.size() "<< cuts_to_be_added.size() << endl;



				for (int i = 0; i < parallellizzation_level; ++i){
					inp_list[i].OVPR_to_OBJECTIVE();
					inp_list[i].set_input_type(2);

					inp_list[i].set_cut_type(t);
					inp_list[i].set_serialization_size(0);
					inp_list[i].set_apply_cuts(true);

					if(best_incumbent_changed){
						inp_list[i].set_new_rhs_for_problem(best_integer_objective - IMP);
						//inp_list[i].set_input_type(2);
					}
					else {
						inp_list[i].set_new_rhs_for_problem(best_integer_objective + IMP);

					}

					inp_list[i].set_rowcutlist(cuts_to_be_added);
// 					if(print) out << "inp_list[" << i <<"].get_rowcutlist().size() " << inp_list[i].get_rowcutlist().size() << endl;
// 					if (inp_list[i].get_perturbation_level()>0) inp_list[i].set_perturbation_type(1);
					#ifdef USING_IOPTIMIZE
					if(inp_list[i].get_obj_type() == OBJECTIVE_TYPE_RANDOM_FROM_HIT_RUN){
						iop_caller.set_walk_type(3);
						cout << "inp_list["<<i<<"]: OBJECTIVE_TYPE_RANDOM_FROM_HIT_RUN "  << endl;
						set_new_random_objective(inp_list[i],iop_caller,out);

					}
					if(inp_list[i].get_obj_type() == OBJECTIVE_TYPE_RANDOM_FROM_LONG_DIKIN){
						iop_caller.set_walk_type(2);
						cout << "inp_list["<<i<<"]: OBJECTIVE_TYPE_RANDOM_FROM_LONG_DIKIN "  << endl;

						set_new_random_objective(inp_list[i],iop_caller,out);

					}
					if(inp_list[i].get_obj_type() == OBJECTIVE_TYPE_RANDOM_FROM_SHORT_DIKIN){
						iop_caller.set_walk_type(1);
						cout << "inp_list["<<i<<"]: OBJECTIVE_TYPE_RANDOM_FROM_SHORT_DIKIN"  << endl;
						set_new_random_objective(inp_list[i],iop_caller,out);

					}
					#endif
					my_master.add_input1(inp_list[i]);
				}

// 				if(best_incumbent_changed){
// 					new_input.set_new_rhs_for_problem(best_integer_objective);
// 					new_input.set_input_type(2);
// 				}

				my_master.send_breaking_command_to_all_slaves(use_continue_slave_clear_all_but_executable_and_command_list);
//
// 				for (int i = 0; i < my_master.get_input1_list().size(); ++i){
//                     cout << "printing inputlist i = " << i << endl;
//                     my_master.get_input1_list()[i].print();
//                     cout << "---xxx--" << endl;
// 				}
				if(print) out << "distributing inputlist" << endl;
				if(print) {
// 					out << "my_master.get_input1_list().size(): " << my_master.get_input1_list().size()<< endl;


				}
				if (my_master.distribute_input1_list() >0 ){
					out << "problem in distribute_input_list" << endl;
				}

				if(print) out << "restarting all slaves" << endl;
				my_master.restart_all_slaves();



			}

			if(print) out << "collecting all results " << endl;
			my_master.collect_results(0);
			if(print) out << "all results collected.. " << "sort_output1_list" << endl;

// 			my_master.sort_output1_list();

			vector<Pcp_output> outputlist;
			currentsolutionlist.clear();
			currentrowcuts.clear();
			failed_to_find_counter++;
			best_incumbent_changed = false;
			int size = my_master.get_output1_list().size();
			for (int i = 0; i < my_master.get_output1_list().size();++i){

// 				if((!my_master.get_output1_list()[i].get_is_incumbent()) && (!my_master.get_output1_list()[i].get_is_cut()) ){
// 					problem_infeasible = true;
// 					cout << "PROBLEM IS INFEASIBLE FROM OUTPUT " << i << " COUNTER " << counter << " ... breaking." << endl;
// 					break;
// 				}
				if((!my_master.get_output1_list()[i].get_is_feasible())){
					problem_infeasible = true;
					cout << "PROBLEM IS INFEASIBLE FROM OUTPUT " << i << " COUNTER " << counter << " ... breaking." << endl;
					break;
				}

// 			    cout << "printing outputlist for i = "<< i <<endl;
// 				my_master.get_output1_list()[i].print();
// 				cout << "--" <<endl;

				if (my_master.get_output1_list()[i].get_is_incumbent()){
					cout << "-- " << my_master.get_output1_list()[i].get_is_cut() << " " << my_master.get_output1_list()[i].get_is_incumbent() << " --" << endl;
					for (int j = 0; j < my_master.get_output1_list()[i].get_solutionlist().size();++j){
						mySerializableSolution as = my_master.get_output1_list()[i].get_solutionlist()[j];
						if (as.get_original_objective() < best_integer_objective){
							best_incumbent_changed = true;
							best_integer_objective = as.get_original_objective();
						}
					}

				}
				if (my_master.get_output1_list()[i].get_serialization_size_rowcutlist()>0){


					vector<mySerializableRowCut> asd = my_master.get_output1_list()[i].get_rowcutlist();
					currentrowcuts.insert(currentrowcuts.end(),asd.begin(),asd.end());


				failed_to_find_counter= 0;
				}


			}


			if (problem_infeasible){
				cout << "counter: " << counter << "\t Best solution so far: "  << best_integer_objective << endl;
				fcontinue = false;
				break;
			}
			vector<bool> temp(currentrowcuts.size(),false);
			allrowcuts.insert(allrowcuts.end(),currentrowcuts.begin(),currentrowcuts.end());
			cutapplied.insert(cutapplied.end(),temp.begin(),temp.end());


			my_master.clear_output1_list();

			if(bwait){
				cout << "asd currentrowcuts.size : " <<currentrowcuts.size()<< endl;
				sort(currentrowcuts.begin(),currentrowcuts.end());
				cout << "asd currentrowcuts.size : " <<currentrowcuts.size()<< endl;
				cout.flush();
				std::vector<mySerializableRowCut>::iterator it;
				cout << "asd currentrowcuts.size : " <<currentrowcuts.size()<< endl;
				cout.flush();

				it = unique(currentrowcuts.begin(),currentrowcuts.end());
				cout << "asd currentrowcuts.size : " <<currentrowcuts.size()<< endl;
				cout.flush();

// 				currentrowcuts.erase();
				currentrowcuts.resize( std::distance(currentrowcuts.begin(),it) ); // 10 20 30 20 10

				cout << "asd currentrowcuts.size : " <<currentrowcuts.size()<< endl;
				cout.flush();
				for (int i = 0; i < currentrowcuts.size(); ++i){
 					cout << i + 1 <<" : ";
 					currentrowcuts[i].print();
				}

				int tasd;
				std::cin >> tasd;

			}


		}
		else{
// 			Pcp_executable pcp_ex(filename,my_seed_generator.nextSeed_lcm());
// 			Pcp_input pcp_in;
// 			pcp_in.set_change_obj(change_to_random_obj);
// 			Cut_type t(1);
// 			pcp_in.set_cut_type(t);
// 			pcp_in.set_perturbation_level(pert_);
// 			pcp_in.set_serialization_size(0);


		}
// 		cout << "counter "<< counter << endl;
		counter++;

		cout << "counter: " << counter << "\t Best solution so far: "  << best_integer_objective << "\t nCuts: " << allrowcuts.size() << endl;

// 		if (counter > 5) break;

	}

	return best_output;

}


int main(int argc, char* argv[]){
	struct sigaction sigIntHandler;
	double iteration_time_limit = -1;
	double global_time_limit = -1;

	sigIntHandler.sa_handler = my_handler;
	sigemptyset(&sigIntHandler.sa_mask);
	sigIntHandler.sa_flags = 0;

	sigaction(SIGINT, &sigIntHandler, NULL);

	int opt;
	int failed_to_find_limit = 1;
	int perturbation_level = 0;
	string filename = "/home/utkukoc/mps_files/";
	string p_name;
	if (argc >1){

		filename += string(argv[1]);
		p_name = string(argv[1]);
	}
	double tmp = 0;
	while ((opt = getopt(argc, argv, "r:p:o:v:t:i:f:d:b:wS:1:2:cas:l:h:z:m:u:O:k")) != -1){
		switch (opt){
			case 'r':
				list_creator.set_n_random_objective(atoi(optarg));
				break;
			case 'p':
				list_creator.set_n_fixed_perturbation(atoi(optarg));
				break;
			case 'o':
				list_creator.set_n_original_objective(atoi(optarg));
				break;
			case 't':
				global_time_limit = atof(optarg);
				break;
			case 'i':
				iteration_time_limit = atof(optarg);
				break;
			case 'f':
				failed_to_find_limit = atoi(optarg);
				break;
			case 'd':
				list_creator.set_perturbation_level(atof(optarg));
				list_creator.set_perturbation_percentage(atof(optarg));
				break;
			case 'v':
				list_creator.set_n_variable_perturbation(atoi(optarg));
				break;
			case 'b':
				list_creator.set_perturbation_variation(atof(optarg));
				break;
			case 'w':
				bwait= true;
				break;
			case 'k':
				solve_original_centering_problem = true;
				break;
			case 'S':
				tmp = atof(optarg);
				tmp = starting_objective(true, tmp);
				break;
			case 'm':
				list_creator.set_n_random_objective_from_dikin(atoi(optarg));
				using_analytic_center = true;
				break;
			case 'u':
				list_creator.set_n_perturbed_objective_from_dikin(atoi(optarg));
				using_analytic_center = true;
				break;
			case 's':
				list_creator.set_n_short_dikin(atoi(optarg));
				using_analytic_center = true;
				// 				#define USING_IOPTIMIZE
				break;
			case 'l':
				list_creator.set_n_long_dikin(atoi(optarg));
				using_analytic_center = true;
				// 				#define USING_IOPTIMIZE
				break;
			case 'h':
				list_creator.set_n_hit_and_run(atoi(optarg));
				using_analytic_center = true;
				// 				#define USING_IOPTIMIZE
				break;
			case 'O':
// 				output_file_number = atoi(optarg);

				break;
			default:
				cerr << "Usage: " << argv[0] << " <test#>  -r: -p: -o: -v: -t: -i: -f: -d: -b: -w" << endl;
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

	run (filename);

	if (master_already_spawn){
		my_master.send_breaking_command_to_all_slaves(use_breaking_command_kill_all);
	}
	return 0;

}

// kate: indent-mode cstyle; indent-width 1; replace-tabs on; 
