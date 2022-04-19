
//#define USING_MPI

#include "p2_definitions.h"
#include <sys/stat.h>
#include <IpTNLP.hpp>
#include <sys/resource.h>
#include <OsiConicSolverInterface.hpp>
#include <OsiIpoptSolverInterface.hpp>
//#include "start_point_generator6.h"

#define  history_size_limit_l 2000
#define  history_size_limit_u 2001

// 

/** MACRO NEEDED FOR calculate AC, create-distribute objective function or calculate and distribute AC */


// #define DEBUG_ON_THE_FLY

// // #define REDUCE_ALPHA alpha = floor(alpha*1e6*alpha_reduction)*1e-6;
#define REDUCE_ALPHA alpha = alpha*alpha_reduction;

// #define slave_optimality (slave_optimality_learned || slave_optimality_proven)


using namespace std;

// #define PPOF_ne

#define PRINT_LB_UB(runner)\
cout <<"line " << __LINE__ << " lb:" << runner.get_objective_cutoff_constraint_lb() << " ub: " << runner.get_objective_cutoff_constraint_ub() <<endl;


#ifdef PPOF_ne

#define PPOF(out, line, iteration) \
	out << "iteration: "<< iteration << " line: "<< line <<endl;
#else
	#define PPOF(out, line, iteration)
#endif


#define DEBUG_VIA_PRINTING_IN_EFFECT

#ifdef DEBUG_VIA_PRINTING_IN_EFFECT

	#define LP_WRITE(runner, stage, iter, line)\
		stringstream z;\
		z<< "stage_"<< stage << "_iter_" << iter << "_line_" << line ;\
		cout << "FILE_WRITE exporting " << z.str() <<endl;\
		runner.puke_lp_file(z.str().c_str());
#else

	#define LP_WRITE(runner, stage, iter, line)

#endif


#define USE_COMM_V2 \
if ((FP_total_iteration_counter - latest_communication_iteration) > communicate_every_x_iterations){ \
	runner_updated_from_other_slaves = slave_FP_communicate_with_master(runner, will_send, best_incumbent_solution_at_slave_FP, best_incumbent_objective_at_slave_FP);\
	latest_communication_iteration = FP_total_iteration_counter;\
}\
if(runner_updated_from_other_slaves == 2){\
	/** need to generate a new starting solution if not original objective*/\
	external_iteration_counter = 0;	\
	/**previous_rounded_solution= rounded_solution;*/\
	runner_updated_from_other_slaves = 0; \
}

#define USE_COMM_V3 \
if (slave_FP_communicate_with_master_v3(runner, FP_total_iteration_counter, best_incumbent_solution_at_slave_FP, best_incumbent_objective_at_slave_FP,will_send) ==0) external_iteration_counter =0;


#define CREATE_OBJ_AND_OPTIMIZE(runner) \
	l1_norm_timer.restart();\
	My_objective_function ob = runner.f_type_stage_1_or_2_create_obj_for_FP(rounded_solution,alpha,current_stage);\
	runner.set_FP_objective_function(ob);\
	runner.unset_bounds_FP_stage_2();\
	if (current_stage == 2) runner.f_type_stage_2_set_bounds(rounded_solution);\
	int res2 = runner.optimize();\
	l1_norm_timer.pause();\
	if (res2<=0){\
		/** to be done*/\
		/** LP not optimal */\
		cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << " LP not optimal res: "<< res2 << endl;\
		break;\
	}

#define UPDATE_FP_STATUS \
if(curFrac_d < minFracDbl){\
	if (curFrac_d/minFracDbl < 0.9)\
		missedDecr = 0;\
	minFracDbl = curFrac_d;\
	minFracIt = FP_total_iteration_counter;\
	best_starting_solution_so_far = starting_solution;\
}\
else\
	missedDecr++;

#define UPDATE_FP_STATUS_r \
if(curFrac_d < minFracDbl){\
	if (curFrac_d/minFracDbl < 0.9)\
		missedDecr = 0;\
		minFracDbl = curFrac_d;\
		minFracIt = FP_total_iteration_counter;\
		best_starting_solution_so_far = starting_solution;\
		best_rounded_solution_so_far = rounded_solution;\
		}\
	else\
		missedDecr++;
		
#define BREAK_STAGE_ONE\
	if(curFrac_i == 0) {\
		if(printing_breaks) cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << " binary variables are feasible at iteration: " << FP_total_iteration_counter << endl;\
		break; /** bnary variables are feasible*/\
	}\
	if (missedDecr > maxMissedDecr && stage2_iterlim>0) {\
		if(printing_breaks) cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << " too many itarations without 10% fractionality improvement stage1 " <<  endl;\
		break; /** too many itarations withour 10% fractionality improvement */\
	}\
	if (stage_1_restarts > 100) {\
		if(printing_breaks) cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << " too many stage_1_restarts" <<endl;\
		break;\
	}

#define BREAK_STAGE_TWO \
	if (missedDecr > maxMissedDecr){ \
		if(printing_breaks) cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << " too many itarations withour 10% fractionality improvement stage 2 " <<endl;\
		break; /** too many itarations withour 10% fractionality improvement */ \
	}\
	if (stage_2_restarts>100) {\
		if(printing_breaks) cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ <<  " too many stage_2_restarts " <<endl;\
		break; /** too many restarts     */\
	}
				
				
#define HISTORY_CHECK_STAGE(stage_i_restarts)\
	My_solution pre;\
	bool insert = history_insert_version1(rounded_solution,alpha,pre) ;\
	while(!insert && stage_i_restarts < 10*local_iterlim){\
		stage_i_restarts++;\
		if (runner.f_type_restart_1_or_2(starting_solution, rounded_solution, previous_rounded_solution, current_stage, FP_total_iteration_counter) ==0){\
			continue;\
		}\
		insert =   history_insert_version1(rounded_solution,alpha,pre) ;\
	}				
				


#define GLOBAL_TIME_CHECKER \
	if (global_break_timer.stop() > max_time_limit){\
		FINALIZE \
		cout << "slave " << SLAVEID <<  " aborting " <<endl;\
		abort();\
	}

#define OPTIMALITY_CHECK(runner)\
	if(runner.get_objective_cutoff_constraint_lb() > runner.get_objective_cutoff_constraint_ub()){\
		slave_optimality = true;\
	}
	
	
#define UPDATE_RUNNER_FP(best_inc, runner)\
	if (best_inc < latest_updated_fp_objective_in_slave){\
		latest_updated_fp_objective_in_slave = best_inc;\
		double objective_imp = ( best_inc  - runner.get_objective_cutoff_constraint_lb() ) * objective_improvement_percentage_FP;\
		if (objective_imp < objective_improvement_coefficient_FP || (!USE_RELATIVE_IMPROVEMENT_IN_FP) ){\
			objective_imp = objective_improvement_coefficient_FP ;\
		}\
		runner.update_objective_cutoff_constraint_ub ( best_inc -objective_imp );\
	}
	
	//cout << "best_inc :" <<best_inc << " obj imp: " << objective_imp << " USE_RELATIVE_IMPROVEMENT_IN_FP " << USE_RELATIVE_IMPROVEMENT_IN_FP << endl;\
	//cout << "runner.get_objective_cutoff_constraint_lb()" << runner.get_objective_cutoff_constraint_lb() << endl;\
	//cout << "objective_improvement_percentage_FP" << objective_improvement_percentage_FP << endl;\
	
	
#define UPDATE_RUNNER_CP(best_inc, LLB, runner)\
	double objective_imp = ( best_inc  - LLB ) * objective_improvement_percentage_CP;\
	if (objective_imp < objective_improvement_coefficient_CP || (!USE_RELATIVE_IMPROVEMENT_IN_CP) ){\
		objective_imp = objective_improvement_coefficient_CP ;\
		}\
	runner.update_objective_cutoff_constraint_ub ( best_inc -objective_imp );
		
	
 	//#define __WAIT
	#ifdef __WAIT

	int tempss = 1;
	#define DEB_WAIT(_line, _call,_iter) \
	cout << "iter " << _iter <<" call " << _call << " line: "<< _line << endl; \
	if (tempss <= 3) cin >> tempss;\
	if (tempss > 3) tempss--;\
	if(tempss ==1)	{starting_solution.shortlineprint_selected_indices(integer_columns_in_slave,false,cout);; cout << "---------" << endl;}\
		if(tempss ==2)	{rounded_solution.shortlineprint_selected_indices(integer_columns_in_slave,false,cout); cout << "---------" << endl;}\
			if(tempss ==3)	{\
				starting_solution.shortlineprint_selected_indices(integer_columns_in_slave,false,cout);\
				cout << "--" <<endl;\
				rounded_solution.shortlineprint_selected_indices(integer_columns_in_slave,false,cout);\
				cout << "---------" << endl;\
				}

	#else
		#define DEB_WAIT(_line, _call,_iter)
	#endif


	
/** global variables */
Timer global_time;
int xxx =0;
int _ux = 0;
int master_iop_n_threads =0;
My_solution best_incumbent_solution;

//OsiTMINLPInterface * nlp = NULL;

// bool temp_opt_tester = false;

Timer global_break_timer;

Timer time_to_optimize_in_master;
Timer time_to_optimize_in_slave;
Timer time_to_generate_cuts;
Timer time_to_apply_cuts_master;
Timer time_to_apply_cuts_slave;
Timer time_to_calculate_AC;
Timer time_to_create_RW_points_master;
Timer time_to_create_RW_points_slave;
Timer time_to_optimize_in_serial;
Timer time_to_apply_cuts_serial;
Timer time_to_create_starting_points_for_FP_slave;
Timer master_level_FP_timer;
Timer l1_norm_timer;
Timer time_for_CP;
Timer time_for_FP;

Timer time_temp;

Timer time_checker_for_FP_temp;

int communicate_every_x_iterations = -1;

string p_name, global_filename;

bool solve_original_centering_problem = false;
int type = 2; /** THIS REPRESENTS THE RANDOM WALK TYPE WE WILL USE*/

int number_of_cuts_per_nonoriginal_slave = -5;
int number_of_iterations = -10;

/** 1: FP 2: CP*/
vector<int> CP_FP_vector; 
int FP_iterations = -1;
int CP_iterations = -1;

int efficiency_threshold = -3;
double timelimit = -60;
double FP_time_limit = -60;
double FP_max_time_limit = 3600*24*360;
double CP_time_limit = -10;
// double SWAP_TIME = 10;

double aggressive_flip_probability = -1;

int experiment_number = -1;
int setting_number = -1;

int FP_min_iter_lim = -500;
int FP_max_iter_lim = 2043514880;
int CP_iter_lim = -1;


bool slave_optimality = false;
// bool slave_optimality_proven = false;
// bool slave_optimality_learned = false;

vector<int> external_column_types;

bool USE_RELATIVE_IMPROVEMENT_IN_CP = false;
bool USE_RELATIVE_IMPROVEMENT_IN_FP = true;

double objective_improvement_coefficient_FP = -0.1;
double objective_improvement_percentage_FP = -10.00;
double objective_improvement_coefficient_CP = -0.1;
double objective_improvement_percentage_CP = -10.00;

bool stop_fp_when_first_solution_is_found = false;

#ifdef DEBUGGING
static int _x_ =0;
#endif
bool blockingWait = false;
bool big_iteration_change = false;

unsigned slave_seed;
double starting_objective = 1e30;

double alpha_reduction = -0.9;
double initial_alpha = -1;
double alpha =-1;
double alpha_dist =-0.005;
int slave_FP_type = 3;
int initial_rounding_type = -8;

int cout_vector_time = 1;

static constexpr unsigned default_seed = 5489u;
unsigned initial_seed = default_seed;

Seeder my_seed_generator ( initial_seed );

Objective_type_distributor type_dist_for_CP;
Objective_type_distributor type_dist_for_FP;

stringstream summarystring;
stringstream summarystring2;


/** 1, use minimum iteration limit
 *  2, use time limit
 */
int use_what_limit = 1;


int solver_type = 1;



bool temporarily_using_original_objective = false;
bool slave_generates_own_functions_for_FP = false;
bool slave_generates_own_functions_for_CP = false;


vector<bool> miniterlim_checker;
vector<int> slave_iteration_numbers;

vector<int> latest_slave_commands;

vector<int> FP_CP_past;

#ifdef USING_IOPTIMIZE
ioptimize_wrapper master_iop;
ioptimize_wrapper slave_iop;


bool using_AC = false;
double master_ub_for_iop= 1e30;
double master_lb_for_iop= -1e30;
double slave_ub_for_iop= 1e30;
double slave_lb_for_iop= -1e30;


shadow_copy_ioptimize_wrapper latest_AC_copy_by_master;
shadow_copy_ioptimize_wrapper latest_AC_copy_recv_by_slave;

vector<shadow_copy_ioptimize_wrapper> shadow_copy_of_all_calculated_ACs;
vector<shadow_copy_ioptimize_wrapper> shadow_copy_of_all_received_ACs;



#endif
double max_time_limit = 0;

bool running_serial = false;

bool AC_v2 = true;

int master_iop_AC_needs_to_calculated = 0;
                                            /** ONLY TO BE UPDATED UNDER LOCK
                                             *  ** -99, not even initialized
                                             *   -1, there is no need to calculate
                                             *    0, being calculated right now
                                             *    1, needs to be calculated
                                             *    2, is being distributed
											 **/
								
//OsiIpoptSolverInterface *ipopt_solver;
OsiConicSolverInterface * ipopt_solver;

My_solution latest_lp_opt_sent_to_slaves;
double latest_updated_fp_objective_in_slave = 1e31;
											
int total_number_of_ACs_computed_in_master = 0;
int total_number_of_ACs_sent_by_master = 0;
int total_number_of_ACs_received_by_slave  = 0;

/** function pointers for CP_FP_decider and FP_termination */

function <bool() > CP_FP_decider;
function <int() > CP_FP_decider_int;

function <bool ( int current_iteration_number, double time_spent ) >  FP_iterate_decider_slave;
function <int ( int reset_count_final, int current_iteration_number, double time_spent ) >  FP_iterate_decider_master;

// pthread_mutex_t lock_center_calculation_tester = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t lock_center_calculation_status = PTHREAD_MUTEX_INITIALIZER;

pthread_mutex_t copy_and_add_shadow_copy = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t complete_and_return_in_distribute_AC = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t recv_AC_mutex= PTHREAD_MUTEX_INITIALIZER;


pthread_mutex_t sr_AC_mutex_send = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t sr_AC_mutex_recv = PTHREAD_MUTEX_INITIALIZER;


pthread_t tth,tth2;
pthread_attr_t att,att2;



pthread_t thread1;
pthread_attr_t attr1;
int global_test_data = 100;

#ifdef PPOF_ne
	ofstream s_ppof;
#endif

vector<My_solution> outside_history_check_solution;
vector<double>      outside_history_check_alpha;

// double alphadist2 = 1.0e-4;
stringstream FP_param_file;
stringstream CP_param_file;
stringstream PWR_param_file;


int MY_SLAVE_ID_from_input = -1;
int WORLD_SIZE_from_input= 1;
string MY_UNIQ_FILENAME = "zzz";


bool using_AC= false;

unsigned reduced_size = history_size_limit_u - history_size_limit_l;
/** when history_size goes over history_size_limit_u reduce it to history_size_limit_l */

bool FP_iterate_use_min_iter_lim ( int current_iteration_number );
bool FP_iterate_use_time_limit( int current_iteration_number, double time );
int time_limit_master_checker(int  reset_count_final);
int min_iterlim_master_checker(int  reset_count_final);



my_signaler signaler;

double get_average_time(Timer &t){
	double s = 0, avg_s = 0;
	IFSLAVE s = t.get_duration();
#ifdef USING_MPI
	MPI_Reduce ( &s, &avg_s, 1, MPI_DOUBLE, MPI_SUM, 0, MPI::COMM_WORLD );
	avg_s /= N_SLAVES;
#else
	avg_s = -1;
#endif
	
	return avg_s;
}

int read_FP_parameters_from_file(cpfp_executable &runner){
	ifstream in (FP_param_file.str().c_str());
	int int_param;
	double dbl_param;
	string str_param;
	int use_what_limit = 0;
	string s = FP_param_file.str();
// 	IFSLAVE_y(1) cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
	
	if (!file_exists(s)){
		IFSLAVE_y(1) cout << "FP parameter file is not available" <<endl;
		return -1;
	}	

	while(!in.eof()){
		string param;
		in >> param;
		
		if(strcmp((param.c_str()),"MAXITER1") == 0) {
			in >> int_param;	
			runner.set_FP_INT_PARAM(ENM_FP_INT_PARAM_MAXITER1, int_param);
		}
		
		if(strcmp((param.c_str()),"MAXITER2") == 0) {
			in >> int_param;	
			runner.set_FP_INT_PARAM(ENM_FP_INT_PARAM_MAXITER2, int_param);
		}
		if(strcmp((param.c_str()),"MAXITERW1") == 0) {
			in >> int_param;				
			runner.set_FP_INT_PARAM(ENM_FP_INT_PARAM_MAXITERW1, int_param);			
		}		
		if(strcmp((param.c_str()),"MAXITERW2") == 0) {
			in >> int_param;	
			runner.set_FP_INT_PARAM(ENM_FP_INT_PARAM_MAXITERW2, int_param);
		}
		if(strcmp((param.c_str()),"MINFLIP") == 0) {
			in >> int_param;	
			runner.set_FP_INT_PARAM(ENM_FP_INT_PARAM_MINFLIP, int_param);
		}
		if(strcmp((param.c_str()),"ROUNDING") == 0) {
			in >> int_param;	
			initial_rounding_type = int_param;
			runner.set_FP_INT_PARAM(ENM_FP_INT_PARAM_ROUNDING, int_param);
		}
		if(strcmp((param.c_str()),"OPTIMIZATION_ALGORITHM") == 0) {
			in >> int_param;	
			runner.set_FP_INT_PARAM(ENM_FP_INT_PARAM_OPTIMIZATION_ALGORITHM, int_param);
		}
		if(strcmp((param.c_str()),"LP_PRESOLVE") == 0) {
			in >> int_param;	
			runner.set_FP_INT_PARAM(ENM_FP_INT_PARAM_LP_PRESOLVE, int_param);
			
		}
		if(strcmp((param.c_str()),"MIP_PRESOLVE") == 0) {
			in >> int_param;	
			runner.set_FP_INT_PARAM(ENM_FP_INT_PARAM_MIP_PRESOLVE, int_param);
		}
		if(strcmp((param.c_str()),"FLIP_TRESHOLD") == 0) {
			in >> dbl_param;	
			runner.set_FP_DOUBLE_PARAM( ENM_FP_DOUBLE_PARAM_FLIP_TRESHOLD, dbl_param);
			
		}
		if(strcmp((param.c_str()),"ALPHA_INIT") == 0) {
			in >> dbl_param;	
			initial_alpha = dbl_param;
			runner.set_FP_DOUBLE_PARAM( ENM_FP_DOUBLE_PARAM_ALPHA_INIT, dbl_param);
			
		}		
		if(strcmp((param.c_str()),"ALPHA_QUOD") == 0) {
			in >> dbl_param;	
			runner.set_FP_DOUBLE_PARAM( ENM_FP_DOUBLE_PARAM_ALPHA_QUOD, dbl_param);
			
		}
		if(strcmp((param.c_str()),"ALPHA_DIST") == 0) {
			in >> dbl_param;	
			runner.set_FP_DOUBLE_PARAM( ENM_FP_DOUBLE_PARAM_ALPHA_DIST, dbl_param);
		}
		if(strcmp((param.c_str()),"USE_FP_LIMIT") == 0) {
			in >> str_param;	
// 			IFSLAVE_y(1) cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
			if (strcmp((str_param.c_str()),"T") == 0){
// 				IFSLAVE_y(1) cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
				use_what_limit = 1;
				FP_iterate_decider_slave = bind ( FP_iterate_use_time_limit,placeholders::_1, placeholders::_2);
				FP_iterate_decider_master = bind(time_limit_master_checker, placeholders::_1);
			}
			if (strcmp((str_param.c_str()),"I") == 0){
// 				IFSLAVE_y(1) cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
				use_what_limit = 2;
				FP_iterate_decider_slave = bind ( FP_iterate_use_min_iter_lim,placeholders::_1 );
				FP_iterate_decider_master = bind(min_iterlim_master_checker, placeholders::_1);
			}
		}
		if(strcmp((param.c_str()),"LIMIT") == 0) {
// 			IFSLAVE_y(1) cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << " use_what_limit: "<<use_what_limit <<  endl;
			if (use_what_limit == 1){
				in >> dbl_param;
				FP_time_limit = dbl_param;
// 				IFSLAVE_y(1) cout << "FP_time_limit: " << FP_time_limit << " dbl_param: " <<dbl_param <<endl;
			}
			if (use_what_limit == 2){
				in >> int_param;
				FP_min_iter_lim = int_param;
			}
		}
		if(strcmp((param.c_str()),"aggressive_flip_probability") == 0) {
			in >> dbl_param;
			aggressive_flip_probability = dbl_param;
		}
	}
	return 0;
}

void read_CP_parameters_from_file(){
// 	ifstream in (CP_param_file.str().c_str());
// 	int int_param;
// 	double dbl_param;
// 	string str_param;
// 	int use_what_limit = 0;
// 	while(!in.eof()){
// 		string param;
// 		in >> param;		
// 	}
 	return;
}

int read_Parallel_parameters_from_file(){
	ifstream in (PWR_param_file.str().c_str());
	int int_param;
	double dbl_param;
	int walk_type = -1;
    
	string s = PWR_param_file.str();
	if (!file_exists(s)){
		IFMASTER {
			cout << "parallel parameter file is not available" <<endl;
			cout << PWR_param_file.str().c_str() <<endl;
		}
		return -1;
	}
	while(!in.eof()){
		string param;
		in >> param;
		if(strcmp((param.c_str()),"slave_generates_own_functions_for_FP") == 0) {
			string str_param;
			in >> str_param;
			if(strcmp((str_param.c_str()),"true") == 0) {
				slave_generates_own_functions_for_FP = true;
			}
			if(strcmp((str_param.c_str()),"false") == 0) {
				slave_generates_own_functions_for_FP = false;
			}
		}
		if(strcmp((param.c_str()),"slave_generates_own_functions_for_CP") == 0) {
			string str_param;
			in >> str_param;
			if(strcmp((str_param.c_str()),"true") == 0) {
				slave_generates_own_functions_for_CP = true;
			}
			if(strcmp((str_param.c_str()),"false") == 0) {
				slave_generates_own_functions_for_CP = false;
			}
		}
		if(strcmp((param.c_str()),"solve_original_centering_problem") == 0) {
			string str_param;
			in >> str_param;
			if(strcmp((str_param.c_str()),"true") == 0) {
				solve_original_centering_problem = true;
			}
			if(strcmp((str_param.c_str()),"false") == 0) {
				solve_original_centering_problem = false;
			}
		}		
		if(strcmp((param.c_str()),"random_walk_type") == 0) {
			string str_param;
			in >> str_param;
			if(strcmp((str_param.c_str()),"short_dikin") == 0) {
				walk_type = 1;
			}
			if(strcmp((str_param.c_str()),"long_dikin") == 0) {
				walk_type = 2;
			}
			if(strcmp((str_param.c_str()),"hit_run") == 0) {
				walk_type = 3;
			}
		}
		if(strcmp((param.c_str()),"number_FP") == 0) {
			string str_param;
			in >> str_param;
			if(strcmp((str_param.c_str()),"original") == 0) {
				in >> int_param;
				type_dist_for_FP.add_instance(OBJECTIVE_TYPE_ORIGINAL, int_param);
			}
			if(strcmp((str_param.c_str()),"random") == 0) {
				in >> int_param;
				if (walk_type == 1)	type_dist_for_FP.add_instance(OBJECTIVE_TYPE_RANDOM_FROM_SHORT_DIKIN, int_param);
				if (walk_type == 2)	type_dist_for_FP.add_instance(OBJECTIVE_TYPE_RANDOM_FROM_LONG_DIKIN, int_param);
				if (walk_type == 3)	type_dist_for_FP.add_instance(OBJECTIVE_TYPE_RANDOM_FROM_HIT_RUN, int_param);
				using_AC = true;
				
			}
			if(strcmp((str_param.c_str()),"perturbed") == 0) {
				in >> int_param;
				if (walk_type == 1)	type_dist_for_FP.add_instance(OBJECTIVE_TYPE_PERTURBED_FROM_SHORT_DIKIN, int_param);
				if (walk_type == 2)	type_dist_for_FP.add_instance(OBJECTIVE_TYPE_PERTURBED_FROM_LONG_DIKIN, int_param);
				if (walk_type == 3)	type_dist_for_FP.add_instance(OBJECTIVE_TYPE_PERTURBED_FROM_HIT_RUN, int_param);
				using_AC = true;
				
			}
		}
		if(strcmp((param.c_str()),"number_CP") == 0) {
			string str_param;
			in >> str_param;
			if(strcmp((str_param.c_str()),"original") == 0) {
				in >> int_param;
				type_dist_for_CP.add_instance(OBJECTIVE_TYPE_ORIGINAL, int_param);
			}
			if(strcmp((str_param.c_str()),"random") == 0) {
				in >> int_param;
				if (walk_type == 1)	type_dist_for_CP.add_instance(OBJECTIVE_TYPE_RANDOM_FROM_SHORT_DIKIN, int_param);
				if (walk_type == 2)	type_dist_for_CP.add_instance(OBJECTIVE_TYPE_RANDOM_FROM_LONG_DIKIN, int_param);
				if (walk_type == 3)	type_dist_for_CP.add_instance(OBJECTIVE_TYPE_RANDOM_FROM_HIT_RUN, int_param);
				using_AC = true;
			}
			if(strcmp((str_param.c_str()),"perturbed") == 0) {
				in >> int_param;
				if (walk_type == 1)	type_dist_for_CP.add_instance(OBJECTIVE_TYPE_PERTURBED_FROM_SHORT_DIKIN, int_param);
				if (walk_type == 2)	type_dist_for_CP.add_instance(OBJECTIVE_TYPE_PERTURBED_FROM_LONG_DIKIN, int_param);
				if (walk_type == 3)	type_dist_for_CP.add_instance(OBJECTIVE_TYPE_PERTURBED_FROM_HIT_RUN, int_param);
				using_AC = true;
				
			}
		}
		if(strcmp((param.c_str()),"perturbation_for_FP") == 0) {
			in >> dbl_param;
			type_dist_for_FP.set_perturbation(dbl_param);
		}
		if(strcmp((param.c_str()),"perturbation_for_CP") == 0) {
			in >> dbl_param;
			type_dist_for_CP.set_perturbation(dbl_param);
		}
		if(strcmp((param.c_str()),"USE_RELATIVE_IMPROVEMENT_IN_FP") == 0) {
			string str_param;
			in >> str_param;
			if(strcmp((str_param.c_str()),"true") == 0) {
				USE_RELATIVE_IMPROVEMENT_IN_FP = true;
			}
			if(strcmp((str_param.c_str()),"false") == 0) {
				USE_RELATIVE_IMPROVEMENT_IN_FP = false;
			}
		}
		if(strcmp((param.c_str()),"USE_RELATIVE_IMPROVEMENT_IN_CP") == 0) {
			string str_param;
			in >> str_param;
			if(strcmp((str_param.c_str()),"true") == 0) {
				USE_RELATIVE_IMPROVEMENT_IN_CP = true;
			}
			if(strcmp((str_param.c_str()),"false") == 0) {
				USE_RELATIVE_IMPROVEMENT_IN_CP = false;
			}
		}
		if(strcmp((param.c_str()),"objective_improvement_coefficient_FP") == 0) {
			in >> dbl_param;
			objective_improvement_coefficient_FP = dbl_param;
		}
		if(strcmp((param.c_str()),"objective_improvement_coefficient_CP") == 0) {
			in >> dbl_param;
			objective_improvement_coefficient_CP = dbl_param;
		}
		if(strcmp((param.c_str()),"objective_improvement_percentage_FP") == 0) {
			in >> dbl_param;
			objective_improvement_percentage_FP = dbl_param;
		}
		if(strcmp((param.c_str()),"objective_improvement_percentage_CP") == 0) {
			in >> dbl_param;
			objective_improvement_percentage_CP = dbl_param;
		}
		if(strcmp((param.c_str()),"global_time_limit") == 0) {
			in >> dbl_param;
			timelimit = dbl_param;
		}
		if(strcmp((param.c_str()),"number_of_iterations") == 0) {
			in >> int_param;
			number_of_iterations = int_param;
		}
		if(strcmp((param.c_str()),"iteration_FP") == 0) {
			in >> int_param;
			FP_iterations = int_param;
			for (int i =0; i < int_param;++i)
				CP_FP_vector.push_back(1);
		}
		if(strcmp((param.c_str()),"iteration_CP") == 0) {
			in >> int_param;
			CP_iterations = int_param;
			for (int i =0; i < int_param;++i)
				CP_FP_vector.push_back(2);
		}
		
	}

	// 	cout << "line: " << __LINE__ <<endl;
	
	
	assert (walk_type>0);
	assert (walk_type<4);
	
	assert(type_dist_for_FP.get_perturbation()>=0);
	assert(type_dist_for_FP.get_perturbation()<=1);
	
	assert(type_dist_for_CP.get_perturbation()>=0);
	assert(type_dist_for_CP.get_perturbation()<=1);
	
	assert(objective_improvement_percentage_CP>=0);
	assert(objective_improvement_percentage_FP>=0);
	
	assert(objective_improvement_coefficient_CP>=0);
	assert(objective_improvement_coefficient_FP>=0);
	assert(timelimit>0);
	type = walk_type;
	return 0;
	
}

void *threadrun_slave_global_time_checker(void* arg){
// 	cout <<  "\t\t\t\t 111111111111111111 " <<endl;
	int trt=0;
	while (1){
		trt++;
		GLOBAL_TIME_CHECKER
		usleep(max_time_limit);
	}
	pthread_exit((void*) NULL);
	return NULL;
}

int external_frac_tester(My_solution s, double &zikkim, int line ){

	int retval =0;
	vector<int> ind =s.get_indices();
	vector<double> ele = s.get_elements();
	zikkim = 0;
	
	for (unsigned i =0; i< ind.size();++i){
		if((external_column_types[ind[i]]) == 0) continue;
		else{
			double r = floor(ele[i] + 0.5) ;
			double frac = fabs(r-ele[i]);
			if (double_inequality(frac,0)){
				zikkim += frac;
				retval++;
			}
		}

	}
	if (retval == 0 || double_equality(zikkim,0)){
// 		cout << "XXXXXXXXXXXXXXXXXXXXXXXXX SLAVE: " << SLAVEID << " line: " << line << " rerval: " << retval << " d_frac :"<< zikkim << endl;
	}
	else{
		cout << "XXXXXXXXXXXXXXXXXXXXXXXXX SLAVE: " << SLAVEID << " line: " << line << " rerval: " << retval << " d_frac :"<< zikkim << endl;
		
	}
	return retval;
}

vector<int> integer_columns_in_slave;

bool history_insert_version1(const My_solution s, double alpha, My_solution &pre){ /**returns true if insert is successful
                                                                            false if already exists    */
    /** quick and dirty*/


	static bool printing = false;
	static int counter=0;
	static int false_returns = 0;
	static int true_returns = 0;
//

    int int_size = (int)outside_history_check_solution.size();


	if(printing) {
		counter++;
		if (counter%100 ==0){
		cout << "SLAVE: " << SLAVEID << " call size: " <<counter << " hist_size: " << int_size << " f: "<< false_returns << " t: "<<true_returns  << " integer_columns_in_slave: "<< integer_columns_in_slave.size() <<endl;
		printing = false;
		}
	}

// 	cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;

    for (int i = int_size-1; i >= 0;--i){

// 		cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << " i: "<< i <<" s: " << int_size << endl;

        if (double_equality(alpha,outside_history_check_alpha[i],alpha_dist)){
// 			cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << " i: "<< i <<" s: " << int_size << endl;

            if (s.equality_on_selected_indices(outside_history_check_solution[i],integer_columns_in_slave)) {
// 				cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << " i: "<< i <<" s: " << int_size << endl;
                pre = outside_history_check_solution[i];
// 				cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
				false_returns++;

                return false;
            }

        }
    }
//     cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;

    outside_history_check_solution.push_back(s);
    outside_history_check_alpha.push_back(alpha);
// 	cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;

    if (outside_history_check_solution.size()>=history_size_limit_u){
// 		cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
        outside_history_check_solution.erase(outside_history_check_solution.begin(), outside_history_check_solution.begin() + reduced_size);
// 		cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
        outside_history_check_alpha.erase(outside_history_check_alpha.begin(), outside_history_check_alpha.begin() + reduced_size);
// 		cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
    }

    true_returns++;

    return true;
}

void history_clear_version1(){
    outside_history_check_solution.clear();
    outside_history_check_alpha.clear();
    return;
}
int rounding_decider(int restarts, int stage = 0){
	if (initial_rounding_type == 8 && restarts <  5) return 0;
	if (initial_rounding_type == 9 && restarts <  5) return 0;
	if (initial_rounding_type == 8 && restarts >= 5) return 6;
	if (initial_rounding_type == 9 && restarts >= 5) return 6;
	return initial_rounding_type;
}

int slave_FP_communicate_with_master(cpfp_executable &arg_runner, bool &arg_will_send,  My_solution &arg_best_incumbent_at_slave_FP, double arg_best_incumbent_objective_at_slave_FP){
//     static int latest_communication_iteration = 0;
		if (running_serial) return 0;

            /** inform master about latest objective if better found from the latest iteration */
        int retvalpp =20;
        if ( arg_will_send ) {
// 			IFSLAVE_y(2) cout << "arg_best_incumbent_objective_at_slave_FP" << arg_best_incumbent_objective_at_slave_FP << endl;
            double latest_objective_sent = arg_best_incumbent_objective_at_slave_FP;
//  			cout << "Slave " << SLAVEID << " SENDing a value of "  << latest_objective_sent << endl;
            #ifdef USING_MPI
	            MPI_Send ( &latest_objective_sent,1,MPI_DOUBLE,0,FOUND_A_FEASIBLE_SOLUTION_IN_FP,MPI_COMM_WORLD );
	            req_stat_vectors r;
                arg_best_incumbent_at_slave_FP._send_mpi(0,TAG_SENDING_SOLUTION,r.v_req,r.v_stat);
	        #else
	            signaler.send_signal_d(0,FOUND_A_FEASIBLE_SOLUTION_IN_FP,latest_objective_sent);
	            arg_best_incumbent_at_slave_FP._send_alt1(0,TAG_SENDING_SOLUTION, signaler);
	        #endif

//                 if(slave_FP_type == 0){
 
//                 }
// 			cout << "Slave " << SLAVEID << " MUST HAVE SENT MPISENDing a value of "  << latest_objective_sent << endl;
				
            retvalpp =10;
            arg_will_send = false;
        }
        /** probe master if it send a new objective value*/

        int flag = -99;
        PROBE(0,TAG_SENDING_OBJECTIVE_VALUE_IN_FP)
        
		if ( flag == TAG_SIGNALER_PROBE_SUCCESSFUL ) {
          //  cout << "in slave "<< SLAVEID << " line: " << __LINE__ << endl;
            double r_d = 1e99;
            RECV_D(0,TAG_SENDING_OBJECTIVE_VALUE_IN_FP)
//             cout << "in slave "<< SLAVEID << " line: " << __LINE__ << endl;
// 			DEBON(r_d, flag)
// 			cout << "Slave " << SLAVEID << " RECEIVING a value of "  << r_d << endl;
			
            if ( r_d < arg_best_incumbent_objective_at_slave_FP ) {
// 				cout << "Slave " << SLAVEID << " USING a value of "  << r_d << " instead of " << arg_best_incumbent_objective_at_slave_FP <<  endl;
				arg_best_incumbent_objective_at_slave_FP = r_d;
				UPDATE_RUNNER_FP(r_d, arg_runner);
				OPTIMALITY_CHECK(arg_runner);
				
				if(!slave_optimality)arg_runner.optimize();
				return retvalpp + 2 ; 
				
            }
            return retvalpp + 1;
        }
    	return retvalpp;

    }

int slave_FP_communicate_with_master_v3(cpfp_executable &arg_runner, int current_iter, My_solution &best_sol, double &best_inc_obj, bool &will_send){
			
	int retval = 1;
	static int latest_communication_iteration = 0;
	int temp = 0;
 		//cout << "slave " << SLAVEID  << " current_iter " << current_iter << " latest_communication_iteration " << latest_communication_iteration <<endl;
	if (((current_iter- latest_communication_iteration) > communicate_every_x_iterations) || (slave_optimality) || (communicate_every_x_iterations <= 1)){
 			//if (will_send) cout << "slave "<< SLAVEID << " sending value " << best_inc_obj << endl;
		temp  = slave_FP_communicate_with_master(arg_runner, will_send, best_sol, best_inc_obj);
		latest_communication_iteration = current_iter;
 			//if (will_send) cout << "slave "<< SLAVEID << " MUST HAVE SENT " << best_inc_obj << endl;
	}
	
 	//cout << "line_ " << __LINE__ << " slave_opt:" << slave_optimality << endl;
	OPTIMALITY_CHECK(arg_runner)

 	//cout << "line_ " << __LINE__ << " slave_opt:" << slave_optimality << endl;
	
	if((temp  == 2)) {
		return 0;
	}
	/** need to generate a new starting solution if not original objective*/
	return retval;
}

#ifdef USING_IOPTIMIZE
bool AC_compare(shadow_copy_ioptimize_wrapper lhs, shadow_copy_ioptimize_wrapper rhs){
	return (lhs.get_number_of_times_AC_calculated_or_received() < rhs.get_number_of_times_AC_calculated_or_received());
}
#endif
/**master sends and slaves recv*/ 
/**this function only sends latest_AC_copy_send_by_master or receives latest_AC_copy_recv_by_slave it does not generate/replace latest_AC_copy_send_by_master and it does not create/reset a ioptimize_wrapper from a shadow copy*/

int AC_comm_v2(){
	// 	static int nt_sr_AC = 0;
	// 	nt_sr_AC++;
	static int latest_AC_sent = -1;
	static int latest_AC_received = -1;
	
	IFSLAVE{
		
		int latest_recvd = -1;
		// 		cout <<" latest_AC_received : " << latest_AC_received <<  endl;
		int flag = -99;
		
        PROBE(0,TAG_SENDING_ANALYTIC_CENTER)
        
		if (flag == TAG_SIGNALER_PROBE_SUCCESSFUL){
            #ifdef USING_IOPTIMIZE
			req_stat_vectors r;
			shadow_copy_ioptimize_wrapper to_receive;
			to_receive._recv_mpi ( 0, TAG_SENDING_ANALYTIC_CENTER,r.v_req, r.v_stat );
			pthread_mutex_lock(&copy_and_add_shadow_copy);
			
			if (to_receive.get_number_of_times_AC_calculated_or_received() > latest_AC_received){
				latest_AC_copy_recv_by_slave = to_receive;
				latest_AC_received = latest_AC_copy_recv_by_slave.get_number_of_times_AC_calculated_or_received();
				total_number_of_ACs_received_by_slave++;
				pthread_mutex_unlock(&copy_and_add_shadow_copy);
				return latest_AC_received;
			}
			else{
				pthread_mutex_unlock(&copy_and_add_shadow_copy);
			}
            #endif
			
		}
		return 0;
	}
	
	IFMASTER{
		#ifdef USING_IOPTIMIZE
        pthread_mutex_lock(&copy_and_add_shadow_copy);
        
		int latest = latest_AC_copy_by_master.get_number_of_times_AC_calculated_or_received();
// 		pthread_mutex_unlock(&copy_and_add_shadow_copy);
		if (latest > latest_AC_sent){
// 			pthread_mutex_lock(&copy_and_add_shadow_copy);			
			shadow_copy_ioptimize_wrapper to_send = latest_AC_copy_by_master;
			latest_AC_sent = latest;
			pthread_mutex_unlock(&copy_and_add_shadow_copy);
			for ( int i = 1; i < WORLDSIZE; ++i ) {
				req_stat_vectors r;
				to_send._send_mpi ( i, TAG_SENDING_ANALYTIC_CENTER,r.v_req, r.v_stat );
			}
			return latest_AC_sent;
		}
		else {
			pthread_mutex_unlock(&copy_and_add_shadow_copy);
			return 0;
		}
#endif
	}
	
	return -2;
}


/** threadrun AC calculation */
void *AC_calc_v2_thread(void *arg){
	
	
	#ifdef USING_IOPTIMIZE
	shadow_copy_ioptimize_wrapper new_AC;
	
	if ( !master_iop.get_problem_initialized() ) {
		
		master_iop.set_seed ( my_default_seed );
		master_iop.setCenProbType ( solve_original_centering_problem );
		master_iop.initialize2 ( global_filename );
		if ( type > 3 ) {
			master_iop.set_walk_type ( type-3 );
		} else {
			master_iop.set_walk_type ( type );
		}
	}
	master_iop.update_sides ( master_lb_for_iop, master_ub_for_iop );
	if ( master_iop.get_analytic_center_changed() ) {
		time_to_calculate_AC.restart();
		master_iop.calculate_analytic_center_self(master_iop_n_threads);
		new_AC = shadow_copy_ioptimize_wrapper(master_iop);
		time_to_calculate_AC.pause();
		
	}
	
	/** be careful locking something else within a lock*/
	pthread_mutex_lock(&copy_and_add_shadow_copy);
	// 			shadow_copy_of_all_calculated_ACs.push_back(new_AC);
	latest_AC_copy_by_master = new_AC;
	total_number_of_ACs_computed_in_master++;
	pthread_mutex_unlock(&copy_and_add_shadow_copy);
	
	pthread_mutex_unlock(lock_center_calculation_status);
	
	pthread_exit ( ( void* ) NULL );
    #endif
	return NULL;
	
}


/** part of this needs to be done in a thread */

int AC_calc_v2(bool complete_and_return = false){
	if(!using_AC){
		return 0;
	}
	bool calculation_starts = false;
	if (complete_and_return){
		pthread_mutex_lock(&lock_center_calculation_status);
		calculation_starts = true;
	}
	else {
		if (pthread_mutex_trylock(&lock_center_calculation_status) == 0){
			/** trylock successful */ 
			calculation_starts = true;
		}
	}
	
	
	if (calculation_starts){
		/** calculate AC in a thread this part in a thread */
		pthread_create ( &thread1,&attr1,AC_calc_v2_thread, ( void* ) ( NULL ) );
		return 2;
	}
	if (complete_and_return){
		/** this makes sure that we wait completion of center calculation*/
		pthread_mutex_lock(&lock_center_calculation_status);
		pthread_mutex_unlock(&lock_center_calculation_status);
		return 1;
	}
	
	return 0;
	
	
}


void *threadrun_receive_analytic_center_in_slave ( void* arg ){ 
	static int number_of_times_AC_receive_called = 0;
	number_of_times_AC_receive_called++;
    
// 	static int latest_AC_received = -1;
	int t;
// // 	pthread_mutex_trylock(recv_AC_mutex);
// // 	pthread_mutex_trylock ( &recv_AC_mutex ) == 0 
// // 	if (pthread_mutex_trylock ( &recv_AC_mutex ) == 0) 
// 		
// 	while(!pthread_mutex_trylock ( &recv_AC_mutex )){
// 		pthread_mutex_unlock ( &recv_AC_mutex );
// 		shadow_copy_ioptimize_wrapper to_receive;
// 		/** receive AC if necc */
// 		int flag = -99;
// 		MPI_Status status;
// // 		cout <<" latest_AC_received : " << latest_AC_received <<  endl;
// 		t = MPI_Iprobe ( 0,TAG_SENDING_ANALYTIC_CENTER,MPI::COMM_WORLD, &flag, &status );
// 		req_stat_vectors r;
// 		if ( flag == 1 ) {
// 			to_receive._recv_mpi ( 0, TAG_SENDING_ANALYTIC_CENTER,r.v_req, r.v_stat );
// 			pthread_mutex_lock(&copy_and_add_shadow_copy);
// 			shadow_copy_of_all_received_ACs.push_back(to_receive);
// 			sort(shadow_copy_of_all_received_ACs.begin(), shadow_copy_of_all_received_ACs.end(),AC_compare);
// 			pthread_mutex_unlock(&copy_and_add_shadow_copy);
// 			usleep(SLEEP_TIME_BETWEEN_AC_PROBES);
// 		}
// 		else{
// 			usleep(SLEEP_TIME_BETWEEN_AC_PROBES*10);
// 		}
// 	}
// 	
// 	
// 	pthread_exit ( ( void* ) NULL );
	return NULL;
	t++;
		
	
}

void *threadrun_distribute_analytic_center_in_master ( void* arg ){

    static int number_of_times_AC_distribute_called = 0;
    number_of_times_AC_distribute_called ++;


    static int latest_AC_sent = -1;

//     //cout << "threadrun_distribute_analytic_center_in_master line: "<<  __LINE__  << " distribute complete total_number_of_ACs_sent_by_master: " << total_number_of_ACs_sent_by_master << " number_of_times_AC_distribute_called: " << number_of_times_AC_distribute_called << " shadow_cp_size: "<< shadow_copy_of_all_calculated_ACs.size() << endl;
// 
//     pthread_mutex_lock(&copy_and_add_shadow_copy);
//     
// //     int latest = (int)shadow_copy_of_all_calculated_ACs.size()-1;
// // 	shadow_copy_ioptimize_wrapper to_send = shadow_copy_of_all_calculated_ACs[latest];
// //     
//     pthread_mutex_unlock(&copy_and_add_shadow_copy);
// //     bool sending =false;
//     if (latest > latest_AC_sent){
//         /** new AC needs to be sent*/
// 
// //         sending =true;
//         latest_AC_sent = latest;
//         for ( int i = 1; i < MPI::COMM_WORLD.Get_size(); ++i ) {
//             //         cout << "threadrun_distribute_analytic_center_in_master line: "<<  __LINE__  << " sending to "<< i << endl;
// 
//             req_stat_vectors r;
//             to_send._send_mpi ( i, TAG_SENDING_ANALYTIC_CENTER,r.v_req, r.v_stat );
//             //         cout << "threadrun_distribute_analytic_center_in_master line: "<<  __LINE__  << " sent to "<< i << endl;
// 
//         }
//         total_number_of_ACs_sent_by_master++;
//     }
// 
// 
//     //cout << "threadrun_distribute_analytic_center_in_master line: "<<  __LINE__  << " distribute complete total_number_of_ACs_sent_by_master: " << total_number_of_ACs_sent_by_master << " number_of_times_AC_distribute_called: " << number_of_times_AC_distribute_called << " latest AC sent: " << latest_AC_sent << " AC_number: " << shadow_copy_of_all_calculated_ACs[latest].get_number_of_times_AC_calculated_or_received();
//    // if (sending) cout << " SENT NOW" << endl;
//    // if (!sending) cout << " DID NOT SEND NOW" << endl;
//     pthread_mutex_unlock(&complete_and_return_in_distribute_AC);
//     pthread_exit ( ( void* ) NULL );
    return NULL;

}

void *threadrun_calculate_analytic_center_in_master ( void* arg ){
    if(!using_AC){
        pthread_exit ( ( void* ) NULL );
        return NULL;
    }
    static int called_times = 0;
    called_times ++;
    /** since the thread is locked we are sure that analytic center will be calculated */
//     cout << " threadrun_calculate_analytic_center_in_master called  " << called_times << " total_number_of_ACs_computed_in_master " << total_number_of_ACs_computed_in_master << endl;

#ifdef USING_IOPTIMIZE
//     master_iop_AC_needs_to_calculated = 0;

    if ( !master_iop.get_problem_initialized() ) {

        master_iop.set_seed ( my_default_seed );
        master_iop.setCenProbType ( solve_original_centering_problem );
        master_iop.initialize2 ( global_filename );
        if ( type > 3 ) {
            master_iop.set_walk_type ( type-3 );
        } else {
            master_iop.set_walk_type ( type );
        }
    }
//     cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
	
    master_iop.update_sides ( master_lb_for_iop, master_ub_for_iop );
	
// 	cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
    bool recalculation_occured = false;
    if ( master_iop.get_analytic_center_changed() ) {
// 		cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
        recalculation_occured = true;
        time_to_calculate_AC.restart();
// 		cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
        master_iop.calculate_analytic_center_self(master_iop_n_threads);
// 		cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
        shadow_copy_ioptimize_wrapper new_AC(master_iop);
		time_to_calculate_AC.pause();
		
		// 		cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;

        pthread_mutex_lock(&copy_and_add_shadow_copy);
        shadow_copy_of_all_calculated_ACs.push_back(new_AC);
		latest_AC_copy_by_master = new_AC;
// 		cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
        total_number_of_ACs_computed_in_master++;
//         cout << " shadow_cp_size: "<< shadow_copy_of_all_calculated_ACs.size() << endl;
// 		cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
        pthread_mutex_unlock(&copy_and_add_shadow_copy);
// 		cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;

    }
//     cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
	pthread_mutex_unlock ( &lock_center_calculation_status );
	
//     master_iop_AC_needs_to_calculated = -1;
// 	cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
    if (!recalculation_occured){
//         cout << " threadrun_calculate_analytic_center_in_master called  " << called_times << " total_number_of_ACs_computed_in_master " << total_number_of_ACs_computed_in_master << endl;
    }
    else{
//        cout << " threadrun_calculate_analytic_center_in_master called  " << called_times << " total_number_of_ACs_computed_in_master " << total_number_of_ACs_computed_in_master << endl;
        threadrun_distribute_analytic_center_in_master(arg);
		
    }
//     cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
   
	
#endif
    pthread_exit ( ( void* ) NULL );
    return NULL;
}

/**  if this function is called, that means there is need to calculate AC*/
int master_calculate_AC_if_possible()  /** returns 0, -1  if failed to calculate*/{
    /** try to lock */
    int retval =0;
    if ( pthread_mutex_trylock ( &lock_center_calculation_status ) == 0 ) {
        /** success means that ac calculation will stat in a new thread */

        pthread_create ( &thread1,&attr1,threadrun_calculate_analytic_center_in_master, ( void* ) ( NULL ) );
        retval = total_number_of_ACs_computed_in_master;
        return retval;

    } else {
        /** failed to lock. Previous analytic center calculation is continuing or AC is being distributed */
        return 0;
    }
    return -1;
}

/** input what_to_do = 	0, calculate_analytic_center
 * 						1, create_objective_function_list NOT IMPLEMENTED use 3 and 4 instead
 * 						2, disrtribute_AC
 * 						3, just lock for me (must have complete_and_return == true)
 * 						4, just unlock for me (must have complete_and_return == true)
 *                      5, update_sides;
 *
 *
 * if complete_and_return == false, it it tries to run in a thread. If thread creation is successful, returns 0, else returns the currnet master_iop_status.
 * else (complete_and_return == true) wait until locks are cleared run (possibly in this thread) and return the result. */


int master_all_iop_related_function_calls ( int what_to_do, bool complete_and_return = false ){

    if (!using_AC) return -9;

//     static int master_iop_status = -2;
/**    master_iop_AC_needs_to_calculated
                                            * -99, not even initialized
											*	-1, there is no need to calculate
											*	 0, being calculated right now
											*	 1, needs to be calculated
											*	 2, is being distributed
											*/

// 	static int master_iop_status_secdond = 0;
//     cout << __LINE__  << " master_all_iop_related_function_calls( " << what_to_do << ", " << complete_and_return << ")" << endl;

    switch ( what_to_do ) {
    case 0:
//         if ( (master_iop_status == -1) && (master_iop_status_second=0)) {
//             cout << "master_all_iop_related_function_calls line: "<<  __LINE__  << " returns " << -1 << endl;
//             return -1;
//         }

        if ( complete_and_return ) {

            /** wait until we can lock*/

            pthread_mutex_lock ( &lock_center_calculation_status );


//             master_iop_status = 0;
            pthread_create ( &thread1, &attr1, threadrun_calculate_analytic_center_in_master, ( void* ) NULL );
            /** */
            pthread_mutex_lock ( &lock_center_calculation_status );
//             if (master_iop_status ==0) master_iop_status = -1;
            pthread_mutex_unlock ( &lock_center_calculation_status );
//             cout << "master_all_iop_related_function_calls line: "<<  __LINE__  << " returns " << 0 << endl;
            return 0;

        } 
        else {
			if ( pthread_mutex_trylock ( &lock_center_calculation_status ) == 0 ) {
//             master_iop_status = 0;
			master_iop_n_threads = 1;
// 			cout << "master_all_iop_related_function_calls line: "<<  __LINE__  << " Will return " << 0 << endl;
// 			cout.flush();
            pthread_create ( &thread1, &attr1, threadrun_calculate_analytic_center_in_master, ( void* ) NULL );
//             cout << "master_all_iop_related_function_calls line: "<<  __LINE__  << " returns " << 0 << endl;
// 			cout.flush();
			return 0;
        	}
        	else{
// 				cout << " failed to lock:" <<endl;
        	}
		}
        
        break;
// 		case 1:
// 			/** since we need to create_objective_functionlist we may have to wait until AC calculation is done*/
// 			if (complete_and_return){
// 				pthread_mutex_lock(lock_center_calculation_status);
// 				pthread_mutex_lock(lock_center_calculation_status);
// 				pthread_mutex_unlock(lock_center_calculation_status);
// 			}
// 			else if (pthread_mutex_trylock(lock_center_calculation_status) == 0){
// 			}
// 			break;
    case 2:
        if ( complete_and_return ) {
            pthread_mutex_lock ( &complete_and_return_in_distribute_AC);
//             master_iop_status = 2;
            pthread_create ( &thread1, &attr1, threadrun_distribute_analytic_center_in_master, ( void* ) NULL );

            /** we must have locked the mutex, need to wait until all cleared to return*/
            pthread_mutex_lock ( &complete_and_return_in_distribute_AC );

			pthread_mutex_unlock ( &complete_and_return_in_distribute_AC );
        } else {
			if ( pthread_mutex_trylock ( &complete_and_return_in_distribute_AC ) == 0 ) {
				//             master_iop_status = 0;
				master_iop_n_threads = 1;
// 				cout << "master_all_iop_related_function_calls line: "<<  __LINE__  << " Will return " << 0 << endl;
// 				cout.flush();
				pthread_create ( &thread1, &attr1, threadrun_distribute_analytic_center_in_master, ( void* ) NULL );
				
// 				cout << "master_all_iop_related_function_calls line: "<<  __LINE__  << " returns " << 0 << endl;
// 				cout.flush();
				return 0;
			}
			else{
// 				cout << " failed to lock:" <<endl;
			}
			
			
            return 2;
        }
        break;
    case 3:
        if ( complete_and_return ) {
            pthread_mutex_lock ( &lock_center_calculation_status );
            return 0;
        }
        break;
    case 4:
        if ( complete_and_return ) {
            pthread_mutex_unlock ( &lock_center_calculation_status );
//             cout << "master_all_iop_related_function_calls line: "<<  __LINE__  << " returns " << 0 << endl;
            return 0;
        }
        break;
    case 5:
//         master_iop.update_sides(master_lb_for_iop, master_ub_for_iop);
        break;
    }
//     cout << "master_all_iop_related_function_calls line: "<<  __LINE__  << " returns " << master_iop_status << endl;
    return 0;

}

int update_slave_iop (bool reset = false){
	 #ifdef USING_IOPTIMIZE 
// 	cout << " shadow_copy_of_all_received_ACs.size(): "	<< shadow_copy_of_all_received_ACs.size() << endl;
	
	
	if (!reset){
		if (slave_iop.get_number_of_times_AC_calculated_or_received() < (int)shadow_copy_of_all_received_ACs.size() ){
// 			cout << " shadow_copy_of_all_received_ACs.size(): "	<< shadow_copy_of_all_received_ACs.size() << endl;
			return -1;
		}
	}
	/** receive AC if necc */
	int flag = -99;
	MPI_Status status;
// 	cout <<" latest_AC_received : " << latest_AC_received <<  endl;
	int t = MPI_Iprobe ( 0,TAG_SENDING_ANALYTIC_CENTER,MPI::COMM_WORLD, &flag, &status );
	req_stat_vectors r;
	if ( flag == TAG_SIGNALER_NEW_SIGNAL) {
		shadow_copy_ioptimize_wrapper to_receive;
		
		to_receive._recv_mpi ( 0, TAG_SENDING_ANALYTIC_CENTER,r.v_req, r.v_stat );

		pthread_mutex_lock(&copy_and_add_shadow_copy);
		shadow_copy_of_all_received_ACs.push_back(to_receive);
		sort(shadow_copy_of_all_received_ACs.begin(), shadow_copy_of_all_received_ACs.end(),AC_compare);
		pthread_mutex_unlock(&copy_and_add_shadow_copy);
		++t;
	}
	
	shadow_copy_ioptimize_wrapper to_recv;
	int latest_update_index = -1;
	int latest_version = -9;
	
	pthread_mutex_lock(&copy_and_add_shadow_copy);
	for (unsigned i = 0; i < shadow_copy_of_all_received_ACs.size();++i){
		
		if (latest_version < shadow_copy_of_all_received_ACs[i].get_number_of_times_AC_calculated_or_received()){
			latest_version = shadow_copy_of_all_received_ACs[i].get_number_of_times_AC_calculated_or_received();
			latest_update_index = i;
		}
	}
	
	if (latest_version>0 && latest_update_index >=0 ){
		assert (latest_update_index < (int)shadow_copy_of_all_received_ACs.size());
		to_recv = shadow_copy_of_all_received_ACs[latest_update_index];
	}
	
	pthread_mutex_unlock(&copy_and_add_shadow_copy);
	
	if (latest_version > 0) {
		if ((slave_iop.get_number_of_times_AC_calculated_or_received() < latest_version ) || reset ){
			/** a newer version is available*/
			slave_iop.set_from_shadow_copy(to_recv);	
			slave_iop.iop_set_seed ( slave_seed );
			
// 			slave_iop.iop_set_seed ( 2 );
			
			// 		usleep(10000000);
			
			// 		cout << "init sampler: " <<
			slave_iop.init_sampler();
			return 1;
		}
		
	}
#endif
	
	return 0;
	
}

int slave_create_objective_general ( cpfp_executable &runner, My_solution &arg_current_lp_optimum_at_slave, bool it_is_CP ){

    static int local_counter =0;
    static int received_AC_counter = 0;
    int obj_type = 0;
    double pert = 0;

    local_counter++;
    if ( it_is_CP ) {
        obj_type = type_dist_for_CP.I_am_running_what_objective ( SLAVEID );
        pert = type_dist_for_CP.get_perturbation();
    } else {
        obj_type = type_dist_for_FP.I_am_running_what_objective ( SLAVEID );
        pert = type_dist_for_FP.get_perturbation();

    }
//     int retval = -1;
    My_objective_function new_obj = runner.get_original_objective_function();
	/** NOTE: THIS NEEDS TO CHANGE*/
	IFMASTER{
		runner.set_objective_function ( new_obj, 1 );
		return 0;
	}
    if ( obj_type==0 ) {
		if (it_is_CP)  runner.set_objective_function ( new_obj, 3 );
		else  runner.set_objective_function ( new_obj, 1 );
        /** do not return here or else other slaves cannot get theirs ACs*/
//         retval = 0;
//         return 0;
    }
    /** initialize if necc */
    #ifdef USING_IOPTIMIZE 
    if ( !(slave_iop.get_problem_initialized()) ) {
        if (slave_iop.get_analytic_center_changed() && (received_AC_counter>0)){cout << " line " << __LINE__ <<  "-- WHAT THE HACK " << endl; }
        slave_iop.set_seed ( my_default_seed );
        slave_iop.setCenProbType ( solve_original_centering_problem );
        slave_iop.initialize2 ( global_filename );
        slave_seed = my_seed_generator.lth_seed_lcm ( SLAVEID+1 );
        slave_iop.set_seed(slave_seed);
        if ( type > 3 ) {
            slave_iop.set_walk_type ( type-3 );
        } else {
            slave_iop.set_walk_type ( type );
        }
        if (slave_iop.get_analytic_center_changed() && (received_AC_counter>0)){cout << " line " << __LINE__ <<  "-- WHAT THE HACK " << endl; }
    }
    /** receive AC if necc */
	
	
// 	IFSLAVE{
// 		
// 		
// 		int flag = -99;
// 		MPI_Status status;
// 
// 		int t = MPI_Iprobe ( 0,TAG_SENDING_ANALYTIC_CENTER,MPI::COMM_WORLD, &flag, &status );
// 		req_stat_vectors r;
// 		if ( flag == 1 ) {
// 	//         cout << " slave " << SLAVEID << " receiving AC at line " << __LINE__ << endl;
// 			if (slave_iop.get_analytic_center_changed() && (received_AC_counter>0)){cout << " line " << __LINE__ <<  "-- WHAT THE HACK " << endl; }
// 			slave_iop._recv_mpi ( 0, TAG_SENDING_ANALYTIC_CENTER,r.v_req, r.v_stat );
// 	//         cout << " slave " << SLAVEID << " received AC at line " << __LINE__ << " AC_number: "<< slave_iop.get_number_of_times_AC_calculated_or_received() << endl;
// 			slave_iop.after_receive_update();
// 			slave_iop.iop_set_seed ( slave_seed );
// 			slave_iop.init_sampler();
// 			++t;
// 			received_AC_counter++;
// 			if (slave_iop.get_analytic_center_changed() && (received_AC_counter>0)){cout << " line " << __LINE__ <<  "-- WHAT THE HACK " << endl; }
// 		}
// 	}


	 if (update_slave_iop() == 0){ /*cout << "slave_iop for slave "<<SLAVEID <<" not changed " <<endl; */}
// 	 cout << " update_slave_iop(): returns " << update_slave_iop() <<endl;
//     if (slave_iop.get_analytic_center_changed() && (received_AC_counter>0)){cout << " line " << __LINE__ <<  "-- WHAT THE HACK " << endl; }

    vector<double> current_center_as_v = slave_iop.get_current_center_point();

    vector<double> slave_random_point;

    /** create random points */
	received_AC_counter = slave_iop.get_number_of_times_AC_calculated_or_received();
	
	//     bool rand_point_generated = false;
    if((received_AC_counter>0) && (obj_type != OBJECTIVE_TYPE_ORIGINAL)){

        for (int i =0; i < 10;++i) {
//             time_to_create_RW_points.restart();
			time_to_create_RW_points_slave.restart();
          
            int re = slave_iop.get_next_random_point ( slave_random_point );

			time_to_create_RW_points_slave.pause();
			
            if (re  != 0){
                cout << " NO SAMPLE POINT GENERATED attempt "<< i+1 << " re: " << re << endl;
				if (update_slave_iop(true) == 0) { cout << "slave_iop NEEDS TO RESET for slave "<<SLAVEID << " but not changed " <<endl; }
				
// 				slave_iop.iop_set_seed(slave_seed);
// 				slave_iop.init_sampler();
                continue;
            }
            else
            {
//                 rand_point_generated= true;

                i =10;
                break;
            }
        }

    }
#endif
    else if ((obj_type != OBJECTIVE_TYPE_ORIGINAL)){
//         cout << " SLAVE " << SLAVEID << " line " << __LINE__ << " slave_seed " <<slave_seed << endl;
        return 1;
    }

    vector<double> current_lp_optimum_as_v;

    switch ( obj_type ) {
    case OBJECTIVE_TYPE_ORIGINAL:
        new_obj = runner.get_original_objective_function();
        break;
    case OBJECTIVE_TYPE_RANDOM_FROM_SHORT_DIKIN:
    case OBJECTIVE_TYPE_RANDOM_FROM_LONG_DIKIN:
    case OBJECTIVE_TYPE_RANDOM_FROM_HIT_RUN:

#ifdef USING_IOPTIMIZE
        new_obj = My_objective_function ( vector_subtract<double> ( slave_random_point, current_center_as_v ) );
#endif
        break;
    case OBJECTIVE_TYPE_PERTURBED_FROM_SHORT_DIKIN:
    case OBJECTIVE_TYPE_PERTURBED_FROM_LONG_DIKIN:
    case OBJECTIVE_TYPE_PERTURBED_FROM_HIT_RUN:

#ifdef USING_IOPTIMIZE
        current_lp_optimum_as_v = arg_current_lp_optimum_at_slave.get_solution_vector();
        new_obj= My_objective_function (
                     vector_subtract<double> (
                         vector_multiply_sum<double> (
                             pert, slave_random_point,
                             1-pert, current_lp_optimum_as_v
                         ),
                         current_center_as_v ) );

#endif

        /** a new objective function such that
         * d_p   = perturbation * d_r + (1-perturbation) * d_o
         *       = perturbation * (x_r - x_ac) + (1-perturbation) * (x*_lp - x_ac)
         *       = perturbation * x_r + (1-perturbation) * x*_lp - x_ac
         *       = {perturbation * x_r + (1-perturbation) * x*_lp} - x_ac
         *       = vector_multiply_sum{p,x_r,(1-p),x_lp} - x_ac
         *       = vector_subtract{ vector_multiply_sum{p,x_r,(1-p),x_lp}, x_ac}
         */
        break;

    default:
        new_obj = runner.get_original_objective_function();
        break;
    }

    if (new_obj.is_empty()) {
//         new_obj = runner.get_auxilary_objective_function();
        cout << " SLAVE " << SLAVEID << " line " << __LINE__ << " new_obj is empty " << slave_seed << endl;
		new_obj = runner.get_original_objective_function();
		if (new_obj.is_empty()) return -1; 
		
    }
    else{
//         stringstream asd;
//         asd << "to_be_deleted" << SLAVEID << ".txt";
//         ofstream zikkim(asd.str().c_str(), ios_base::app);
//         zikkim << "slave " << SLAVEID << " iteration " << local_counter << " ";
//         new_obj.shortlineprint(zikkim );
    }
    //new_obj.print();
    if (it_is_CP)  runner.set_objective_function ( new_obj, 3 );
	if (!it_is_CP)  runner.set_objective_function ( new_obj, 1 );
	
    return 0;

}

int slave_create_objective_FP ( cpfp_executable &runner, My_solution &arg_current_lp_optimum_at_slave ){
    return slave_create_objective_general ( runner, arg_current_lp_optimum_at_slave,false );
}

int slave_create_objective_CP ( cpfp_executable &runner, My_solution &arg_current_lp_optimum_at_slave ){
    return slave_create_objective_general ( runner, arg_current_lp_optimum_at_slave,true );
}

bool CP_FP_three_rounds_of_FP_than_a_single_CP(){
    static int call_count = 0;
    call_count++;
    if ( call_count%4 == 0 ) {
        return true;
    }
    return false;
}

int CP_FP_n_rounds_of_FP_than_m_CP(){
	static int call_count = 0;
	call_count++;
	if (call_count >= (int)CP_FP_vector.size()) call_count = 0;
	return CP_FP_vector[call_count];

}


bool FP_iterate_use_min_iter_lim ( int current_iteration_number ){
    static bool sent = false;
    static bool asd = true;
    static int p_time = 1;
    cout << "FP_iterate_use_min_iter_lim enter at line " << __LINE__ << endl; 
    if (current_iteration_number <0 || time < 0){
        sent = false;
    }
//     if (current_iteration_number %100 == 0)     {
// 	    cout << " SLAVE " << SLAVEID << " iteration " << current_iteration_number << endl;
// 	}

// 	IFSLAVE_y(6)    cout << " SLAVE " << SLAVEID << " iteration " << current_iteration_number << endl;
    if (asd){
        time_temp.start();
    }

    if (time_temp.stop() > p_time){
        cout << " SLAVE " << SLAVEID << " current_iteration_number: " << current_iteration_number << endl;
        p_time+=1;
    }

    if ( current_iteration_number < FP_min_iter_lim && !slave_optimality) {

        sent = false;
        cout << "FP_iterate_use_min_iter_lim returns " << true << " at line " << __LINE__ << endl; 
        return true;
    }
	if (running_serial) return true;
//     if (current_iteration_number %100 == 0)  cout << " SLAVE " << SLAVEID << " iteration " << current_iteration_number << endl;

    if ( !sent ) {
// 		cout << " SLAVE " << SLAVEID <<" SENDING I AM DONE " << endl;
   	#ifdef USING_MPI
    	MPI_Send ( &current_iteration_number, 1, MPI_INT, 0,  TAG_SENDING_DONE_MINIMUM_NUMBER_OF_ITERATIONS, MPI::COMM_WORLD );
       #else
 		signaler.send_signal_i(0,TAG_SENDING_DONE_MINIMUM_NUMBER_OF_ITERATIONS,current_iteration_number);
       #endif
        
        sent = true;
    }

    int flag = -99;
    
    PROBE(0,TAG_ALL_DONE_MINIMUM_NUMBER_OF_ITERATIONS)

	if ( flag == TAG_SIGNALER_PROBE_SUCCESSFUL ) {
        int n_done;
        int r_i;
        RECV_I(0,TAG_ALL_DONE_MINIMUM_NUMBER_OF_ITERATIONS)
        n_done = r_i;
//information  to clear the recv buffer for the next iteration, if we do not receive next IPROBE call will result TRUE */
//         cout << " SLAVE " << SLAVEID <<" RECEIVED ALL DONE sending final number of iterations: " << current_iteration_number << endl;
        int s_i = current_iteration_number;
        SEND_I(0,TAG_SENDING_ALL_DONE_FINAL_NUMBER_OF_ITERATIONS)
        cout << " SLAVE " << SLAVEID << " sent " << endl;
        asd = true;
        cout << "FP_iterate_use_min_iter_lim returns " << false << " at line " << __LINE__ << endl; 
        return false;

    }
    cout << "FP_iterate_use_min_iter_lim returns " << true << " at line " << __LINE__ << endl; 
    return true;

}

bool  FP_iterate_use_time_limit(int current_iteration_number, double time){
        static bool sent = false;
// 		static bool check = false;
	
    //DEBON("FP_iterate_use_time_limit enter at line ") 
    if (current_iteration_number <0 || time < 0){
        sent = false;
// 		if (check) cout << " SLAVE " << SLAVEID << " line:  " << __LINE__ << " sent = false again time" << time <<" iter: "<< current_iteration_number <<  endl;
		
    }

//     static double last_print_time = 0;
//     if (time > last_print_time + 0.99 ){
//         cout << "slave " << SLAVEID << " FP_ time " << time << endl;
//         last_print_time = time;
//     }
    if (time < FP_time_limit && !slave_optimality){
        sent = false;
		//DEBON("time:", time , " FP_time_limit: ", FP_time_limit);
// 		if (check) cout << " SLAVE " << SLAVEID << " line:  " << __LINE__ << " sent = false again" << endl;
        //DEBON("FP_iterate_use_time_limit returns true",0); 

        return true;
    }
	if(0){//running_serial) {
        DEBOFF ("FP_iterate_use_time_limit returns false",0); 
        return false;
    }
    if (!sent){
        int s_i = current_iteration_number;
        //SEND_I(0,TAG_SENDING_DONE_WITH_TIME_LIMIT)
    	 #ifdef USING_MPI
    	MPI_Send ( &current_iteration_number, 1, MPI_INT, 0,  TAG_SENDING_DONE_WITH_TIME_LIMIT, MPI::COMM_WORLD );
       #else
		 signaler.send_signal_i(0,TAG_SENDING_DONE_WITH_TIME_LIMIT,current_iteration_number);
       #endif
        
// 		check = true;
        sent = true;
		cout << " SLAVE " << SLAVEID << " informs master that it is done by time  " << endl;
        
	}
    int flag = -99;
    PROBE(0,TAG_ALL_DONE_WITH_TIME_LIMIT)
    
	if ( flag == TAG_SIGNALER_PROBE_SUCCESSFUL ) {
//         cout << " SLAVE " << SLAVEID <<" RECEIVED ALL DONE sending final number of iterations: " << current_iteration_number << endl;

        int r_i;
        RECV_I(0,TAG_ALL_DONE_WITH_TIME_LIMIT);
        
        //MPI_Recv ( &n_done,1,MPI_INT,0,TAG_ALL_DONE_WITH_TIME_LIMIT,MPI::COMM_WORLD,&status ); 
        
        /** altough not really needed we are receiving this  information  to clear the recv buffer for the next iteration, if we do not receive next IPROBE call will result TRUE */

//         cout << " SLAVE " << SLAVEID <<" RECEIVED ALL DONE sending final number of iterations: " << current_iteration_number << " n_done:" <<n_done << endl;
        
    
 #ifdef USING_MPI
        MPI_Send ( &current_iteration_number, 1, MPI_INT, 0,  TAG_SENDING_ALL_DONE_FINAL_NUMBER_OF_ITERATIONS, MPI::COMM_WORLD );
        #else
		 signaler.send_signal_i(0,TAG_SENDING_ALL_DONE_FINAL_NUMBER_OF_ITERATIONS,current_iteration_number);
       #endif
        sent = false;
//         cout << " SLAVE " << SLAVEID << " sent " << endl;
//         asd = false;
        DEBOFF("FP_iterate_use_time_limit returns false",0); 

        return false;
        
    }

    return true;
    /** inform master that this slave is done */
//     if ( !sent ) {
//         //      cout << " SLAVE " << SLAVEID <<" SENDING I AM DONE " << endl;
//         MPI_Send ( &current_iteration_number, 1, MPI_INT, 0,  TAG_SENDING_DONE_MINIMUM_NUMBER_OF_ITERATIONS, MPI::COMM_WORLD );
//         sent = true;
//     }

//     return false;


}

int min_iterlim_master_checker ( int  reset_count_final){
    /** 0 means reset,
    													*   1 means count,
    													*   2 means final,*/
// 	static int first_time = 0;
    static int retval = 0;
    static int p_Time = cout_vector_time;

//     bool temp = false;
    bool pr = false;
    if ( reset_count_final ==0 ) {
        time_temp.start();
        miniterlim_checker.resize ( WORLDSIZE,false );
        slave_iteration_numbers.resize ( WORLDSIZE,0 );
        for ( int source=1; source < WORLDSIZE; ++source ) {
            miniterlim_checker[source] = false;
            slave_iteration_numbers[source] = 0;
        }


        retval = 0;
        DEBOFF("min_iterlim_master_checker returns 0"); 

        return 0;
    }
    if ( reset_count_final == 1 ) {
        if(time_temp.stop() >= p_Time){
            pr = true;
        }
        for ( int source=1; source < WORLDSIZE; ++source ) {

            if (pr && source == 1){

//                 vector_print<bool>(miniterlim_checker, cout);
//                 vector_print<int>(slave_iteration_numbers,cout);
//                 for ( unsigned i = 1; i < slave_iteration_numbers.size(); ++i ) {
//
//                     cout << " ["<<             miniterlim_checker[i]<<
//                     "]:"<< slave_iteration_numbers[i];
//                 }
//                 cout << endl;
                p_Time +=cout_vector_time;
                pr = false;
//                 temp = true;
            }

            if ( /*(retval < (int)miniterlim_checker.size() -2 ) &&*/ ( miniterlim_checker[source] ) ) {
                continue;
            }
            int flag = -99;
            
            PROBE(source,TAG_SENDING_DONE_MINIMUM_NUMBER_OF_ITERATIONS)
            //MPI_Status status;
            //int t = MPI_Iprobe ( source,TAG_SENDING_DONE_MINIMUM_NUMBER_OF_ITERATIONS,MPI::COMM_WORLD, &flag, &status );

			if ( flag == TAG_SIGNALER_PROBE_SUCCESSFUL ) {
                int r_i ;
                RECV_I(source, TAG_SENDING_DONE_MINIMUM_NUMBER_OF_ITERATIONS)
                //MPI_Recv ( &slave_iteration_numbers[source],1, MPI_INT, source, TAG_SENDING_DONE_MINIMUM_NUMBER_OF_ITERATIONS,MPI::COMM_WORLD, &status );
                slave_iteration_numbers[source] = r_i;

                if ( slave_iteration_numbers[source] >= FP_min_iter_lim ) {
                    if ( !miniterlim_checker[source] ) {
                        retval++;
                    }
                    miniterlim_checker[source] = true;
                }
            }
           

        }
//         if (temp){
//             cout << " returning " << retval;
//             temp = false;
//         }
        DEBOFF("min_iterlim_master_checker returns ", retval );
        return retval;
    }
//  	cout <<" min_iterlim_master_checker retval: " << retval << endl;
    if ( reset_count_final == 2 ) {

        p_Time = cout_vector_time;
        pr = false;

        /** call everyone that we are done */
        cout <<  "call everyone that we are done " << endl;
        #ifdef USING_MPI
        for ( int dest=1; dest <WORLDSIZE; ++dest ) {
            
            MPI_Send ( &retval,1,MPI_INT,dest,TAG_ALL_DONE_MINIMUM_NUMBER_OF_ITERATIONS,MPI::COMM_WORLD );

        }
        #else 
        	 
		signaler.send_signal_i(DESTINATION_ALL_SLAVES, TAG_ALL_DONE_MINIMUM_NUMBER_OF_ITERATIONS,retval);
       #endif
//         usleep(1000000);
        /** now receive final iteration counts */



        int n_rec= 0;
        vector<bool> received;
        received.resize(WORLDSIZE, false);
        while (n_rec < WORLDSIZE - 1){
            for ( int dest=1; dest < WORLDSIZE; ++dest ) {
                if (received[dest]) continue;
                int flag = -99;
                
                PROBE(dest,TAG_SENDING_ALL_DONE_FINAL_NUMBER_OF_ITERATIONS)
//                MPI_Status status;
  //              int t = MPI_Iprobe(dest,TAG_SENDING_ALL_DONE_FINAL_NUMBER_OF_ITERATIONS,MPI::COMM_WORLD, &flag, &status );
                
				if (flag == TAG_SIGNALER_PROBE_SUCCESSFUL){
                    int r_i; 
                    RECV_I(dest, TAG_SENDING_ALL_DONE_FINAL_NUMBER_OF_ITERATIONS)
//                    MPI_Recv ( &slave_iteration_numbers[dest],1, MPI_INT, dest, TAG_SENDING_ALL_DONE_FINAL_NUMBER_OF_ITERATIONS,MPI::COMM_WORLD, &status );
                    slave_iteration_numbers[dest] = r_i;
                    n_rec++;
                    received[dest] = true;
                }

            }
        }

        /**   */
//         for ( int dest=1; dest < MPI::COMM_WORLD.Get_size(); ++dest ) {
//             cout << " receiving from " << dest << endl;
//             MPI_Status status;
//             MPI_Recv ( &slave_iteration_numbers[dest],1, MPI_INT, dest, TAG_SENDING_ALL_DONE_FINAL_NUMBER_OF_ITERATIONS,MPI::COMM_WORLD, &status );
//             cout << " received " << slave_iteration_numbers[dest]<< " from dest " << dest <<  endl;
//         }


        /** print iteration counts */
		cout << "line "<<__LINE__ << "  ";

        for ( unsigned i = 1; i < slave_iteration_numbers.size(); ++i ) {

            cout << " ["<<             miniterlim_checker[i]<<
                 "]:"<< slave_iteration_numbers[i];
        }
        cout << endl;
        cout << "min_iterlim_master_checker returns " << retval << " at line " << __LINE__ << endl; 
        return retval;
    }
    cout << "min_iterlim_master_checker returns " << -1 << " at line " << __LINE__ << endl; 
    return -1;

}

int time_limit_master_checker(int  reset_count_final){
    static int retval = 0;
    static int p_Time = cout_vector_time;
    bool pr = false;

    if (reset_count_final == 0){
        time_temp.start();

        master_level_FP_timer.start();
        miniterlim_checker.resize ( WORLDSIZE,false );
        slave_iteration_numbers.resize ( WORLDSIZE,0 );
        for ( int source=1; source < WORLDSIZE; ++source ) {
            miniterlim_checker[source] = false;
            slave_iteration_numbers[source] = 0;
        }


        retval = 0;
        //cout << "time_limit_master_checker returns " << retval << " at line " << __LINE__ << endl; 
        return 0;
    }

    if ( reset_count_final == 1 ) {
        if(time_temp.stop() >= p_Time){
            pr = true;
        }
        for ( int source=1; source < WORLDSIZE; ++source ) {
            if (pr && source == 1){



                //                 vector_print<bool>(miniterlim_checker, cout);
                //                 vector_print<int>(slave_iteration_numbers,cout);
// 				cout << "\t\t\t\t\t\t\t NOT DONE LIST: ";
//                 for ( unsigned i = 1; i < slave_iteration_numbers.size(); ++i ) {
// 					if(!miniterlim_checker[i]) cout << i-1  << " ";
// //                     cout << " ["<<             miniterlim_checker[i]<<
// //                     "]:"<< slave_iteration_numbers[i];
//                 }
//                 cout << endl;
//                 cout << endl;
                p_Time +=cout_vector_time;
                pr = false;
                //                 temp = true;
            }
            if ( /*(retval < (int)miniterlim_checker.size() -2 ) &&*/ ( miniterlim_checker[source] ) ) {
                continue;
            }
            int flag = -99;
            
            PROBE(source,TAG_SENDING_DONE_WITH_TIME_LIMIT)
//            MPI_Status status;
//            int t = MPI_Iprobe ( source,TAG_SENDING_DONE_WITH_TIME_LIMIT,MPI::COMM_WORLD, &flag, &status );

			if ( flag == TAG_SIGNALER_PROBE_SUCCESSFUL ) {
                int r_i;
                RECV_I(source,TAG_SENDING_DONE_WITH_TIME_LIMIT)
                
                //MPI_Recv ( &slave_iteration_numbers[source],1, MPI_INT, source, TAG_SENDING_DONE_WITH_TIME_LIMIT,MPI::COMM_WORLD, &status );
                slave_iteration_numbers[source] = r_i;
                
//                 if ( slave_iteration_numbers[source] >= FP_min_iter_lim ) {
                    if ( !miniterlim_checker[source] ) {
                        retval++;
                    }
                    miniterlim_checker[source] = true;
//                 }
            }


        }
        //         if (temp){
//                     cout << " returning " << retval;
        //             temp = false;
        //         }
        if (retval == N_SLAVES) {

//             for ( int dest=1; dest < MPI::COMM_WORLD.Get_size(); ++dest ) {
//
//                 cout << dest<< ":" << miniterlim_checker[dest ] << " ";
//
//                 }
//                 cout << endl;
        }
        //cout << "time_limit_master_checker returns " << retval << " at line " << __LINE__ << endl;
        return retval;
    }
//     cout << "line "<<__LINE__ << "  " << endl;

    if ( reset_count_final == 2 ) {
//         cout << "line "<<__LINE__ << "  " << endl;

        p_Time = cout_vector_time;
        pr = false;
        /** call everyone that we are done */
        #ifdef USING_MPI
        for ( int dest=1; dest < WORLDSIZE; ++dest ) {
            MPI_Send ( &retval,1,MPI_INT,dest,TAG_ALL_DONE_WITH_TIME_LIMIT,MPI::COMM_WORLD );
        }
        #else 
        	 signaler.send_signal_i(DESTINATION_ALL_SLAVES,TAG_ALL_DONE_WITH_TIME_LIMIT,retval);
        #endif
//         cout << "line "<<__LINE__ << "  " << endl;
//         usleep(50000);

        /** now receive final iteration counts */
        int n_rec= 0;
        vector<bool> received;
        received.resize(WORLDSIZE, false);
        while (n_rec < WORLDSIZE - 1){
            for ( int dest=1; dest < WORLDSIZE; ++dest ) {
                if (received[dest]) continue;

                int flag = -99;
                
                
                PROBE(dest,TAG_SENDING_ALL_DONE_FINAL_NUMBER_OF_ITERATIONS)
//                MPI_Status status;
  //              int t = MPI_Iprobe(dest,TAG_SENDING_ALL_DONE_FINAL_NUMBER_OF_ITERATIONS,MPI::COMM_WORLD, &flag, &status );
                
				if (flag ==TAG_SIGNALER_PROBE_SUCCESSFUL){
                    int r_i;
//                     cout << "receiving from  dest: " << dest << endl;
                    RECV_I(dest,TAG_SENDING_ALL_DONE_FINAL_NUMBER_OF_ITERATIONS)
                    slave_iteration_numbers[dest] = r_i;
//                    MPI_Recv ( &slave_iteration_numbers[dest],1, MPI_INT, dest, TAG_SENDING_ALL_DONE_FINAL_NUMBER_OF_ITERATIONS,MPI::COMM_WORLD, &status );

                    n_rec++;
                    received[dest] = true;
//                     cout << " dest: " << dest << " returns " << slave_iteration_numbers[dest] << " n_rec: " << n_rec <<endl;

                    
                }
                else {

                }

            }
        }

//         for ( int dest=1; dest < MPI::COMM_WORLD.Get_size(); ++dest ) {
//             cout <<" receiving from " << dest << endl;
//             MPI_Status status;
//             MPI_Recv ( &slave_iteration_numbers[dest],1, MPI_INT, dest, TAG_SENDING_ALL_DONE_FINAL_NUMBER_OF_ITERATIONS,MPI::COMM_WORLD, &status );
//
//             cout << " dest: " << dest << " returns " << slave_iteration_numbers[dest] << endl;
//         }
        /** print iteration counts */
//         cout << "line "<<__LINE__ << "  " << endl;
//
//         for ( unsigned i = 1; i < slave_iteration_numbers.size(); ++i ) {
//
//             cout << " ["<<             miniterlim_checker[i]<<
//             "]:"<< slave_iteration_numbers[i];
//         }
//
//         cout << endl;
        cout << "time_limit_master_checker returns " << retval << " at line " << __LINE__ << endl;
        return retval;
    }
    cout << "time_limit_master_checker returns " << -1 << " at line " << __LINE__ << endl;
    return -1;


}


Generate_cut_type cut_type_generator(){
    Generate_cut_type t;
    t.generate_all();
    t.unset_generate_type ( CUT_TYPE_SIMPLE_ROUNDING );
    t.unset_generate_type ( CUT_TYPE_GOMORY );
    t.unset_generate_type ( CUT_TYPE_ALL_DIFFERENT );

    //     t.unset_generate_type( CUT_TYPE_KNAPSACK);
    //     t.unset_generate_type( CUT_TYPE_REDSPLIT);
    t.unset_generate_type ( CUT_TYPE_L_AND_P );
    t.unset_generate_type ( CUT_TYPE_LIFT_AND_PROJECT );
    return t;
}

int compare_two_cuts ( mySerializableRowCut &f, mySerializableRowCut &s ){
    /** assumes efficiencies are calculated */
    int retval = 0;

    if ( f.get_efficiency ( CUT_EFFICIENCY_RELATIVE_VIOLATION ) > s.get_efficiency ( CUT_EFFICIENCY_RELATIVE_VIOLATION ) ) {
        retval++;
    } else {
        retval--;
    }

    if ( f.get_efficiency ( CUT_EFFICIENCY_ADJUSTED_DISTANCE ) > s.get_efficiency ( CUT_EFFICIENCY_ADJUSTED_DISTANCE ) ) {
        retval++;
    } else {
        retval--;
    }

    if ( f.get_efficiency ( CUT_EFFICIENCY_DISTANCE_VARIANT ) > s.get_efficiency ( CUT_EFFICIENCY_DISTANCE_VARIANT ) ) {
        retval++;
    } else {
        retval--;
    }

    if ( f.get_efficiency ( CUT_EFFICIENCY_OBJECTIVE_PARALLELISM_WRT_AUX ) > s.get_efficiency ( CUT_EFFICIENCY_OBJECTIVE_PARALLELISM_WRT_AUX ) ) {
        retval++;
    } else {
        retval--;
    }

    if ( f.get_efficiency ( CUT_EFFICIENCY_EXPECTED_IMPROVEMENT_WRT_AUX ) > s.get_efficiency ( CUT_EFFICIENCY_EXPECTED_IMPROVEMENT_WRT_AUX ) ) {
        retval++;
    } else {
        retval--;
    }

    if ( f.get_efficiency ( CUT_EFFICIENCY_SUPPORT ) > s.get_efficiency ( CUT_EFFICIENCY_SUPPORT ) ) {
        retval++;
    } else {
        retval--;
    }

    if ( f.get_efficiency ( CUT_EFFICIENCY_INTEGRAL_SUPPORT ) > s.get_efficiency ( CUT_EFFICIENCY_INTEGRAL_SUPPORT ) ) {
        retval++;
    } else {
        retval--;
    }


    return retval;
}

bool cut_comparison_relative_violation ( const mySerializableRowCut &lhs,const mySerializableRowCut &rhs ){
    return ( lhs.get_efficiency ( CUT_EFFICIENCY_RELATIVE_VIOLATION ) < rhs.get_efficiency ( CUT_EFFICIENCY_RELATIVE_VIOLATION ) );
}

vector<mySerializableRowCut> rowCutSelectionProcessAtSlave ( vector<mySerializableRowCut> &initial_list, int size = 1 ){
    sort ( initial_list.begin(),initial_list.end() );

    std::vector<mySerializableRowCut>::iterator it;

    it = unique ( initial_list.begin(),initial_list.end() );
    initial_list.resize ( std::distance ( initial_list.begin(),it ) );


    if ( size == 0 ) {
        return initial_list;
    }
    vector<mySerializableRowCut> final_list;
    vector<bool> to_be_added;
    to_be_added.resize ( initial_list.size(),true );
    int t;
    for ( unsigned i = 0; i < initial_list.size(); ++i ) {
        for ( unsigned j = i+1; j < initial_list.size(); ++j ) {
            t = compare_two_cuts ( initial_list[i], initial_list[j] );

            if ( t > efficiency_threshold ) {
                to_be_added[j]= false;
//                 cout << "slave " << SLAVEID << " comparison " << t << endl;

            }
            if ( t < -efficiency_threshold ) {
                to_be_added[i]= false;
//                 cout << "slave " << SLAVEID << " comparison " << t << endl;

            }
        }
        if ( to_be_added[i] ) {
            final_list.push_back ( initial_list[i] );
        }
    }
//     cout << " final_list.size(): "<< final_list.size() << " size: "<<size <<endl;
    if ( ( final_list.size() > ( unsigned ) size ) && ( size >0 ) ) {
        sort ( final_list.begin(),final_list.end(),cut_comparison_relative_violation );
        final_list.resize ( size );
//         cout << "--<< =--" <<endl;
    }


    return final_list;
}


vector<mySerializableRowCut> rowCutSelectionProcessAtMaster ( vector<mySerializableRowCut> &initial_list, int size = 1 ){
    vector<mySerializableRowCut> final_list;
    return initial_list;
}


int send_recv_to_from_master ( double &d_send_recv, bool send_recv, int to_from ){
    /** returns 0 on success
    send_recv true = > send; false => recv */
    //EASY_SEND_AND_RECEIVE mpi_signaler;
    if ( send_recv ) { /** sending */
      
        signaler.send_signal_as_tag_and_double ( TAG_SENDING_OBJECTIVE_VALUE,d_send_recv, 0 );

    } else {
        int flag = -1;
        signaler.receive_signal_as_tag_and_double ( 0 , flag);
        if (flag == TAG_SIGNALER_RETURNS_NOERROR) d_send_recv = signaler.get_double();
    }
    return 0;

}


bool create_log ( int id  ){
    int l[] = {LOG_SET};
    int size = sizeof ( l ) / sizeof ( l[0] );
    for ( int i =0; i< size ; ++i ) {
        if ( l[i] == id ) {
            return true;
        }
    }
    return false;
}

bool create_cout ( int id ){
    int l[] = {COUT_SET};
    int size = sizeof ( l ) / sizeof ( l[0] );
    for ( int i =0; i< size ; ++i ) {
        if ( l[i] == id ) {
            return true;
        }
    }
    return false;
}

int print_index ( int id  ){
    int retval = 0;
    if ( create_cout ( id ) ) {
        retval+=1;
    }
    if ( create_log ( id ) ) {
        retval+=2;
    }
    return retval;
}

int filename_generator ( stringstream &newfilename,int id = -1 ){
    newfilename.clear();
    newfilename << FILENAME_FOR_SLAVES << id << ".txt";
    return 0;
}
/*
void testing ( int type = 0,char* argv1 = NULL ){
    {
        Generate_cut_type d;
        d.generate_all();
        test_<Generate_cut_type> ( d );
        return;
    }

    if ( type >= 1 ) {
        char *f_name_lp = argv1;

        OsiClpSolverInterface *clp = new OsiClpSolverInterface;

        CglGomory gomory_cut_generator;
        CglKnapsackCover knapsack_cut_generator;
        CglRedSplit redsplit_cut_generator;
        CglSimpleRounding simplerounding_cut_generator;



        if ( strcmp ( & ( f_name_lp[strlen ( f_name_lp )-3] ), ".lp" ) == 0 ) {
            clp->readLp ( f_name_lp );
        } else {
            if ( strcmp ( & ( f_name_lp[strlen ( f_name_lp )-4] ), ".mps" ) == 0 ) {
                clp->readMps ( f_name_lp );
            } else {
                printf ( "### ERROR: unrecognized file type\n" );
                exit ( 1 );
            }
        }


        clp->initialSolve();

        const double * x = clp->getStrictColSolution();
        OsiCuts cuts;
        cout << "generating cuts" << endl;

        gomory_cut_generator.generateCuts ( *clp, cuts );
        knapsack_cut_generator.generateCuts ( *clp, cuts );

        redsplit_cut_generator.generateCuts ( *clp, cuts );
        simplerounding_cut_generator.generateCuts ( *clp, cuts );

        vector<mySerializableRowCut> my_cut_list;
        int number_of_cuts = cuts.sizeCuts() ;


        for ( int i = 0; i < number_of_cuts; ++i ) {
            mySerializableRowCut cut3 ( cuts.rowCut ( i ) );
            cut3.calculate_efficiency ( x );
            my_cut_list.push_back ( cut3 );
        }


        cout << "---------------------"<< endl;
        test_<mySerializableRowCut> ( my_cut_list[type - 1] );

    }



    else {
        vector<double> asd;
        asd.push_back ( 0 );
        for ( int i =0; i < 10; ++i ) {
            asd.push_back ( i*i/ ( double ) 100 );
        }
        My_solution sol ( asd );
        My_objective_function obj ( asd );
        test_<My_objective_function> ( obj );
    }
    return;
}
*/
int slave_receive_objective_function_CP ( cpfp_executable &runner ){

    if ( type_dist_for_CP.I_am_running_what_objective ( SLAVEID ) > 0 ) {
        My_objective_function new_objective_function;
        /** THIS IS RECV # (3) in SLAVE
         * THIS RECV MATCHES SEND # (3) in MASTER
         * ******************************/
        OBJECT_RECV ( 0,TAG_SENDING_OBJECTIVE_FUNCTION,new_objective_function )
        /** UPDATE */
        runner.set_objective_function ( new_objective_function, 3 );
    }

    return 0;
}

int slave_receive_objective_function_FP ( cpfp_executable &runner ){
    if ( type_dist_for_FP.I_am_running_what_objective ( SLAVEID ) > 0 ) {
        My_objective_function new_objective_function;
        /** THIS IS RECV # (3) in SLAVE
         * THIS RECV MATCHES SEND # (3) in MASTER
         * ******************************/
        OBJECT_RECV ( 0,TAG_SENDING_OBJECTIVE_FUNCTION,new_objective_function )
        /** UPDATE */
        runner.set_objective_function ( new_objective_function, 1);
    }
    return 0;
}

int slave_receive_rowcuts ( cpfp_executable &runner, vector<mySerializableRowCut> &allreceivedcuts ){
    //EASY_SEND_AND_RECEIVE mpi_signaler;
    /** THIS IS SIGNAL # (signal_2) from master to slave
     * RECEIVING
     * ******************************/
    int flag = -1;
    
    int n_rowcuts = signaler.receive_signal_i ( 0, TAG_SENDING_ROWCUT, flag );
    if (flag == TAG_SIGNALER_RETURNS_ERROR) n_rowcuts = 0;
    
    __SRS ( ++_x_ )


    if ( signaler.get_tag() == TAG_SENDING_ROWCUT ) {
        /** receive rowcuts if available*/
        if ( signaler.get_number() >0 ) {
            vector<mySerializableRowCut> cuts_to_be_added;

            cuts_to_be_added.resize ( signaler.get_number() );

            for ( int i =0; i < signaler.get_number(); ++i ) {
                /** THIS IS RECV # (4) in SLAVE
                 * THIS RECV MATCHES SEND # (4) in MASTER
                 * ******************************/
                OBJECT_RECV ( 0,TAG_SENDING_ROWCUT,cuts_to_be_added[i] )
            }
            /** UPDATE */

            allreceivedcuts.insert ( allreceivedcuts.end(), cuts_to_be_added.begin(), cuts_to_be_added.end() );
            time_to_apply_cuts_slave.restart();
            runner.apply_cuts ( cuts_to_be_added );
            time_to_apply_cuts_slave.pause();
        }
    }
    return 0;

}

int slave_receive_cut_type ( Generate_cut_type &cut_t ){
    //EASY_SEND_AND_RECEIVE mpi_signaler;
    static int dzz= 0;
    /** THIS IS SIGNAL # (signal_3) from master to slave
     * RECEIVING
     * ******************************/
    int flag = -1;
    int ttt = signaler.receive_signal_i ( 0, TAG_SENDING_GENERATE_CUT_TYPE, flag );
    
    if (flag == TAG_SIGNALER_RETURNS_ERROR) {
        cout << " error in receiving signal " << endl;
        return 0;
    
    }
    __SRS ( ++_x_ )

    if ( signaler.get_tag() == TAG_SENDING_GENERATE_CUT_TYPE ) {
        if ( signaler.get_number() >0 ) {
            /** receive cut type if necessary*/

            /** THIS IS RECV # (5) in SLAVE
             * THIS RECV MATCHES SEND # (5) in MASTER
             * ******************************/
            OBJECT_RECV ( 0,TAG_SENDING_GENERATE_CUT_TYPE,cut_t )
        }
    }
    IFSLAVE_y ( 1 ) {

        if ( dzz%10 == 0 ) {
            cut_t.set_generate_type ( CUT_TYPE_SIMPLE_ROUNDING );
        }
        dzz++;
    }
    return 0;
}

int slave_receive_lb ( cpfp_executable &runner ){
    /** get LB for objective cutoff constarint */
	static int lb_counter = 0;
	
    static double LB_for_original_objective = -1e99;
    /**MPI_Recv(); */
    //EASY_SEND_AND_RECEIVE my_mpi_signaler;

    /** THIS IS RECV # (1)  in SLAVE
     * THIS RECV MATCHES SEND # (1) in MASTER
     * ******************************/
    int flag = -1;
	double r_d = -1e99;
	bool update = false;
	if (lb_counter == 0){
		//DEBON ("YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYRDDDDDDDDDDDDDDDDDDDDDDD")
		r_d = signaler.receive_signal_d(0,TAG_SENDING_LOWERBOUND,flag, ALT1_RECV_BLOCKING);
		DEBON ("FIRST LB RECEIVED", r_d)
		update = true;
	}
 	else{
		PROBE(0,TAG_SENDING_LOWERBOUND)
		if(flag == TAG_SIGNALER_PROBE_SUCCESSFUL){
			RECV_D(0,TAG_SENDING_LOWERBOUND)
			if (double_equality(r_d,TAG_SIGNALER_RETURNS_ERROR)){
				cout << " couldn't receive LB " << endl; 
			}
			else{
				update = true;
			}
		}
	}
	
	
	if (update){
		//cout << " LB received as " << r_d << endl;
		lb_counter ++;
		LB_for_original_objective = r_d;			
		double current_lb = runner.get_objective_cutoff_constraint_lb();
		if (double_significantly_more(LB_for_original_objective, current_lb)){
			__SRD ( ++_x_ )
			
			runner.update_objective_cutoff_constraint_lb ( LB_for_original_objective );
			//cout << " LB updated as  " << LB_for_original_objective <<" for the lb_counter = "<< lb_counter<< "th time. IT WAS " <<current_lb << endl;
		}
	}
	else{
			if (lb_counter == 0 && !update){
				cout << " FLAG " << flag << " LB needs to be received at least once  " <<endl;
				}
			}
			
			
			
		
		
		
// 	}
    //double LB = signaler.receive_signal_d ( 0, TAG_SENDING_LOWERBOUND, flag );
    
    //if (flag == TAG_SIGNALER_RETURNS_ERROR) return  TAG_SIGNALER_RETURNS_ERROR;
    
    //DEBOFF(LB_for_original_objective)
    // LB_for_original_objective = LB;

    //runner.update_objective_cutoff_constraint_lb ( LB_for_original_objective );
    return 0;

}

int slave_receive_best_objective_so_far_FP ( cpfp_executable &runner ){
    stringstream f_out;
    if ( create_log(SLAVEID) ) {
        f_out << "slave_FP_id_" << SLAVEID;
    }
    ofstream s_out ( f_out.str().c_str(),ios_base::app );

    //EASY_SEND_AND_RECEIVE mpi_signaler;
    int flag = -1;
	
	PROBE(0,TAG_SENDING_OBJECTIVE_VALUE_IN_FP);
	
	if (flag == TAG_SIGNALER_PROBE_SUCCESSFUL){
		double r_d = 1e99;
		RECV_D(0, TAG_SENDING_OBJECTIVE_VALUE_IN_FP);
		if (double_equality(r_d,TAG_SIGNALER_RETURNS_ERROR)){
			cout << " couldn't receive UB " << endl; 
		}
		else{
			
			double best_incumbent_objective = r_d;
			
			if (best_incumbent_objective < 1e30){
// 				cout << "slave " << SLAVEID << " received best_incumbent_objective as " << best_incumbent_objective << " it was " << latest_updated_fp_objective_in_slave << endl; 
// 				cout << "slave " << SLAVEID << "runner.get_objective_cutoff_constraint_ub() BEFORE: "<< runner.get_objective_cutoff_constraint_ub() << endl;
				UPDATE_RUNNER_FP(best_incumbent_objective,runner);
// 				cout << "slave " << SLAVEID << "runner.get_objective_cutoff_constraint_ub() AFTER : "<< runner.get_objective_cutoff_constraint_ub() << endl;
			}
		}
			
	}
    
    
    return 0;
}

int slave_receive_best_objective_so_far_CP ( cpfp_executable &runner ){
    //EASY_SEND_AND_RECEIVE mpi_signaler;
	DEBON ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAS")
	//// NEEDS TO BE REDONE WRT slave_receive_best_objective_so_far_FP
    int flag = -1;
    double OBJS = signaler.receive_signal_d(0, TAG_SENDING_OBJECTIVE_VALUE, flag);
    
    if (flag == TAG_SIGNALER_RETURNS_ERROR) return  TAG_SIGNALER_RETURNS_ERROR;
    
    if ( signaler.get_tag() == TAG_SENDING_OBJECTIVE_VALUE ) {
        
		double best_incumbent_objective = signaler.get_double();
		double llb = runner. get_objective_cutoff_constraint_lb();
		UPDATE_RUNNER_CP(best_incumbent_objective,llb,runner);
		/** TODO: check if it is done correctly*/
//         runner.update_objective_cutoff_constraint_ub ( best_incumbent_objective - objective_improvement_coefficient_CP );
    }

    return 0;
}

int slave_CP ( cpfp_executable &runner, My_solution &arg_current_lp_optimum_at_slave, vector<mySerializableRowCut> &allreceivedcuts, int &previous_number_of_cuts ){

    static int iteration_counter_slave = 0;
    iteration_counter_slave++;
    //EASY_SEND_AND_RECEIVE mpi_signaler;

    Generate_cut_type cut_t;
    if ( !slave_generates_own_functions_for_CP ) {
        slave_receive_objective_function_CP ( runner );
        /** receive new objective function if necessary */

    } else {
        slave_create_objective_CP ( runner, arg_current_lp_optimum_at_slave );
    }

    /** receive rowcuts if available */
    slave_receive_rowcuts ( runner,allreceivedcuts );

    /** receive cut_type if available */
    slave_receive_cut_type ( cut_t );


    /** receive best_incumbent_objective so far if available */
    slave_receive_best_objective_so_far_CP ( runner );

    /** OPTIMIZE */
    time_to_optimize_in_slave.restart();
    int res = runner.optimize();
    time_to_optimize_in_slave.pause();
    stringstream s;


    switch ( res ) {
//         case -1:
//             CPFP_ERROR ( "lp neither optimal nor infeasible ", cout )
//             break;
    case 0: /**infeasible: CAN THIS BE THE CASE ?*/
        CPFP_ERROR ( "lp infeasible ", cout )

        cout << " slave " << SLAVEID << " lb: " << runner.get_objective_cutoff_constraint_lb() <<" ub: " << runner.get_objective_cutoff_constraint_ub() << endl;
        s <<  "slave_" << SLAVEID;
        s << "_iteration_" << iteration_counter_slave<< "_infeas";

//         runner.puke_lp_file ( s.str().c_str() );

        break;
    case 1: /** feasible */
        break;
    case -6 :
        CPFP_ERROR ( "lp neither optimal nor infeasible ", cout )
        cout << "res is " << res << endl;


        BWAIT
//         runner.puke_lp_file ( s.str().c_str() );
        break;
    default :

        CPFP_ERROR ( "lp neither optimal nor infeasible ", cout )
        cout << "res is " << res << endl;
        break;
    }

    My_solution sol;
    vector<mySerializableRowCut> slave_rowcutlist;
    double integer_infeasibility;


    int int_feasibility = runner.get_solution ( sol,integer_infeasibility );
    if ( ( res == 1 ) && ( int_feasibility == 0 ) ) {
        /** FOUND a feasible solution
         * no need to generate cuts
         */

    } else {
        if ( res == 1 ) {
            /** NOT a feasible solution
            *  generate row cut list
            */
            time_to_generate_cuts.restart();
            slave_rowcutlist = runner.generate_rowcutlist ( cut_t, allreceivedcuts, previous_number_of_cuts );
            time_to_generate_cuts.pause();

            for ( unsigned i =0; i<slave_rowcutlist.size(); ++i ) {
                slave_rowcutlist[i].set_efficiency ( CUT_EFFICIENCY_BGN, MACRO_SLAVE_ID );
            }
        }
    }

    /** IT IS TIME TO INFORM MASTER ABOUT SOLUTION OR ROW CUT LIST
     */
    if ( ( res == 1 ) && ( int_feasibility == 0 ) ) {

        /** since it is feasible sent the solution */

        /** THIS IS SIGNAL # (signal_5) from slave to master
         * SENDING
         * NOTE : This is only the IF part
         * NOTE : so ther is another signal_5 in ELSE
         * ******************************/

        signaler.send_signal_i ( 0, TAG_SENDING_SOLUTION,1 );



        /** THIS IS SEND # (6) in SLAVE
         * THIS SEND MATCHES RECV # (6) in MASTER
         * ******************************/
        OBJECT_SEND ( 0,TAG_SENDING_SOLUTION,sol )

        /** also tell the master that no rowcut will be sent
         */

        /** THIS IS SIGNAL # (signal_6) from slave to master
         * SENDING
         * NOTE : This is only the IF part
         * NOTE : so ther is another signal_6 in ELSE
         * ******************************/
        signaler.send_signal_i (0, TAG_SENDING_ROWCUT,1 );


    } else {
        if ( res == 1 ) {
            /** since it is not feasible send the generated cuts
            * as there are no feasible solution in the current setting
            * first tell the master that there is no solution to be sent
            * */

            /** THIS IS SIGNAL # (signal_5) from slave to master
            * SENDING
            * NOTE : This is ELSE of the previous if
            * NOTE : so ther is another signal_5
            * ******************************/

            

            /** Now it is time to send the set of cuts
            */
            vector<mySerializableRowCut> cuts_to_be_sent;

            if ( type_dist_for_CP.I_am_running_what_objective ( SLAVEID ) == OBJECTIVE_TYPE_ORIGINAL ) {
                cuts_to_be_sent = rowCutSelectionProcessAtSlave ( slave_rowcutlist, 0 );
            } else {
                cuts_to_be_sent = rowCutSelectionProcessAtSlave ( slave_rowcutlist, number_of_cuts_per_nonoriginal_slave );


            }
            int cutsize = ( int ) cuts_to_be_sent.size();
            if (cutsize >0){
                
                signaler.send_signal_i(0, TAG_SENDING_SOLUTION, cutsize);
                
                /** THIS IS SIGNAL # (signal_6) from slave to master
                * SENDING
                * NOTE : This is ELSE of the previous if
                * NOTE : so ther is another signal_6
                * ******************************/

                for ( int i = 0; i < cutsize; ++i ) {
                    /** THIS IS SEND # (7) in SLAVE
                    * THIS SEND MATCHES RECV # (7) in MASTER
                    * ******************************/
                    OBJECT_SEND ( 0,TAG_SENDING_ROWCUT, cuts_to_be_sent[i] )
                }
                
            }

        } else {
            /** THIS IS SIGNAL # (signal_5) from slave to master
             * SENDING
             * NOTE : This is ELSE ELSE of the previous if
             * NOTE : so ther is another signal_5
             * *******************************/
            signaler.send_signal_i(0, TAG_SENDING_SOLUTION,0 );
            /** THIS IS SIGNAL # (signal_6) from slave to master
             * SENDING
             * NOTE : This is ELSE ELSE of the previous if
             * NOTE : so ther is another signal_6
             * ******************************/
            signaler.send_signal_i(0, TAG_SENDING_ROWCUT,0);


        }
    }


    return 0;
}



int slave_generate_FP_starting_solution(cpfp_executable &runner, My_solution &arg_current_lp_optimum_at_slave, My_solution &fill_in_this_starting_solution, int &int_frac, double &double_frac){
	if (type_dist_for_FP.I_am_running_what_objective(SLAVEID) == 0){
		My_objective_function new_obj = runner.get_original_objective_function();
		runner.set_objective_function ( new_obj , 1);
// // 		runner.optimize();

// 		cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
// 		fill_in_this_starting_solution = arg_current_lp_optimum_at_slave;
//
// 		int_frac = runner.get_fractionality_of_solution(fill_in_this_starting_solution,double_frac);
// 		return 1;
	}
	else{
		if ( !slave_generates_own_functions_for_FP ) {
			slave_receive_objective_function_FP ( runner );
			/** receive new objective function if necessary */
// 			cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
		} else {
			/** create an objective function*/
// 			cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
			slave_create_objective_FP ( runner, arg_current_lp_optimum_at_slave );
		}
	}
// 	cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
    int res = runner.optimize();

    switch(res){
        case 1:
            int_frac = runner.get_solution(fill_in_this_starting_solution,double_frac);
// 			cout << "line: " << __LINE__ << endl;
// 			fill_in_this_starting_solution.print();
// 			cout << "line: " << __LINE__ << endl;

            break;
        case -5:
//             break;
        default:
            stringstream s;
            s<< "problematic_lp_file_"<<p_name<<"_slave_" << SLAVEID << "_res_" << res;
            runner.puke_lp_file(s.str().c_str());
            cout << " res is :" << res << endl;
        break;

    }

    return res;

}


int slave_FP_v4(cpfp_executable &runner, My_solution &arg_current_lp_optimum_at_slave){
	
// 	IFSLAVE_y(1) cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
	read_FP_parameters_from_file(runner);
 	//IFSLAVE cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
 	//IFSLAVE runner.print_FP_parameters(std::cout);
	
	time_for_FP.start();

	My_solution starting_solution;
	My_solution rounded_solution;
	My_solution previous_rounded_solution;
	My_solution best_incumbent_solution_at_slave_FP;
	My_solution best_starting_solution_so_far;
	My_solution best_rounded_solution_so_far;
	
	int stage1_iterlim;
	int stage2_iterlim;
	int stage1_resetlim;
	int stage2_resetlim;
	int tot_found=0;
	int latest_find_iteration = -1;
	if(0)latest_find_iteration++;
	double temp = latest_updated_fp_objective_in_slave;
	slave_receive_best_objective_so_far_FP ( runner );

	//DEBON();
	if (runner.get_objective_cutoff_constraint_lb() > runner.get_objective_cutoff_constraint_ub()){
		cout <<  " slave " << SLAVEID << " .. LB > UB" <<endl;
	}
	//DEBON("slave_generates_own_functions_for_FP", slave_generates_own_functions_for_FP);
	
	/** INITIAL OBJECTIVE SETTING IS DONE HERE*/
	if ( !slave_generates_own_functions_for_FP ) {
        DEBON();
		slave_receive_objective_function_FP ( runner );
        DEBON();
		/** receive new objective function if necessary */
		
	} else {
		/** create an objective function*/
        //arg_current_lp_optimum_at_slave.print();
        //IFSLAVE_y(2) {DEBWAIT()}
		time_to_create_starting_points_for_FP_slave.restart();
		slave_create_objective_FP ( runner, arg_current_lp_optimum_at_slave );
        
		time_to_create_starting_points_for_FP_slave.pause();
		//IFSLAVE_y(2) {DEBWAIT()}
	}
	
	
	//DEBON();
	
	
	runner.get_FP_INT_PARAM(ENM_FP_INT_PARAM_MAXITER1, stage1_iterlim);
	runner.get_FP_INT_PARAM(ENM_FP_INT_PARAM_MAXITER2, stage2_iterlim);
	runner.get_FP_INT_PARAM(ENM_FP_INT_PARAM_MAXITERW1, stage1_resetlim);
	runner.get_FP_INT_PARAM(ENM_FP_INT_PARAM_MAXITERW2, stage2_resetlim);
	runner.get_FP_INT_PARAM(ENM_FP_INT_PARAM_ROUNDING, initial_rounding_type);
	
 	runner.get_FP_DOUBLE_PARAM(ENM_FP_DOUBLE_PARAM_ALPHA_INIT, initial_alpha);
	runner.get_FP_DOUBLE_PARAM(ENM_FP_DOUBLE_PARAM_ALPHA_QUOD, alpha_reduction);
	runner.get_FP_DOUBLE_PARAM(ENM_FP_DOUBLE_PARAM_ALPHA_DIST, alpha_dist);

	int rounding_type = initial_rounding_type;
	
	int ncols = runner.get_ncols();
	
	int nBinInt = runner.get_number_of_binaries();
	int nGenInt = runner.get_number_of_general_integers();
	int nCont   = ncols - nBinInt - nGenInt;
	runner.init_FP_stage2();
	
	int current_stage = 1;
	int missedDecr = 0;
	int maxMissedDecr = stage1_resetlim;
	int minFracIt = 0; if(0) minFracIt++;
	double minFracDbl = std::numeric_limits< double >::max();
	
	double best_incumbent_objective_at_slave_FP = std::numeric_limits< double >::max();
	
	starting_solution.set_original_size(ncols);
	rounded_solution.set_original_size(ncols);
	previous_rounded_solution.set_original_size(ncols);
	best_incumbent_solution_at_slave_FP.set_original_size(ncols);
	
	best_starting_solution_so_far.set_original_size(ncols);
	
	external_column_types = runner.get_column_types();
	
	int big_iteration_counter = 0;
	int external_iteration_counter =0;
	int FP_total_iteration_counter = 0;
	bool printing_breaks = false;
	
	bool fp_continue_flag = true;
	
	int stage_1_restarts = 0;
	int stage_2_restarts = 0;
	int total_s1_restarts = 0;
	int total_s2_restarts = 0;
	bool original_continues = false;
	
	
	bool will_send = false;
	bool found = false;
	
	alpha = initial_alpha;
	
// 	int aggressive_restart_count =0;
// 	cout << "initial_alpha "<< initial_alpha <<endl; 
	while (fp_continue_flag){
		//IFSLAVE_y(2) usleep(1000);
		
		external_iteration_counter++;
		big_iteration_counter++;
		double curFrac_d = std::numeric_limits< double >::max();
		int curFrac_i = ncols;
		if (big_iteration_counter > 1 && type_dist_for_FP.I_am_running_what_objective(SLAVEID) == OBJECTIVE_TYPE_ORIGINAL) {
			original_continues = true;
		}
		if (external_iteration_counter == 1 && !original_continues){
			/** create initial starting solution*/
			int opt_res = slave_generate_FP_starting_solution(runner, arg_current_lp_optimum_at_slave, starting_solution, curFrac_i,curFrac_d);
			if (opt_res<=0) {
				cout << "problem in finding an initial solution line: "<< __LINE__ << " opt_res: "<< opt_res << " external_iteration_counter: " << external_iteration_counter << " big_iteration_counter: "<< big_iteration_counter << endl;
				/** TO DO*/ 
// 				break;
			}
			stage_1_restarts =0;
			stage_2_restarts =0;
			best_starting_solution_so_far.clear();
			best_rounded_solution_so_far.clear();
			minFracDbl = curFrac_d;
			minFracIt = FP_total_iteration_counter;
			if(curFrac_i ==0){ /** int sol*/
				double ext_test_dbl;
				int ext_test_int = 			external_frac_tester( starting_solution, ext_test_dbl, __LINE__);
				if (ext_test_int > 0) {
					cout << "line " << __LINE__<< " external_frac_tester returns: " << ext_test_int << " " << ext_test_dbl << endl;
				}
				
				double new_obj;
				if (runner.calculate_original_objective_value_of_solution(new_obj, starting_solution) != 0) {
					cout << "problem in line: " << __LINE__ << endl;
				}
				cout << "slave_FP v3 finds a solution in line "<<__LINE__ << " with obj: " << starting_solution.get_original_objective() << endl;
				if (current_stage ==2){
					runner.unset_bounds_FP_stage_2();
				}
				
				if (nCont >0){
					My_solution polished_sol;
					double polished_obj;
					int pol_res  = runner.polish_integer_solution(starting_solution, polished_sol, polished_obj);
					
					if (pol_res ==1){
						starting_solution = polished_sol;
					}
				}
				best_starting_solution_so_far = starting_solution;
				best_rounded_solution_so_far = rounded_solution;
// 				cout << "slave " << SLAVEID<< " FP  finds a solution with obj: " << starting_solution.get_original_objective() << endl;
				if(best_incumbent_objective_at_slave_FP> starting_solution.get_original_objective()){
					best_incumbent_objective_at_slave_FP= starting_solution.get_original_objective();
					best_incumbent_solution_at_slave_FP = starting_solution;
					found = true;
				}
			}
			else{
				previous_rounded_solution =  rounded_solution;
				rounding_type = rounding_decider(stage_1_restarts,current_stage);
				runner.round_sol(starting_solution, rounded_solution, rounding_type,current_stage);
			}
		}
		else{
			alpha = initial_alpha;
			history_clear_version1();
			best_rounded_solution_so_far.clear();
// 			if (tot_found==0 || latest_find_iteration < (FP_total_iteration_counter - 1000)) {
// 				aggressive_restart_count++;
// // 				cout << "slave "<< SLAVEID  << "  aggressive_flip for the " << aggressive_restart_count << "th time n_changed: " <<  runner.aggressive_flip(rounded_solution, aggressive_flip_probability)<< endl;
// 			}
// 			rounded_solution.clear();
// 			previous_rounded_solution.clear();
//  			rounded_solution = best_rounded_solution_so_far;
		}
		
		/** starting solution and rounded solution are generated*/
		/** STAGE 1*/
		if(!found && (nBinInt>0) && stage1_iterlim >0 && fp_continue_flag  &&(!slave_optimality)){
			missedDecr =0;
			minFracDbl = std::numeric_limits< double >::max();
			minFracIt = FP_total_iteration_counter;
			
			int local_iteration = 0;
			int local_iterlim = stage1_iterlim;
			maxMissedDecr = stage1_resetlim;
			current_stage = 1;
			stage_1_restarts =0;
			
			curFrac_d = std::numeric_limits< double >::max();
			curFrac_i = ncols;
			while ((local_iteration < local_iterlim) && (fp_continue_flag) && (!slave_optimality)){
				FP_total_iteration_counter++;
				local_iteration++;
				USE_COMM_V3
				
				REDUCE_ALPHA
				//DEBON()
				CREATE_OBJ_AND_OPTIMIZE(runner)
				//DEBON()
				
				curFrac_i = runner.get_solution_stage(starting_solution, curFrac_d, current_stage) ;
				//DEBON()
				best_rounded_solution_so_far = rounded_solution;
				rounding_type = rounding_decider(stage_1_restarts, current_stage);
				previous_rounded_solution =  rounded_solution;
				int n_changed = runner.round_sol(starting_solution, rounded_solution, rounding_type ,current_stage);
				
				UPDATE_FP_STATUS_r
				
				BREAK_STAGE_ONE
				
				if (n_changed ==0){
					n_changed += runner.getNextIntegerPoint(starting_solution, rounded_solution, current_stage);
				}
				HISTORY_CHECK_STAGE(stage_1_restarts)
				fp_continue_flag = FP_iterate_decider_slave( FP_total_iteration_counter, time_for_FP.stop());				
			}
		
		
		int res_ext = runner.optimize_extra_mip(starting_solution, current_stage);
		if (res_ext == 1){
			My_solution extra_sol = runner.get_solution_of_extra_mip(current_stage);
			starting_solution = extra_sol;
		}
		curFrac_i = runner.get_fractionality_of_solution(starting_solution,curFrac_d);
		
		if (curFrac_i == 0) {
			double new_obj;
			if (runner.calculate_original_objective_value_of_solution(new_obj,starting_solution) != 0) {
				cout << "problem in line: " << __LINE__ << endl;
			}
			if (current_stage ==2){
				runner.unset_bounds_FP_stage_2();
			}
			
			if (nCont > 0){
				My_solution polished_sol;
				double polished_obj;
				int pol_res  = runner.polish_integer_solution(starting_solution,polished_sol, polished_obj);
				if (pol_res ==1){
					starting_solution = polished_sol;
				}
			}
// 			best_starting_solution_so_far = starting_solution;
			best_rounded_solution_so_far = rounded_solution;
			/** found a feasible solution*/
			if(best_incumbent_objective_at_slave_FP> starting_solution.get_original_objective()){
				best_incumbent_objective_at_slave_FP = starting_solution.get_original_objective();
				best_incumbent_solution_at_slave_FP = starting_solution;
				found = true;
			}
		}
		history_clear_version1();
		}
		/** STAGE 2*/
		if (!found && stage2_iterlim>0 && fp_continue_flag && (!slave_optimality)){
			missedDecr =0;
			minFracDbl = std::numeric_limits< double >::max();
			minFracIt = FP_total_iteration_counter;
			
			int local_iteration = 0;
			int local_iterlim = stage2_iterlim;
			maxMissedDecr = stage2_resetlim;
			current_stage = 2;
			stage_2_restarts =0;
			curFrac_d = std::numeric_limits< double >::max();
			curFrac_i = ncols;

			rounded_solution = best_rounded_solution_so_far;
// 			rounded_solution = best_starting_solution_so_far;
			while ((local_iteration < local_iterlim) && (fp_continue_flag) && !slave_optimality){
				
				FP_total_iteration_counter++;
				local_iteration++;
				USE_COMM_V3
				REDUCE_ALPHA
				CREATE_OBJ_AND_OPTIMIZE(runner)
				
				curFrac_i = runner.get_solution_stage(starting_solution, curFrac_d, current_stage);
				
				rounding_type = rounding_decider(stage_2_restarts,current_stage);
				previous_rounded_solution =  rounded_solution;
				int n_changed = runner.round_sol(starting_solution, rounded_solution,rounding_type, current_stage);
	
				UPDATE_FP_STATUS_r
// 				best_rounded_solution_so_far = rounded_solution;
				
				if(curFrac_i == 0) { /** we have a feasible solution */
// 					double zikkim;
// 					int ext_res = external_frac_tester(starting_solution,zikkim,__LINE__);
// 
// 					if (ext_res > 0){
// 						cout << "external_frac_tester returns " << ext_res << " " << zikkim <<endl;
// // 						starting_solution.print_selected_indices(integer_columns_in_slave,false,cout);
// 						
// 					}
					double new_obj = 1e98;
					if (runner.calculate_original_objective_value_of_solution(new_obj, starting_solution) != 0) {
						cout << "problem in line: " << __LINE__ << endl;
					}
					if (current_stage ==2){
						runner.unset_bounds_FP_stage_2();
					}
					if (nCont > 0){
						My_solution polished_sol;
						double polished_obj;
						int pol_res  = runner.polish_integer_solution(starting_solution,polished_sol, polished_obj);
						if (pol_res ==1){
							starting_solution = polished_sol;
						}
					}
					best_starting_solution_so_far = starting_solution;
					best_rounded_solution_so_far = rounded_solution;
					
					/** found a feasible solution*/
					if(best_incumbent_objective_at_slave_FP> starting_solution.get_original_objective()){
						best_incumbent_objective_at_slave_FP= starting_solution.get_original_objective();
						best_incumbent_solution_at_slave_FP = starting_solution;
						found = true;
					}
					break;
				}
				BREAK_STAGE_TWO
				
				if (n_changed ==0){
					n_changed += runner.getNextIntegerPoint(starting_solution, rounded_solution, current_stage);
				}
				
				HISTORY_CHECK_STAGE(stage_2_restarts)
				fp_continue_flag = FP_iterate_decider_slave(FP_total_iteration_counter, time_for_FP.stop());
			}
			curFrac_i =runner.get_fractionality_of_solution(starting_solution,curFrac_d);
			if (curFrac_i == 0) {
				double new_obj = 1e30;
				if (runner.calculate_original_objective_value_of_solution(new_obj,starting_solution) != 0) {
					cout << "problem in line: " << __LINE__ << endl;
				}
				if(best_incumbent_objective_at_slave_FP > starting_solution.get_original_objective()){
					best_incumbent_objective_at_slave_FP = starting_solution.get_original_objective();
					best_incumbent_solution_at_slave_FP = starting_solution;
					found = true;
				}
				
			}
			history_clear_version1();
		}
		if(found){
			tot_found++;
			
// 			cout << "slave " << SLAVEID << " tot_found: " << tot_found <<endl; 
			latest_find_iteration = FP_total_iteration_counter;
			UPDATE_RUNNER_FP(best_incumbent_objective_at_slave_FP,runner)
			OPTIMALITY_CHECK(runner)

			found = false;
			will_send = true;
		}
		
		total_s1_restarts += stage_1_restarts;
		total_s2_restarts += stage_2_restarts;
		
		USE_COMM_V3
		
		
		if(external_iteration_counter>0){
// 			static int tttt=0; 
// 			tttt++;
// 			cout << " slave "<< SLAVEID << " goes here " <<tttt <<"th time " <<endl;
			best_rounded_solution_so_far.clear();
			alpha = initial_alpha;
			history_clear_version1();
// 			previous_rounded_solution.clear
		}
		
		
		if (slave_optimality) {
			fp_continue_flag = FP_iterate_decider_slave(FP_total_iteration_counter, time_for_FP.stop());
			// 			fp_continue_flag = FP_iterate_decider_slave(FP_max_iter_lim, FP_max_time_limit);
		}
		else {
			if (fp_continue_flag) fp_continue_flag = FP_iterate_decider_slave(FP_total_iteration_counter, time_for_FP.stop());
		}
		
		
	}
	
	cout << " slave "<< SLAVEID <<  " finalizes FP_V3 at line "<< __LINE__ << " with obj: " << best_incumbent_objective_at_slave_FP; 
	
	if (slave_optimality) cout <<" optimality proven ";
	cout <<endl;
	
	FP_iterate_decider_slave(-1,-1);
	runner.undo_FP_stage2();
	return 0;
}
int slave_run ( int argc, char* argv[], string argfilename ="__NULL" ){

	DEBON();
// 	cout << " line " <<__LINE__ <<endl;

// 	struct rlimit rl;

// 	cout << getrlimit(RLIMIT_NPROC , &rl) << endl;

// 	cout << "Default value is : "<< rl.rlim_cur << endl;
// 	cout << "Default value is : "<< rl.rlim_max << endl;
	pthread_create(&tth,&att,threadrun_slave_global_time_checker, (void*)NULL);
	
// 	if ((type_dist_for_CP.I_am_running_what_objective(SLAVEID)!= OBJECTIVE_TYPE_ORIGINAL) || (type_dist_for_FP.I_am_running_what_objective(SLAVEID)!= OBJECTIVE_TYPE_ORIGINAL)){
// 	pthread_create(&tth2,&att2,threadrun_receive_analytic_center_in_slave, (void*)NULL);
// 	}
		
// 	cout << pthread_create(&tth,&att,threadrun_slave_global_time_checker, (void*)NULL) << endl;
// 	cout << " line " <<__LINE__ <<endl;

    cpfp_executable my_runner ( argfilename );
    slave_seed = my_seed_generator.lth_seed_lcm ( MACRO_SLAVE_ID);
	
    my_runner.set_seed ( slave_seed );
	//DEBON();
	if ( my_runner.initialize(solver_type,ipopt_solver) == 1 ) {
		objective_improvement_coefficient_FP = 1;
		objective_improvement_coefficient_CP = 1;
        
    }
    //DEBON();
    vector<int> column_types = my_runner.get_column_types();
    for (unsigned i =0; i < column_types.size();++i){
        if (column_types[i]>0) integer_columns_in_slave.push_back(i);
    }


    double sb = my_runner.calculate_pseudo_bound();
    my_runner.set_pseudo_bounds(sb);
    stringstream sname;
    filename_generator ( sname, MACRO_SLAVE_ID );
    bool printing_to_file = create_log(SLAVEID);
// 	bool printing_to_scrn = create_cout();
    // int print_place = print_index();
    if ( printing_to_file ) {
        ofstream out ( sname.str().c_str(),ios_base::app );
    }
    //EASY_SEND_AND_RECEIVE my_mpi_signaler;

    //DEBON()
    slave_receive_lb ( my_runner );
    //DEBON()
    int signal_coming_from_master = 1;
#ifdef USING_MPI
    req_stat_vectors r;
#endif
    Generate_cut_type cut_type;
    int iteration_counter_slave = 0;
    int prev_n_cuts = 0;
    vector<mySerializableRowCut> allreceivedcuts;

    My_solution current_lp_optimum_at_slave;
    stringstream f_out;
    if ( create_log(SLAVEID) ) {
        f_out << "slave_FP_id_" << SLAVEID;
    }
    ofstream s_sout ( f_out.str().c_str(),ios_base::app );

    int iteration_count_on_slaves_main_code =0;
    if ( create_log(SLAVEID) ) s_sout << "slaves log count " << iteration_count_on_slaves_main_code  << " line " << __LINE__<< endl;

// 	cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
    //DEBON()
    while ( signal_coming_from_master >0 ) {
        iteration_count_on_slaves_main_code++;
        if ( create_log(SLAVEID) ) s_sout << "slaves log count " << iteration_count_on_slaves_main_code  << " line " << __LINE__<< endl;
//        DEBON("signal_coming_from_master: ", signal_coming_from_master)
        iteration_counter_slave++;

        slave_receive_lb ( my_runner );
        //DEBON("iteration_counter_slave ", iteration_counter_slave, "signal_coming_from_master: ", signal_coming_from_master)
        if ( create_log(SLAVEID) ) s_sout << "slaves log count " << iteration_count_on_slaves_main_code  << " line " << __LINE__<< endl;
 		//cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
        
        OBJECT_RECV(0,TAG_SENDING_SOLUTION,current_lp_optimum_at_slave)
        
        //current_lp_optimum_at_slave._recv_mpi ( 0,TAG_SENDING_SOLUTION,r.v_req, r.v_stat );
        /** THIS IS SIGNAL # (signal_1) from master to slave
         * RECEIVING
         * ******************************/
        if ( create_log(SLAVEID) ) s_sout << "slaves log count " << iteration_count_on_slaves_main_code  << " line " << __LINE__<< endl;
        int flag = -1;
		
		int  last_command = COMMAND_LIST_DO_OBJ_FP_OR;
        int new_command = -1;
        PROBE(0,TAG_SIGNAL_COMMAND)
		
		if ( flag == TAG_SIGNALER_PROBE_SUCCESSFUL ) {
			int r_i; 
			RECV_I(0,TAG_SIGNAL_COMMAND);
			new_command = r_i;
			last_command=  new_command;
			signal_coming_from_master = last_command;
			//DEBON("signal_coming_from_master: ", signal_coming_from_master, " flag: ", flag)
			//DEBON(r_i)
			if (r_i == TAG_SIGNALER_RETURNS_ERROR) {
				DEBON("signal_coming_from_master: ", signal_coming_from_master, " flag: ", flag)
			}
			
			
		}
        
        
        if ( create_log(SLAVEID) ) s_sout << "slaves log count " << iteration_count_on_slaves_main_code  << " line " << __LINE__<< endl;
        
//         if ( signaler.get_number() <0 ) {
//             cout << "slaves log count " << iteration_count_on_slaves_main_code  << " line " << __LINE__<< endl;
//             if ( create_log(SLAVEID) ) s_sout << "slaves log count " << iteration_count_on_slaves_main_code  << " line " << __LINE__<< endl;
// 
//             break;
//         }
        
       /* 
        if( !(signaler.get_tag() == TAG_CHANGE_STATUS_TO_CP || signaler.get_tag() == TAG_CHANGE_STATUS_TO_FP)){
            
            signaler.set_tag(TAG_CHANGE_STATUS_TO_FP);
        }*/
//         int sigtag = signaler.get_tag();
//         cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
        
        //DEBON("last_command", last_command)
        
        if ( last_command == TAG_CHANGE_STATUS_TO_CP || last_command == COMMAND_LIST_DO_CP) {
            prev_n_cuts = 0;
            if ( create_log(SLAVEID) ) s_sout << "slaves log count " << iteration_count_on_slaves_main_code  << " line " << __LINE__<< endl;
			DEBON("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ this should not be written" )
            slave_CP ( my_runner, current_lp_optimum_at_slave, allreceivedcuts, prev_n_cuts );
            if ( create_log(SLAVEID) ) s_sout << "slaves log count " << iteration_count_on_slaves_main_code  << " line " << __LINE__<< endl;

        }

        //cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl; 
        // works for original FP only 21,02,2022
        
        if ( last_command == TAG_CHANGE_STATUS_TO_FP || last_command == COMMAND_LIST_DO_OBJ_FP_OR) {
            
            if ( create_log(SLAVEID) ) s_sout << "slaves log count " << iteration_count_on_slaves_main_code  << " line " << __LINE__<< endl;
 				
            //cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
			
				slave_FP_v4 ( my_runner, current_lp_optimum_at_slave);
			if ( create_log(SLAVEID) ) s_sout << "slaves log count " << iteration_count_on_slaves_main_code  << " line " << __LINE__<< endl;
			
        }


        /** NOW all data is sent
         * it is time to WAIT master for signal it will either
         * a) wait for next iteration
         * b) terminate
         * c) run FP()
         *
         * Consider use of barriers at this point
         */
//  		cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
    }
//     cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
	pthread_mutex_lock ( &recv_AC_mutex );
    return 0;
}

int master_cut_selection ( vector<mySerializableRowCut> &biglist, vector<int> &application_status, vector<mySerializableRowCut> &retlist ){

    retlist.clear();
    int retval = 0;
    for ( unsigned i =0; i < application_status.size(); ++i ) {
        if ( application_status[i] == 0 ) {
            retval++;
            retlist.push_back ( biglist[i] );
            application_status[i] = 1;
        }
    }
    return retval;
}

int master_dist_obj_function_list ( vector<My_objective_function> &obj_func_list, vector<int> &slave_types, int number_of_originals, int number_of_random_points_needed ){
    for ( int i = 0; i <number_of_random_points_needed; ++i ) {
        if ( slave_types[i+number_of_originals] == OBJECTIVE_TYPE_ORIGINAL ) {
            CPFP_ERROR ( "trying to send a function to a slave that uses original objective",cout );
            /** will not send any objective function*/
            continue;
        }

        /** THIS IS SEND # (3) in MASTER
         * THIS SEND MATCHES RECV # (3) in SLAVE
         * ******************************/
//         if (obj_func_list[i].is_empty()) {cout << "\n\n\n\n\nsending empty function\n\n\n\n\n" << endl;}
        OBJECT_SEND ( i+number_of_originals+1, TAG_SENDING_OBJECTIVE_FUNCTION, obj_func_list[i] )
    }

    return 0;
}

int master_dist_rowcutlist ( vector<mySerializableRowCut> &cuts_to_be_added,cpfp_executable &runner, int iteration ){
    
//    EASY_SEND_AND_RECEIVE master_mpi_signaler;
    
    if ( iteration == 1 ) {
        /** THIS IS SIGNAL # (signal_2) from master to all slaves
         * SENDING
         * ******************************/
        SIGNAL_ALL_SLAVES_INT ( TAG_SENDING_ROWCUT,0 )
        __MSS ( ++_x_ )

    } else {
        /** current version: adding all cuts that has not been used before */


        int tttl = cuts_to_be_added.size();

        /** THIS IS SIGNAL # (signal_2) from master to all slaves
         * SENDING
         * ******************************/

        SIGNAL_ALL_SLAVES_INT ( TAG_SENDING_ROWCUT,tttl )

        for ( unsigned i = 0; i < cuts_to_be_added.size(); ++i ) {
            /** THIS IS SEND # (4) in MASTER
             * THIS SEND MATCHES RECV # (4) in SLAVE
             * ******************************/
            OBJECT_SEND_ALL_SLAVES ( TAG_SENDING_ROWCUT,cuts_to_be_added[i] )

        }
    }
    return 0;
}

int master_create_and_dist_cut_types(){
    //EASY_SEND_AND_RECEIVE master_mpi_signaler;


    /** THIS IS SIGNAL # (signal_3) from master to all slaves
     * SENDING
     * ******************************/

    Generate_cut_type t = cut_type_generator();


    SIGNAL_ALL_SLAVES_INT ( TAG_SENDING_GENERATE_CUT_TYPE,1 )
    __MSS ( ++_x_ )

    /** currently all slaves generate all cut types*/

    /** THIS IS SEND # (5) in MASTER
     * THIS SEND MATCHES RECV # (5) in SLAVE
     * ******************************/
    OBJECT_SEND_ALL_SLAVES ( TAG_SENDING_GENERATE_CUT_TYPE, t )


    return 0;
}

int master_CP_collect_results ( double &best_inc_obj, double &pre_best_inc_obj, My_solution &best_sol, vector<mySerializableRowCut> &itcuts ){
    stringstream pr;

    /** first wait for signal_5 to get the number of solution each slave will send
     */
    itcuts.clear();
    //EASY_SEND_AND_RECEIVE mpi_signaler;

    for ( int source=1; source < WORLDSIZE; ++source ) {
//         IFMASTER cout << __LINE__ << " source " << source << endl;

        /** THIS IS SIGNAL # (signal_5) from all slaves to master
         * RECEIVING
         * ******************************/
        int flag = -1;
        int signal1 = signaler.receive_signal_i( source, TAG_SENDING_SOLUTION , flag);
        
        //if (flag == TAG_SIGNALER_RETURNS_ERROR) signal1 = TAG_SIGNALER_RETURNS_ERROR;
    
        /** check if the signal is a solution IT SHOULD BE A SOLUTION even if the slave did not provide one*/
        if ( signaler.get_tag() != TAG_SENDING_SOLUTION || flag == TAG_SIGNALER_RETURNS_ERROR) {
            cout << "WARNING: expecting solution with TAG_SENDING_SOLUTION: " << TAG_SENDING_SOLUTION << " but received TAG " << signaler.get_tag() << endl;
        } else {

            if ( signaler.get_number() >0 ) {
                pr << signaler.get_number() << "S ";
                /** that source has a number of solutions to send
                 * as for the first implementation we assume there is only one solution
                 * but to be on the safe side receive the number of solutions stored in
                 * the signaler
                 * */
                for ( int i =0; i < signaler.get_number(); ++i ) {
                    My_solution sol;
                    /** THIS IS RECV # (6) in MASTER
                     * THIS RECV MATCHES SEND # (6) in SLAVE
                     * ******************************/
                    OBJECT_RECV ( source,TAG_SENDING_SOLUTION,sol )

                    if ( sol.get_original_objective() < best_inc_obj ) {
                        pre_best_inc_obj = best_inc_obj;
                        best_inc_obj = sol.get_original_objective();
                        best_sol = sol;
//                         if ( pre_best_inc_obj - best_inc_obj > 0 ) {
//                             /** NOTE:to be replace with something significant */
//                             #ifdef USING_IOPTIMIZE
//                             master_iop.update_rhs ( best_incumbent_objective_at_master );
//                             #endif
//
//                         }
                    }
                }
                /** supposedly received all solutions and updated the best incumbent at hand
                 * */

            }

        }

        /** Now collecting all the cuts from all the slave
         * wait for the signal ROWCUT
         */

        /** THIS IS SIGNAL # (signal_6) from all slaves to master
         * RECEIVING
         * ******************************/
        signaler.receive_signal_i(source, TAG_SENDING_ROWCUT, flag);

    
        if ( signaler.get_tag() != TAG_SENDING_ROWCUT || flag == TAG_SIGNALER_RETURNS_ERROR) {
            cout << "WARNING: expecting solution with TAG_SENDING_ROWCUT: " << TAG_SENDING_ROWCUT << " but received TAG " << signaler.get_tag() << endl;
        } else {
            if ( signaler.get_number() >0 ) {
                pr << signaler.get_number() << "r ";

                /** there is a set of cuts that needs to be received from that slave */
//                 cout << "SLAVE "<< source -1<< " (mast recv) " <<  mpi_signaler.get_number() << " cuts " << endl;
//                 cout << "_____" << endl;
                for ( int i =0; i< signaler.get_number(); ++i ) {
                    mySerializableRowCut cut1;
                    /** THIS IS RECV # (7) in MASTER
                     * THIS RECV MATCHES SEND # (7) in SLAVE
                     * ******************************/
                    OBJECT_RECV ( source,TAG_SENDING_ROWCUT,cut1 )
                    itcuts.push_back ( cut1 );
                }
            }

        }

        /** done for slave = source */

    }


    return 0;

}

int master_create_objective_function_list ( int n_orig, int n_random_p_needs,vector<int> &slave_types, vector<My_objective_function> &obj_func_list, vector<vector< double > > &rand_point_list, vector<double> &current_center_as_v, vector<double> &current_lp_optimum_as_v, bool rest= false ){

    int size_r = rand_point_list.size();
    static int counter = 0;
    if ( rest ) {
        counter = 0;
    }
    obj_func_list.resize ( n_random_p_needs );

    static int i_start = 0;
    static int center_s = -1;
    vector<double> center_substitude;

    if ( center_s <0 ) {
        center_substitude = current_center_as_v;
    } else {
        center_substitude = rand_point_list[center_s];
    }

    for ( int i = 0; i < n_random_p_needs; ++i ) {
        switch ( slave_types[i+n_orig] ) {
        case OBJECTIVE_TYPE_ORIGINAL:
            break;
        case OBJECTIVE_TYPE_RANDOM_FROM_SHORT_DIKIN:
        case OBJECTIVE_TYPE_RANDOM_FROM_LONG_DIKIN:
        case OBJECTIVE_TYPE_RANDOM_FROM_HIT_RUN:

#ifdef USING_IOPTIMIZE
            obj_func_list[i] = My_objective_function ( vector_subtract<double> ( rand_point_list[i_start], center_substitude ) );
            if ( obj_func_list[i].is_empty() ) {
                cout << "---------\n cannot create an objective_function "<< "istart: " << i_start << " counter " << counter << " center_s " << center_s <<   "\n----------" << endl;

                cout << "__________________" <<endl;
                cout << rand_point_list[i_start].size() << endl;
                cout << center_substitude.size() << endl;
                cout << "__________________" <<endl;
                BWAIT
            }

#endif
            i_start++;
            counter++;
            break;
        case OBJECTIVE_TYPE_PERTURBED_FROM_SHORT_DIKIN:
        case OBJECTIVE_TYPE_PERTURBED_FROM_LONG_DIKIN:
        case OBJECTIVE_TYPE_PERTURBED_FROM_HIT_RUN:

#ifdef USING_IOPTIMIZE

            obj_func_list[i] = My_objective_function (
                                   vector_subtract<double> (
                                       vector_multiply_sum<double> (
                                           type_dist_for_CP.get_perturbation(), rand_point_list[i_start],
                                           1-type_dist_for_CP.get_perturbation(),current_lp_optimum_as_v
                                       ),
                                       center_substitude ) );

#endif
            i_start++;
            counter++;
            /** a new objective function such that
             * d_p   = perturbation * d_r + (1-perturbation) * d_o
             *       = perturbation * (x_r - x_ac) + (1-perturbation) * (x*_lp - x_ac)
             *       = perturbation * x_r + (1-perturbation) * x*_lp - x_ac
             *       = {perturbation * x_r + (1-perturbation) * x*_lp} - x_ac
             *       = vector_multiply_sum{p,x_r,(1-p),x_lp} - x_ac
             *       = vector_subtract{ vector_multiply_sum{p,x_r,(1-p),x_lp}, x_ac}
             */
            break;

        default
                :
            break;
        }

        if ( i_start == center_s ) {
            i_start++;
        }

        if ( i_start>=size_r ) {
            i_start = 0;
            center_s ++;
            if ( center_s < size_r ) {
                center_substitude = rand_point_list[center_s];
            } else {
                center_s = -1;
                center_substitude = current_center_as_v;
            }
        }
        if ( big_iteration_change ) {
            center_s = -1;
            i_start = 0;
            big_iteration_change = false;

        }

    }
    return 0;
}

int master_CP ( vector<My_objective_function> &obj_func_list, vector<int> &slave_types, vector<mySerializableRowCut> &allcuts, vector<int> &application_status, vector<mySerializableRowCut> &cuts_to_be_added, My_solution &best_inc_sol,  double &best_inc_at_master,  double &pre_inc_at_master, cpfp_executable &runner, int n_orig, int n_random_p_needs, int iteration_c, ostream &out ){


    //EASY_SEND_AND_RECEIVE master_mpi_signaler;

    vector<mySerializableRowCut> it_row_cuts;
    /** SIGNAL SLAVES SO THAT THEY CONTINUE WITH CP*/

    /** THIS IS SIGNAL # (signal_1) from master to all slaves
     * SENDING
     * ******************************/
    SIGNAL_ALL_SLAVES_INT ( TAG_CHANGE_STATUS_TO_CP,1 )


    /** distribute objective functions */
    if ( !slave_generates_own_functions_for_CP ) {
        master_dist_obj_function_list ( obj_func_list, slave_types, n_orig, n_random_p_needs );
    }


    /** objective functions are distributed */



    /** tell the slave that you are going to send X rowcuts */

    master_dist_rowcutlist ( cuts_to_be_added, runner, iteration_c );

    /** distribute cut types*/
    master_create_and_dist_cut_types();


    /** THIS IS SIGNAL # (signal_4) from master to all slaves
     * SENDING
     * ******************************/
    /** signal best objective function value so far*/

    SIGNAL_ALL_SLAVES_DOUBLE ( TAG_SENDING_OBJECTIVE_VALUE,best_inc_at_master )


    /** everything has beed sent
     * NOW WAIT FOR EACH SLAVE TO COMPLETE AND COLLECT RESULTS
     *
     */





    master_CP_collect_results ( best_inc_at_master, pre_inc_at_master, best_inc_sol ,it_row_cuts );



    /** now all results are collected from all slaves
     * it is time to
     * a) evaluate current cuts and add then to cut pool
     * b) evaluate cut pool as a system
     * c) NOTE: iop_already updated in terms of rhs
     * d)
     */

    /** sort the cuts to eliminate equals */

    sort ( it_row_cuts.begin(),it_row_cuts.end() );
    std::vector<mySerializableRowCut>::iterator it;

    it = unique ( it_row_cuts.begin(),it_row_cuts.end() );
    it_row_cuts.resize ( std::distance ( it_row_cuts.begin(),it ) );

    /** iteration cuts are inserted to all cuts */
    vector<int> temp ( it_row_cuts.size(),0 );
    allcuts.insert ( allcuts.end(),it_row_cuts.begin(),it_row_cuts.end() );
    application_status.insert ( application_status.end(),temp.begin(),temp.end() );

    return 0;
}

int master_FP_inner (My_solution &best_inc_sol, double &best_inc_at_master,  double &pre_inc_at_master, double &current_LB, ostream &m_out = std::cout ){

    vector<double> objective_values_from_slaves;
    objective_values_from_slaves.resize ( WORLDSIZE,1e31 );
    bool continue_fp_flag = true;
    Timer time_checker;
    time_checker.start();
    static double latest_printed_time = 0;
    latest_printed_time = 0;
    int local_counter = 0;
    double latest_update_for_slaves = best_inc_at_master;
    FP_iterate_decider_master ( 0,0,0 );
	double best_obj_from_slaves = 1e31;
	
//     min_iterlim_master_checker ( 0);
    while ( continue_fp_flag ) {
        local_counter++;
        //if (local_counter%10 == 0) cout << "local_counter "<< local_counter << endl;
        /**probe all slaves for new incumbent*/
        //bool obj_updated = false;
        for ( int i=1; i < WORLDSIZE; ++i ) {
//             cout << "checking slave " << i << endl;
            int flag = -99;
            PROBE(i,FOUND_A_FEASIBLE_SOLUTION_IN_FP)
			//if (flag == TAG_SIGNALER_PROBE_SUCCESSFUL) {DEBON("flag:", flag)}
			if ( flag == TAG_SIGNALER_PROBE_SUCCESSFUL ) {
//                                 cout << "flag == 1" << endl;
 				//cout << "receiveing from slave " << i << endl;
				//cout << "reprobe" << endl;
				//PROBE(i,FOUND_A_FEASIBLE_SOLUTION_IN_FP)
				//cout << "flag = "<< flag << endl;
				//cout << signaler.get_double() << endl;
				
                double r_d = 1e99;
                RECV_D(i, FOUND_A_FEASIBLE_SOLUTION_IN_FP);
				//cout << "received value " << r_d << endl;
                objective_values_from_slaves[i] = r_d;


//                 }

					//obj_updated = true;					
// 				cout << "slave " << i -1 << " sends " << zik <<  " time passed " << time_checker.stop() << " best was " << best_obj_from_slaves << endl;
					
                if ( best_obj_from_slaves > r_d ) {
                    My_solution asd;

//               
                    OBJECT_RECV(i,TAG_SENDING_SOLUTION, asd);
// 					cout <<" --- UPDATING best_obj_from_slaves (" << i-1 <<") " << endl; 
// 					cout << best_obj_from_slaves <<" " << zik << endl;
                    best_obj_from_slaves = r_d;
// 					cout << best_obj_from_slaves <<" " << zik << endl;
// 					master_ub_for_iop = best_obj_from_slaves;
					
					// 					pre_inc_at_master = best_inc_at_master;
// 					best_inc_at_master = best_obj_from_slaves;
                    best_inc_sol = asd;
					m_out << "latest_objective : " << setprecision ( 10 ) << best_obj_from_slaves <<" " << r_d <<  " time passed__  " << time_checker.stop() << endl;
					cout << "latest_objective : " << setprecision ( 10 ) << best_obj_from_slaves <<" " << r_d <<  " time passed " << time_checker.stop() << endl;
					
// 					if (best_obj_from_slaves < current_LB){
// 						continue_fp_flag =false;
// 					}
//                     cout <<" a feasible solution from slave " << i << " with value " << zik << " is received " <<endl;
// 					cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
                }
//                 cout <<" a feasible solution from slave " << i << " with value " << zik << " is received " <<endl;
// 				cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;


            }
            ////// t++;

        }
        if ( ( latest_update_for_slaves >best_obj_from_slaves ) ) {

//             cout <<" broadcasting new objective value "<< best_obj_from_slaves<<endl;

            /** tell all slaves about the situation */
            pre_inc_at_master = best_inc_at_master;
            best_inc_at_master = best_obj_from_slaves;
			 #ifdef USING_IOPTIMIZE 
            master_ub_for_iop = best_obj_from_slaves;
            #endif
            latest_update_for_slaves = best_obj_from_slaves;
			
			#ifdef USING_MPI
            	for ( int i =1; i < WORLDSIZE; ++i ) {
                	MPI_Send ( &latest_update_for_slaves, 1, MPI_DOUBLE, i ,TAG_SENDING_OBJECTIVE_VALUE_IN_FP, MPI_COMM_WORLD );
            	}
            
            #else
 				signaler.send_signal_d(DESTINATION_ALL_SLAVES,TAG_SENDING_OBJECTIVE_VALUE_IN_FP,latest_update_for_slaves);
       		#endif

            if ( stop_fp_when_first_solution_is_found ) {
                continue_fp_flag =false;
            }
            else{
//                 cout << " calling  master_all_iop_related_function_calls(0);" << endl;
                master_all_iop_related_function_calls(0);
            }
        } 
        else {
            usleep ( SLEEP_TIME_BETWEEN_FP_PROBES );
        }

        if ( ( time_checker.stop() - latest_printed_time ) >1 ) {
			m_out << "latest_objective : " << setprecision ( 10 ) << latest_update_for_slaves <<  " "<< best_inc_at_master << " time passed " << floor(time_checker.stop()) << endl;
			cout << "latest_objective : " << setprecision ( 10 ) << latest_update_for_slaves <<  " "<< best_inc_at_master << " time passed " << floor(time_checker.stop()) << endl;
            latest_printed_time = time_checker.stop();
//             for(unsigned i = 1; i < slave_iteration_numbers.size();++i){
// 				cout <<" ["<<             miniterlim_checker[i]<<
// 				 "]:"<< slave_iteration_numbers[i];
//             }
//             cout << endl;



        }
        /** do I need to send new objective function set */
//         if (time_checker.stop() >= inner_time_lim) continue_fp_flag = false;
        int n_slaves_done = FP_iterate_decider_master ( 1,1,1 );

//         cout << "n_slaves_done: " << n_slaves_done<< endl;

        if ( n_slaves_done == N_SLAVES ) {
//         if ( min_iterlim_master_checker ( 1) == MPI::COMM_WORLD.Get_size() - 1 ) {
            /** all slaves must have been done*/

			cout << " all slaves must have been done" <<endl;
			m_out << " all slaves must have been done" <<endl;
            FP_iterate_decider_master ( 2,2,2 );
//             min_iterlim_master_checker( 2 );
            continue_fp_flag = false;
			cout << "latest: " ;
			m_out << "latest: " ;
            for ( unsigned i = 1; i < slave_iteration_numbers.size(); ++i ) {
				cout << " ["<<miniterlim_checker[i]<< "]:"<< slave_iteration_numbers[i];
				m_out << " ["<<miniterlim_checker[i]<< "]:"<< slave_iteration_numbers[i];
            }
            cout << endl;
			m_out << endl;
        }
        
        
    }

    int done_with_FP = TAG_CHANGE_STATUS_NEXT_FP_ITERATION;

    #ifdef USING_MPI
    for ( int i =1; i < WORLDSIZE; ++i ) {
        MPI_Send ( &done_with_FP,1,MPI_INT,i,TAG_CHANGE_STATUS_NEXT_FP_ITERATION,MPI_COMM_WORLD );
    }
    #else 
     //signaler.send_signal_i(0,TAG_CHANGE_STATUS_NEXT_FP_ITERATION,done_with_FP);
       #endif

    return 0;

}



int master_FP_outer ( vector<My_objective_function> &obj_func_list, vector<int> &slave_types, My_solution &best_inc_sol, vector<double> &current_lp_optimum_as_v, double &best_inc_at_master,  double &pre_inc_at_master,  double &current_LB, int n_orig, int n_random_p_needs, ostream &m_out = std::cout ){

	static int outer_counter = 0;
	outer_counter++;
	static double latest_sent_ub = 1e98;
    //EASY_SEND_AND_RECEIVE master_mpi_signaler;

	

    /** SIGNAL SLAVES SO THAT THEY CONTINUE WITH FP*/

    /** THIS IS SIGNAL # (signal_10) from master to all slaves
     * SENDING
     * ******************************/
	for (int i = 1; i <WORLDSIZE; ++i){
		if (latest_slave_commands[i] !=COMMAND_LIST_DO_OBJ_FP_OR){
			latest_slave_commands[i] =COMMAND_LIST_DO_OBJ_FP_OR;
			signaler.send_signal_i(i,TAG_SIGNAL_COMMAND, COMMAND_LIST_DO_OBJ_FP_OR);
		}
	}
	//SIGNAL_ALL_SLAVES_INT ( TAG_SIGNAL_COMMAND, COMMAND_LIST_DO_OBJ_FP_OR )


    /** THIS IS SIGNAL # (signal_14) from master to all slaves
     * SENDING
     * ******************************/
    /** signal best objective function value so far*/
	if (double_significantly_less(best_inc_at_master,latest_sent_ub)){
		latest_sent_ub = best_inc_at_master;
		SIGNAL_ALL_SLAVES_DOUBLE ( TAG_SENDING_OBJECTIVE_VALUE_IN_FP,best_inc_at_master )
	}

// 	cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;

    /** distribute objective functions */

    if ( !slave_generates_own_functions_for_FP ) {
        master_dist_obj_function_list ( obj_func_list,slave_types,n_orig,n_random_p_needs );
    }
    /** objective functions are distributed */


    /** tell the slave that you are going to send X rowcuts */
    // master_dist_rowcutlist ( cuts_to_be_added, runner, iteration_c);

    /** distribute FP types*/
//     master_create_and_dist_FP_types();




cout << "starting parallel FP for the outer_counter="<< outer_counter << "th time" << endl;
// cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;

    master_FP_inner (best_inc_sol, best_inc_at_master, pre_inc_at_master, current_LB, m_out );




    return 0;
}


bool what_to_run_now ( double gap_n = 0,double gap_o = 0 ){
    static bool latest_run_CP = true;
    if ( !latest_run_CP ) {
        latest_run_CP = true;
        return latest_run_CP;
    }
    double plus = 0;
    if ( double_equality ( gap_o,0 ) ) {
        plus = 1;
    }

    if ( ( gap_o - gap_n ) / ( gap_o + plus ) > 0.01 ) {
        cout << " gap change " <<100* ( gap_o - gap_n ) / ( ( gap_o + plus ) ) << "% continue with CP " <<endl;
        return latest_run_CP ;
    } else {
        cout << " gap change " <<100* ( gap_o - gap_n ) / ( ( gap_o + plus ) ) << "% continue with FP " <<endl;
        latest_run_CP = false;
        return latest_run_CP;
    }

}

int output_id ( int i2=0 ){
    static int int_i = 0;
    bool ftemp = true;
    while ( ftemp ) {
        int_i++;
        stringstream ss1;
        ss1 << "detailed_results/Exp_"<<experiment_number<<"_setting_" <<"/"<< p_name <<"_"<< type_dist_for_CP.string_print() <<"_il"<< number_of_iterations << "_tl" << timelimit<<"_" << int_i;
        string dds = ss1.str();
        ftemp = file_exists ( dds );
    }
    return ( int_i - i2 );
}

int master_send_lb ( double LB_for_original_objective ){
    //EASY_SEND_AND_RECEIVE master_mpi_signaler;
	static double latest_sent_lb = -1e99;
	if (double_significantly_less(latest_sent_lb,LB_for_original_objective)){
		latest_sent_lb = LB_for_original_objective;
		SIGNAL_ALL_SLAVES_DOUBLE ( TAG_SENDING_LOWERBOUND,LB_for_original_objective );
	}
	
    /** THIS IS SEND # (1) in MASTER
     * THIS SEND MATCHES RECV # (1) in SLAVE
     * ******************************/
// 	cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
	//SIGNAL_ALL_SLAVES_DOUBLE ( TAG_SENDING_LOWERBOUND,LB_for_original_objective );
// 	cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
    __MSD ( ++_x_ )

    return 0;

}

int master_send_current_lp_optimum ( My_solution &arg_lp_opt ){
	
	if (latest_lp_opt_sent_to_slaves == arg_lp_opt){
		return 1;
	}
	else{
		OBJECT_SEND_ALL_SLAVES(TAG_SENDING_SOLUTION,arg_lp_opt)
	}
	/*
	#ifdef USING_MPI
    for ( int i = 1; i <WORLDSIZE; ++i ) {
        req_stat_vectors r;
// 		cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
        arg_lp_opt._send_mpi ( i, TAG_SENDING_SOLUTION,r.v_req, r.v_stat );
// 		cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
    }
    #else 
    arg_lp_opt._send_alt1(DESTINATION_ALL_SLAVES,TAG_SENDING_SOLUTION,signaler);
    #endif*/
    return 0;

}


int master_run ( int argc, char* argv[], string argfilename ="__NULL" ){
    int master_FP_counter = 0;
    int master_CP_counter = 0;
    double last_written_FP_obj = 1e31;
    double last_written_CP_lb = -1e31;
    stringstream mname;
    filename_generator ( mname, MACRO_SLAVE_ID );
    cout << "Master " << MACRO_SLAVE_ID << " " << mname.str() << endl;

	latest_slave_commands.resize(WORLDSIZE,0);
	

	
	stringstream ss1;
	ss1 << "detailed_results/Exp_" <<  experiment_number<<"_setting_"<<setting_number;
	// 		mkdir(ss1.str().c_str(), S_IWOTH || S_IROTH); 
	mkdir(ss1.str().c_str(), ACCESSPERMS); 
	
	ss1<<"/"<< p_name <<   "_out.txt";
	
	ofstream master_cout ( ss1.str().c_str(),ios_base::app );
    timestamp ( master_cout );
    master_cout <<endl;
	
	type_dist_for_FP.print ( master_cout );
	type_dist_for_CP.print ( master_cout );
	
    master_cout << "time limit: " << timelimit << endl;
    master_cout << "iter limit: " << number_of_iterations << endl;

    unsigned prev_cut_Size = 0;
    int no_new_cut_counter = 0;

    double best_incumbent_objective_at_master = starting_objective;
    double previous_incumbent_objective_at_master = starting_objective;

    cpfp_executable my_runner ( argfilename );
    my_runner.set_seed ( my_default_seed );

    if ( my_runner.initialize(solver_type,ipopt_solver) == 1 ) {
		objective_improvement_coefficient_FP = 1;
		objective_improvement_coefficient_CP = 1;
		// galiba Cout da yazanlar yanlis
        //cout << " all variables integer or all continious variables has zero objective function coefficient." << endl;
		//master_cout << " all variables integer or all continious variables has zero objective function coefficient." << endl;
    }
    double sb = my_runner.calculate_pseudo_bound();
    my_runner.set_pseudo_bounds(sb);
	read_FP_parameters_from_file(my_runner);
//     double objective_imp = 0;
//     if ( USE_RELATIVE_IMPROVEMENT ) {
// //     cout << " grro " << __LINE__ << endl;
//
//         objective_imp = ( best_incumbent_objective_at_master - my_runner.get_objective_cutoff_constraint_lb() ) * objective_improvement_percentage;
// 		if (objective_imp < objective_improvement_coefficient) objective_imp = objective_improvement_coefficient;
//
//     } else {
//         objective_imp = objective_improvement_coefficient;
//
//     }
//
//     objective_imp = 0;
//     my_runner.update_objective_cutoff_constraint_ub ( best_incumbent_objective_at_master -objective_imp );



	cout << "objective_improvement_coefficient_FP: " << objective_improvement_coefficient_FP <<endl;
	master_cout << "objective_improvement_coefficient_FP: " << objective_improvement_coefficient_FP <<endl;
	cout << "objective_improvement_coefficient_CP: " << objective_improvement_coefficient_CP <<endl;
	master_cout << "objective_improvement_coefficient_CP: " << objective_improvement_coefficient_CP <<endl;
	
	
    time_to_optimize_in_master.start();
    int ref = my_runner.optimize();
    time_to_optimize_in_master.pause();
    cout << "ref: "<< ref << endl;




//     my_runner.print_bounds();
//     usleep(10000000);


    if(ref != 1){
         cout << "ref at master level: "<< ref << endl;
//         usleep(10000000);
    }

    double LB_for_original_objective = my_runner.getObjValue();
    double starting_LB = LB_for_original_objective;

    summarystring << " i_LB: " << std::setprecision ( 10 ) << LB_for_original_objective << " ";
	summarystring2 << " i_LB: " << std::setprecision ( 10 ) << LB_for_original_objective << " ";
	master_cout << " i_LB: " << std::setprecision ( 10 ) << LB_for_original_objective << " ";
	
	cout << " i_LB: " << std::setprecision ( 10 ) << LB_for_original_objective << " " << endl;;
	master_cout << " i_LB: " << std::setprecision ( 10 ) << LB_for_original_objective << " " << endl;;

    const double *current_lp_optimum  = my_runner.getColSolution();
    int ncols = my_runner.get_ncols();

    vector<double> current_lp_optimum_as_vector;
    current_lp_optimum_as_vector.resize ( ncols );
    for ( int i = 0; i < ncols; ++i ) {
        current_lp_optimum_as_vector[i]=current_lp_optimum[i];
    }
    //EASY_SEND_AND_RECEIVE master_mpi_signaler;

    master_send_lb ( LB_for_original_objective );

#ifdef USING_IOPTIMIZE
    master_iop.set_seed ( my_default_seed );
    master_iop.setCenProbType ( solve_original_centering_problem );
    master_iop.initialize2 ( argfilename );
    master_iop.iop_set_seed ( my_default_seed );

#endif

    if ( type > 3 ) {
#ifdef USING_IOPTIMIZE
        master_iop.set_walk_type ( type-3 );
#endif

    } else {
#ifdef USING_IOPTIMIZE
        master_iop.set_walk_type ( type );
#endif
    }

    int number_of_random_points_needed_CP = type_dist_for_CP.get_number_of_instances_that_need_AC();
    int number_of_random_points_needed_FP = type_dist_for_FP.get_number_of_instances_that_need_AC();
    int number_of_originals_CP = type_dist_for_CP.get_number_of_instances_that_does_not_need_AC();
    int number_of_originals_FP = type_dist_for_FP.get_number_of_instances_that_does_not_need_AC();

    vector<vector<double> > random_point_list;
    int iteration_limit = number_of_iterations;

    vector<mySerializableRowCut> allrowcuts;
    vector<int> cutapplied; /** 0 means cut not applied but applicable
							  * 1 means cut applied
							  * -1 means cut not applied and will not be considered for application later
							*/
//     vector<My_objective_function> objective_function_list;
// 	for (int iteration_counter =0; iteration_counter < iteration_limit; ++ iteration_counter)
    bool masterflagcontinue =true;
    int iteration_counter =0;

    vector<mySerializableRowCut> iteration_row_cuts;
    vector<mySerializableRowCut> cuts_to_be_added;
    Timer time_checker;
//     Timer swap_timer;
    time_checker.start();
//     swap_timer.start();
    double gap_before = 1e30 - LB_for_original_objective;
    double gap_now = gap_before;
    bool do_CP = false;
    vector<int> types_per_slave_CP = type_dist_for_CP.get_type_per_slave();
    vector<int> types_per_slave_FP = type_dist_for_FP.get_type_per_slave();

    time_to_optimize_in_master.restart();
    my_runner.optimize();
    time_to_optimize_in_master.pause();
    My_solution lp_solution ( ncols, current_lp_optimum );
    lp_solution = My_solution ( my_runner.getColSolution_asVector() );
//     cout << __LINE__ << " --" << endl;
    LB_for_original_objective = my_runner.getObjValue();
#ifdef USING_IOPTIMIZE
    master_lb_for_iop = LB_for_original_objective;
#endif
    vector<double> current_center_as_v;

    if ( number_of_random_points_needed_CP+number_of_random_points_needed_FP>0 ) {
// 		cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
        if (!AC_v2) master_all_iop_related_function_calls ( 0,true );
// 		cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
    }

    while ( masterflagcontinue ) {
        iteration_counter++;
        cout << "iteration "<< iteration_counter <<endl;
        iteration_row_cuts.clear();
//         /** create random points */


//         cout << " LINE: " << __LINE__ <<endl;
// usleep(1000000);
        master_send_lb ( LB_for_original_objective );
//         cout << " LINE: " << __LINE__ <<endl;
// 		usleep(1000000);
// 		lp_solution.print ();
        master_send_current_lp_optimum ( lp_solution );
// 		cout << " LINE: " << __LINE__ <<endl; 
// 		usleep(1000000);

        /** FIX:  most probably fixed*/
        /** WARNING: size is a problem in the following for loop
         * MPI::COMM_WORLD.Get_rank 0 (master)
         *                          1 (slave 1 with original obj)
         *                          2 (slave 2 random obj 1)
         *                          3 (slave 3 random obj 2)
         * type_dist.get_type_per_slave
         *                          0 (slave 1 with original obj)
         *                          1 (slave 2 random obj 1)
         *                          2 (slave 3 random obj 2)
         * Objective_function_list
         *                          0 (slave 2 random obj 1)
         *                          1 (slave 3 random obj 2)
         */
        /** chance random points to objective directions */
        vector<My_objective_function> objective_function_list;

        stringstream s2;


//         vector<double> current_center_as_v = master_iop.get_current_center_point();


//         /** OBJECTIVE FUNCTIONS ARE CREATED*/

        // do_CP = true;
//         if(do_CP){ /** last call was CP */
//             if(time_for_CP.stop() < CP_time_limit){ /** time is NOT up continue with CP*/
//             }
//             else{
//                 do_CP = false;
//             }
//             if(no_new_cut_counter >=5 ) {
//
//                 do_CP = false;
//             }
//
//         }
//         else{ /**last call was FP*/
//             /** change to CP*/
//             do_CP = true;
//             time_for_CP.start();
//         }
//
//
//         if (iteration_counter == 1) do_CP = false;
//
//         //do_CP = do_CP_given_iter_num(iteration_counter);
        double temp_double = best_incumbent_objective_at_master;
        string CP_FP_identifier;

//         cout << "do_CP: " << do_CP <<endl;
//         cout << " LINE: " << __LINE__ <<endl;
//
//         do_CP = CP_FP_decider();

// 		        cout << " LINE: " << __LINE__ <<endl;
		int do_CP_FP = 1;
		//int do_CP_FP = CP_FP_decider_int();
// 		
        cout << "do_CP_FP: " << do_CP_FP <<endl;
		
        if ( do_CP_FP == 2 ) {
             cout << " doing CP " << endl;
            CP_FP_identifier = "(CP)";

            if ( ( !slave_generates_own_functions_for_CP ) && ( number_of_random_points_needed_CP>0 ) ) {
                /** */
                if (!AC_v2)master_all_iop_related_function_calls ( 3,true );
                time_to_create_RW_points_master.restart();
#ifdef USING_IOPTIMIZE 
                current_center_as_v = master_iop.get_current_center_point();

                master_iop.iop_set_seed(my_default_seed);
                master_iop.init_sampler();
                int r_res = master_iop.get_next_n_random_point ( random_point_list,number_of_random_points_needed_CP );
                if ( r_res != 0 ) {
                    cout <<" failed to generate "<< number_of_random_points_needed_CP <<" random points. res: "<< r_res<<endl;
                    master_cout <<" failed to generate "<< number_of_random_points_needed_CP <<" random points. res: "<< r_res <<  endl;
//                     master_iop.set_seed ( my_default_seed );
                    cout << master_iop.get_next_n_random_point ( random_point_list,number_of_random_points_needed_CP ) << endl;
                }
                time_to_create_RW_points_master.pause();
				if (!AC_v2)master_all_iop_related_function_calls ( 4,true );
                master_create_objective_function_list ( number_of_originals_CP,number_of_random_points_needed_CP,types_per_slave_CP,objective_function_list,random_point_list,current_center_as_v, current_lp_optimum_as_vector,false );
#endif
            } else if ( number_of_random_points_needed_CP>0 ) {
                if ( total_number_of_ACs_sent_by_master < total_number_of_ACs_computed_in_master ) {
					if (!AC_v2)master_all_iop_related_function_calls ( 2 );
                }
            }

            // cout << "running a round of CP "<< endl;
            int CP_res = master_CP ( objective_function_list, types_per_slave_CP, allrowcuts,cutapplied, cuts_to_be_added, best_incumbent_solution, best_incumbent_objective_at_master, previous_incumbent_objective_at_master, my_runner, number_of_originals_CP, number_of_random_points_needed_CP, iteration_counter, master_cout );
            LB_for_original_objective = my_runner.getObjValue();
            lp_solution = My_solution ( my_runner.getColSolution_asVector() );

            if ( CP_res == 1 ) {
                cout << "breaking at line "<< __LINE__ << endl;
                SIGNAL_ALL_SLAVES_INT ( -1,-1 )
                return 1;
            }
            if ( temp_double - best_incumbent_objective_at_master > 0 ) {
//                cout << "updating UB to "<< best_incumbent_objective_at_master - objective_improvement_coefficient << endl;
#ifdef USING_IOPTIMIZE
                master_ub_for_iop = best_incumbent_objective_at_master;
#endif
				UPDATE_RUNNER_CP(best_incumbent_objective_at_master,LB_for_original_objective,my_runner)
   
            }
            int tttl = master_cut_selection ( allrowcuts,cutapplied,cuts_to_be_added );
            if ( tttl >0 ) {
                time_to_apply_cuts_master.restart();
                my_runner.apply_cuts ( cuts_to_be_added );
                time_to_apply_cuts_master.pause();
                time_to_optimize_in_master.restart();
                my_runner.optimize() ;
                time_to_optimize_in_master.pause();
                LB_for_original_objective= my_runner.getObjValue();
#ifdef USING_IOPTIMIZE 
                master_lb_for_iop = LB_for_original_objective;
#endif
                lp_solution = My_solution ( my_runner.getColSolution_asVector() );

            }




        } 
        else if (do_CP_FP == 1){
            cout << " doing FP " << endl;

            CP_FP_identifier = "(_FP)";

            master_FP_counter ++;
            if ( last_written_CP_lb <= LB_for_original_objective + objective_improvement_coefficient_CP ) {
                summarystring << "CP" << master_CP_counter << " : " << std::setprecision ( 10 ) <<LB_for_original_objective << " ";
                last_written_CP_lb = LB_for_original_objective;
            } else {
                summarystring << "CP" <<std::setprecision ( 10 ) << master_CP_counter << " : . ";

            }
//             cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
            if ( ( !slave_generates_own_functions_for_FP ) && ( number_of_random_points_needed_FP>0 ) ) {
 				cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
				if (!AC_v2)master_all_iop_related_function_calls ( 3,true );
                time_to_create_RW_points_master.restart();
#ifdef USING_IOPTIMIZE 
                current_center_as_v = master_iop.get_current_center_point();
                master_iop.set_seed ( my_default_seed );

                int r_res = master_iop.get_next_n_random_point ( random_point_list,number_of_random_points_needed_FP );

                if ( r_res != 0 ) {

                    cout <<" failed to generate "<< number_of_random_points_needed_FP <<" random points. res: "<< r_res<<endl;
                    master_cout <<" failed to generate "<< number_of_random_points_needed_FP <<" random points. res: "<< r_res <<  endl;
                    master_iop.set_seed ( my_default_seed );
                    cout << master_iop.get_next_n_random_point ( random_point_list,number_of_random_points_needed_FP ) << endl;
                }
#endif
                time_to_create_RW_points_master.pause();

				if (!AC_v2)master_all_iop_related_function_calls ( 4,true );

                master_create_objective_function_list ( number_of_originals_FP,number_of_random_points_needed_FP,types_per_slave_FP,objective_function_list,random_point_list,current_center_as_v, current_lp_optimum_as_vector,false );

            } else if ( number_of_random_points_needed_FP>0 ) {
// 				cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
                if ( total_number_of_ACs_sent_by_master < total_number_of_ACs_computed_in_master ) {
					if (!AC_v2)master_all_iop_related_function_calls ( 2 );
                }
            }


//             cout << "running a round of FP "<< endl;
            master_FP_outer ( objective_function_list, types_per_slave_FP, best_incumbent_solution, current_lp_optimum_as_vector, best_incumbent_objective_at_master, previous_incumbent_objective_at_master, LB_for_original_objective, number_of_originals_FP, number_of_random_points_needed_FP, master_cout );
            if ( last_written_FP_obj > best_incumbent_objective_at_master ) {
                summarystring << "FP" << master_FP_counter <<" : " <<std::setprecision ( 10 ) << best_incumbent_objective_at_master<< " ";
                last_written_FP_obj = best_incumbent_objective_at_master;
                if ( master_FP_counter == 1 ) {
                    summarystring2 << "FP1: " <<std::setprecision ( 10 ) << best_incumbent_objective_at_master;
                }
            } else {
                summarystring << "FP" << master_FP_counter <<" : . ";

            }
//             cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
            master_CP_counter++;

        }
        gap_before = gap_now;
        gap_now = best_incumbent_objective_at_master - LB_for_original_objective;

        double gap_percentage = abs ( ( gap_now ) /LB_for_original_objective );
        if ( iteration_counter == 0.5 ) {
            cout << "iteration  "<< iteration_counter << "\t\t\tstarting objective:\t" << best_incumbent_objective_at_master<<  "\tlb: " << LB_for_original_objective << "\t";
            global_time.printTime ( global_time.stop() );
            master_cout << "iteration  "<< iteration_counter << "\t\t\tstarting objective:\t" << best_incumbent_objective_at_master<<  "\tlb: " << LB_for_original_objective<< "\t";
            global_time.printTime ( global_time.stop(),master_cout );
        } 
        else {
            cout << "iteration  "<< iteration_counter << CP_FP_identifier.c_str() <<" END\tallrowcuts.size(): "<< allrowcuts.size() << "\tbset sol:\t" << std::setprecision ( 10 ) << best_incumbent_objective_at_master<< " \tlb: " << LB_for_original_objective<<  "\tgap: " << setprecision ( 4 ) << 100*gap_percentage <<"% \ttime: " ;

            global_time.printTime ( global_time.stop() );
            master_cout << "iteration  "<< iteration_counter << CP_FP_identifier.c_str() << " END\tallrowcuts.size(): "<< allrowcuts.size() << "\tbset sol:\t" << std::setprecision ( 10 ) <<best_incumbent_objective_at_master<< " \tlb: " << LB_for_original_objective<<  "\tgap: " << setprecision ( 4 ) << 100*gap_percentage <<"% \ttime: " ;
            global_time.printTime ( global_time.stop(),master_cout );
        }
        if ( ( do_CP ) && ( prev_cut_Size == allrowcuts.size() ) ) {
            no_new_cut_counter++;
        } else {
            no_new_cut_counter =0;
        }
#ifdef USING_IOPTIMIZE 
        master_ub_for_iop = best_incumbent_objective_at_master;
        master_lb_for_iop = LB_for_original_objective;
#endif
        
		//cout << " callling AC =" << __LINE__ << endl;
		
        if(!AC_v2)master_all_iop_related_function_calls ( 0 ); /** start calculating AC */
        //cout << " calculate AC must have been called" << endl;

        prev_cut_Size = allrowcuts.size();

        if ( ( iteration_counter>= iteration_limit ) && ( iteration_limit>0 ) ) {
            cout << "breaking at line "<< __LINE__ << endl;

            break;
        }
        if ( time_checker.stop() >=timelimit ) {
            cout << " TIME LIMIT REACHED " << endl;
            master_cout << " TIME LIMIT REACHED " << endl;

            break;
        }


        if ( no_new_cut_counter >=5 ) {
        }

        if ( LB_for_original_objective > ( best_incumbent_objective_at_master - min(objective_improvement_coefficient_CP, objective_improvement_coefficient_FP ) )) {

            break;
        }

        

    }

    master_send_lb ( LB_for_original_objective );
    master_send_current_lp_optimum ( lp_solution );

    cout << "iteration  "<< iteration_counter << " END\tallrowcuts.size(): "<< allrowcuts.size() << " bset sol:\t" << std::setprecision ( 10 ) <<best_incumbent_objective_at_master<< "\tlb: " << LB_for_original_objective<<endl;

    master_cout << "iteration  "<< iteration_counter << " END\tallrowcuts.size(): "<< allrowcuts.size() << " bset sol:\t" << std::setprecision ( 10 ) << best_incumbent_objective_at_master<< "\tlb: " << LB_for_original_objective<<endl;


    prev_cut_Size = allrowcuts.size();
    SIGNAL_ALL_SLAVES_INT ( -1,-1 );
    global_time.stop();
//     cout << "CP  time:"; global_time.printTime(global_time.stop());

    if ( best_incumbent_objective_at_master < LB_for_original_objective + min(objective_improvement_coefficient_FP, objective_improvement_coefficient_CP)) {
        cout << "An optimal solution is proven..." << endl;
        master_cout << "An optimal solution is proven..." << endl;

    } else {

    }

    p_name = split_string ( argfilename,'/' ).back();

    cout << argfilename <<" st_obj: "<< starting_objective << " st_LB: "<<starting_LB<<" iter: " << iteration_counter <<" n_cuts: " << allrowcuts.size() << " final_obj: " << best_incumbent_objective_at_master << " final_LB: " << LB_for_original_objective << " time: " ;

    type_dist_for_CP.shortlineprint();
    global_time.printTime ( -1 );

    master_cout << argfilename <<" st_obj: "<< starting_objective << " st_LB: "<<starting_LB<<" iter: " << iteration_counter <<" n_cuts: " << allrowcuts.size() << " final_obj: " << best_incumbent_objective_at_master << " final_LB: " << LB_for_original_objective << " time: " ;
    type_dist_for_CP.shortlineprint ( master_cout );
    global_time.printTime ( -1,master_cout );



    cout << p_name <<"\t"<< starting_objective << "\t"<<starting_LB<<"\t" << iteration_counter <<"\t" << allrowcuts.size() << "\t" << best_incumbent_objective_at_master << "\t" << LB_for_original_objective << "\t" ;
    type_dist_for_CP.shortlineprint();


    ofstream out2 ( "output_of_pcpfp.txt",ofstream::app );

    out2 <<  p_name <<"\t"<< starting_objective << "\t"<<starting_LB<<"\t" << iteration_counter <<"\t" << allrowcuts.size() << "\t" << best_incumbent_objective_at_master << "\t" << LB_for_original_objective << "\t" ;
    type_dist_for_CP.shortlineprint ( out2 );
    global_time.printTime ( -1,out2 );


//     string filename2 ;
//     next_file_ ( filename2,true );
//     cout << " mps file to be written on " << filename2.c_str() << endl;
    //my_runner.puke_mps_file(filename2.c_str());
//     cout << endl<< " ---- " << endl;
//     cout << endl<< " ---- " << endl;
//     my_runner.print_solution_wth_column_types(best_incumbent_solution);
//     cout << endl<< " ---- " << endl;
//     cout << endl<< " ---- " << endl;

//     best_incumbent_solution.print();
//     cout << " ---- " << endl;


    summarystring << " F_obj: " << std::setprecision ( 10 ) <<best_incumbent_objective_at_master << " F_lb: " << LB_for_original_objective ;
	
    summarystring2 << " F_obj: " << std::setprecision ( 10 ) <<best_incumbent_objective_at_master << " F_lb: " << LB_for_original_objective;

	stringstream esy;
	esy << "easycopy_" << experiment_number << ".txt";

	ofstream oe(esy.str().c_str(),ios_base::app);

	oe << p_name << "\t" << best_incumbent_objective_at_master << endl;

      return 0;
}

int post_run(){

    double avg_time_to_optimize = get_average_time(time_to_optimize_in_slave);
    double avg_time_to_generate_cuts = get_average_time(time_to_generate_cuts);
    double avg_time_to_apply_cuts = get_average_time(time_to_apply_cuts_slave);
	double avg_time_for_l1_norm = get_average_time(l1_norm_timer);
	double avg_time_for_RW_point_generation = get_average_time(time_to_create_RW_points_slave);
	double avg_time_for_generating_FP_start = get_average_time(time_to_create_starting_points_for_FP_slave);
	
	
	
    IFMASTER {
// 		stringstream sss;
// 		int int_i = output_id ( 1 );
		
// 		sss << "detailed_results/" << p_name <<"_"<< type_dist_for_CP.string_print() <<"_il"<< number_of_iterations << "_tl" << timelimit<<"_"<< int_i;
		
		stringstream ss1;
		ss1 << "detailed_results/Exp_" <<  experiment_number<<"_setting_"<<setting_number;
// 		mkdir(ss1.str().c_str(), S_IWOTH || S_IROTH); 
		mkdir(ss1.str().c_str(), ACCESSPERMS); 
		
		ss1<<"/"<< p_name <<   "_out.txt";
		
		ofstream master_cout ( ss1.str().c_str(),ios_base::app );
		cout << ss1.str().c_str() <<endl;
        
		cout << "time_to_optimize_in_master        : ";
        time_to_optimize_in_master.printTime();
        cout << "time_to_apply_cuts_master         : ";
        time_to_apply_cuts_master.printTime();
        cout << "time_to_calculate_AC at master    : ";
        time_to_calculate_AC.printTime();
        cout << "time_to_create_RW_points          : ";
        time_to_create_RW_points_master.printTime();
		
		master_cout << "time_to_optimize_in_master        : ";
		time_to_optimize_in_master.printTime ( -1,master_cout );
		master_cout << "time_to_apply_cuts_master         : ";
		time_to_apply_cuts_master.printTime ( -1,master_cout );
		master_cout << "time_to_calculate_AC              : ";
		time_to_calculate_AC.printTime ( -1,master_cout );
		master_cout << "time_to_create_RW_points          : ";
		time_to_create_RW_points_master.printTime ( -1,master_cout );
		
        cout << "AVG TIME to optimize in slaves         : "<< avg_time_to_optimize <<endl;
        cout << "AVG TIME to generateCuts in slaves     : "<< avg_time_to_generate_cuts<<endl;
        cout << "AVG TIME to apply cuts   in slaves     : "<< avg_time_to_apply_cuts <<endl;
		cout << "AVG TIME for l1 norm minimization      : "<< avg_time_for_l1_norm << endl;
		cout << "AVG TIME for RW point generation       : "<< avg_time_for_RW_point_generation<< endl;
		cout << "AVG TIME for FP start point generation : "<< avg_time_for_generating_FP_start<< endl;
		
		master_cout <<  "AVG TIME to optimize in slaves         : "<< avg_time_to_optimize <<endl;
		master_cout <<  "AVG TIME to generateCuts in slaves     : "<< avg_time_to_generate_cuts<<endl;
		master_cout <<  "AVG TIME to apply cuts   in slaves     : "<< avg_time_to_apply_cuts <<endl;
		master_cout <<  "AVG TIME for l1 norm minimization      : "<< avg_time_for_l1_norm << endl;
		master_cout <<  "AVG TIME for RW point generation       : "<< avg_time_for_RW_point_generation<< endl;
		master_cout <<  "AVG TIME for FP start point generation : "<< avg_time_for_generating_FP_start<< endl;
		
		best_incumbent_solution.print(master_cout);
		master_cout <<endl;

        stringstream s1;
		s1 << "detailed_results/Exp_" <<  experiment_number<<"_setting_"<<setting_number;
// 		mkdir(s1.str().c_str(), S_IWOTH || S_IROTH); 
// 		mkdir(s1.str().c_str(), ACCESSPERMS); 
		
		s1<<"_summary.txt";
		
	
        ofstream o1 ( s1.str().c_str(),ios_base::app );
        o1 << summarystring.str().c_str() << endl;
        cout << summarystring.str().c_str() << endl;
		
		master_cout << summarystring.str().c_str() << "\t";

		


    }


    return 0;
}

int serial_run ( int argc, char* argv[], string argfilename ="__NULL" ){
	/** NOTE: DO IT AGAIN */

//     double best_incumbent_objective_at_master = starting_objective;

    cpfp_executable my_runner ( argfilename );
    my_runner.set_seed ( my_default_seed );

	if ( my_runner.initialize(solver_type,ipopt_solver) == 1 ) {

		objective_improvement_coefficient_FP = 1;
		objective_improvement_coefficient_CP = 1;
        //cout << " all variables integer or all continious variables has zero objective function coefficient." << endl;

    }

    cout << "objective_improvement_coefficient_FP: " << objective_improvement_coefficient_FP <<endl;

    double LB_for_original_objective;

    My_solution lp_solution;
    Generate_cut_type t = cut_type_generator();

    Timer time_checker;

    time_checker.start();
    bool flagcontinue = true;
    vector<mySerializableRowCut> all_rowcutlist;
    int number_of_cuts = 0;
    vector<mySerializableRowCut> iter_rowcutlist;
    int iteration_counter = 0;

    while ( flagcontinue ) {
        if ( iteration_counter%10 == 0 ) {
            t.set_generate_type ( CUT_TYPE_SIMPLE_ROUNDING );
        } else {
            t.unset_generate_type ( CUT_TYPE_SIMPLE_ROUNDING );
        }
        iteration_counter++;
        time_to_optimize_in_serial.restart();
        int res =my_runner.optimize();
        time_to_optimize_in_serial.pause();
        LB_for_original_objective = my_runner.getObjValue();
        if ( iteration_counter == 1 ) {
            cout << "iteration "<< iteration_counter << "\t\t           starting \tlb: " << LB_for_original_objective<< endl;
        } else {
            cout << "iteration "<< iteration_counter -1<< " END\t\tallrowcuts.size(): "<< number_of_cuts <<  " \tlb: " << LB_for_original_objective<< endl;
        }
        switch ( res ) {
        case 0:
            cout << "LP infeasible " << endl;
            flagcontinue = false;
            break;

        case 1:
            break;
        case -6:
            CPFP_ERROR ( "lp neither optimal nor infeasible ", cout )
            cout << "res is " << res << endl;
            flagcontinue = false;

            break;
        default:
            CPFP_ERROR ( "lp neither optimal nor infeasible ", cout )
            cout << "res is " << res << endl;
            flagcontinue = false;

            break;
        }
        My_solution sol;
        double integer_infeasibility;

        int zero = 0;
        int int_feasibility = my_runner.get_solution ( sol,integer_infeasibility );

        if ( ( res == 1 ) ) {
            if ( ( int_feasibility == 0 ) ) {
                /** FOUND a feasible solution
                * no need to generate cuts
                */
                cout << "LP generates an optimal INTEGER " <<endl;
                flagcontinue = false;
            } else {
                iter_rowcutlist.clear();

                time_to_generate_cuts.restart();
                iter_rowcutlist = my_runner.generate_rowcutlist ( t,all_rowcutlist,zero );
                time_to_generate_cuts.pause();

                iter_rowcutlist =
                    rowCutSelectionProcessAtSlave ( iter_rowcutlist,0 );


                sort ( iter_rowcutlist.begin(),iter_rowcutlist.end() );
                std::vector<mySerializableRowCut>::iterator it;

                it = unique ( iter_rowcutlist.begin(),iter_rowcutlist.end() );
                iter_rowcutlist.resize ( std::distance ( iter_rowcutlist.begin(),it ) );


                number_of_cuts+= ( int ) iter_rowcutlist.size();
                if ( iter_rowcutlist.size() >0 ) {
                    time_to_apply_cuts_serial.restart();
                    my_runner.apply_cuts ( iter_rowcutlist );
                    time_to_apply_cuts_serial.pause();
                } else {

                    flagcontinue = false;

                }
            }
        } else {
            flagcontinue = false;

        }


        if ( iteration_counter > number_of_iterations ) {
            if ( flagcontinue ) {
                time_to_optimize_in_serial.restart();
                my_runner.optimize();
                time_to_apply_cuts_serial.pause();
                LB_for_original_objective = my_runner.getObjValue();
                cout << "iteration "<< iteration_counter << " END\t\tallrowcuts.size(): "<< number_of_cuts <<  " \tlb: " << LB_for_original_objective<< endl;
            }
            flagcontinue = false;
        }

    }
    time_checker.printTime ( time_checker.stop() );

    return 0;

}

int pre_run ( int argc, char *argv[], string &argfilename ){

    time_to_optimize_in_master.reset();
    time_to_optimize_in_slave.reset();
    time_to_optimize_in_serial.reset();
    time_to_generate_cuts.reset();
    time_to_apply_cuts_master.reset();
    time_to_apply_cuts_slave.reset();
    time_to_apply_cuts_serial.reset();
    time_to_calculate_AC.reset();
	time_to_create_RW_points_master.reset();
	time_to_create_RW_points_slave.reset();
    time_checker_for_FP_temp.reset();
    time_checker_for_FP_temp.start();
	time_to_create_starting_points_for_FP_slave.reset();
	l1_norm_timer.reset();
// 	cout << "SLAVE " << SLAVEID << " pid: " << getpid() <<endl;

    char *userName = getenv ( "LOGNAME" );
	
    global_filename += "/home/";
    
	//global_filename+= userName;
	global_filename+="kocu";
	global_filename+="/mps_files/";
    argfilename += "/home/";
    //argfilename+= userName;
	argfilename+= "kocu";
    argfilename+="/mps_files/";

	
	int opt;

	bool FP_file_set = false;
	bool CP_file_set = false;
	bool PWR_file_set = false;
	bool file_Set = false;
	

	
	
    while ( ( opt = getopt ( argc, argv, "Z:W:I:F:C:P:f:O:s:S:9:c:w" ) ) != -1 ) {
        switch ( opt ) {
			case 'Z':
				MY_UNIQ_FILENAME = string(optarg);
				break;
            case 'I':
                MY_SLAVE_ID_from_input = atoi( optarg );
                break;
            case 'W':
                WORLD_SIZE_from_input = atoi(optarg);
                break;
			case 'F':
                FP_param_file.clear();
				FP_param_file << "params/" <<string ( optarg );
				FP_file_set = true;
				break;
			case 'C':
                CP_param_file.clear();
				CP_param_file << "params/" << string ( optarg );
				CP_file_set = true;
				break;
			case 'P':
                PWR_param_file.clear();
				PWR_param_file << "params/" << string ( optarg );
				PWR_file_set = true;
				break;
			case 'f':
				argfilename += string ( optarg );
				global_filename += string ( optarg );
				p_name = string ( optarg );
				file_Set = true;
				break;
        case 'w':
            blockingWait = true;
            break;
		case 'S':
			starting_objective = atof ( optarg );
			break;       
		case '9':
			initial_seed = atoi ( optarg );
			my_seed_generator = Seeder(initial_seed);
			break;
        case 'c':
            communicate_every_x_iterations = atoi ( optarg );
            break;
        case '1':
            stop_fp_when_first_solution_is_found = true;
            break;
        case 'O':
            experiment_number = atoi ( optarg );
            break;
		case 's':
			setting_number = atoi ( optarg );
			break;
        default :

            cout << "Usage: does not include option " << opt << endl;
            cout << "Continue anyway [Y/N] ? ";
            cout.flush();
            char c;
            cin >> c;
            if ( c != 'Y' || c != 'y' ) {
                return 1;
            }
        }
        
    }
    
    signaler = my_signaler(MY_UNIQ_FILENAME, MY_SLAVE_ID_from_input, WORLD_SIZE_from_input);
	signaler.clean_files();
	IFSLAVE signaler.send_signal_as_tag_and_number(0,99,0);
	//abort();
    if(!PWR_file_set){
        PWR_param_file << "params/pp-param";
		//IFMASTER cout << "please set parallel parameters file using -P option "<< endl;
        //IFMASTER cout << "using default parameters " << PWR_param_file.str().c_str() << endl;
        PWR_file_set = true;
	}    
	if(!FP_file_set){
        FP_param_file << "params/fp-param" ; 

		//IFMASTER cout << "please set FP parameters file using -F option "<< endl;
       // IFMASTER cout << "using default parameters " << FP_param_file.str().c_str() << endl;
        FP_file_set = true;
	}
	if(!CP_file_set){
        CP_param_file << "params/cp-param" ;
		//IFMASTER cout << "please set CP parameters file using -C option "<< endl;
        //IFMASTER cout << "using default parameters " << CP_param_file.str().c_str() << endl;
        CP_file_set = true;
	}
	    if ( !file_Set ) {
			//IFMASTER cout << "please set filename using -f option Running for p0033.mps"<< endl;
	        global_filename += "p0033.mps";
	        argfilename += "p0033.mps";
	        p_name = "p0033.mps";
	    }
	
	
	
//     read_FP_parameters_from_file(runner);
    if (read_Parallel_parameters_from_file() !=0) {
		exit(-1);
    }
    
    alpha = initial_alpha;



	if (FP_time_limit * number_of_iterations > max_time_limit) max_time_limit = FP_time_limit * number_of_iterations ;
	if (CP_time_limit * number_of_iterations > max_time_limit) max_time_limit = CP_time_limit * number_of_iterations ;

	if (timelimit > max_time_limit) max_time_limit = timelimit;

	max_time_limit+=20;

// 	max_time_limit = 20;


    
    if ( (strcmp(&p_name.c_str()[strlen (p_name.c_str())-4], ".mps") !=0 ) &&(strcmp(&p_name.c_str()[strlen (p_name.c_str())-3], ".lp") !=0 )){
		IFMASTER cout << " extension is not mps or lp. adding .mps to filename" << endl;
		p_name += ".mps";
		argfilename += ".mps";
		global_filename += ".mps";
	}


    IFMASTER {
        summarystring  << p_name << " " ;
        summarystring2 << p_name << " " ;
    }


//     switch(use_limit){
//         case 1:
//             FP_iterate_decider_slave = bind ( FP_iterate_use_time_limit,placeholders::_1, placeholders::_2);
//             FP_iterate_decider_master = bind(time_limit_master_checker, placeholders::_1);
// 
//             break;
//         case 2:
//             FP_iterate_decider_slave = bind ( FP_iterate_use_min_iter_lim,placeholders::_1 );
//             FP_iterate_decider_master = bind(min_iterlim_master_checker, placeholders::_1);
//             break;
// 
//     };

    /** this part is setting default arguments and fills type dist*/


//     int nc = type_dist_for_CP.get_number_of_instances_that_does_not_need_AC() +type_dist_for_CP.get_number_of_instances_that_need_AC();
//     if ( nc <=1) {
// 		/** only one original objective is running */
// 		int nrem = MPI::COMM_WORLD.Get_size() - nc -1;
// 		if (nrem>0)
//         type_dist_for_CP.add_instance ( OBJECTIVE_TYPE_PERTURBED_FROM_LONG_DIKIN, nrem);
//     }
// 
//     int nf =type_dist_for_FP.get_number_of_instances_that_does_not_need_AC() +type_dist_for_FP.get_number_of_instances_that_need_AC();
// 
//     if ( nf <=1 ) {  
//         /** only one original objective is running */
//         int n_originals_for_FP = 1* ( MPI::COMM_WORLD.Get_size()-1 ) /4;
//         int n_random_for_FP = 3* ( MPI::COMM_WORLD.Get_size()-1 ) /4;
// 
//         if ( n_originals_for_FP < type_dist_for_FP.get_number_of_instances_that_does_not_need_AC() ) {
//             n_originals_for_FP = type_dist_for_FP.get_number_of_instances_that_does_not_need_AC();
//         }
// 
//         int n_rem = MPI::COMM_WORLD.Get_size() - 1 - n_originals_for_FP - n_random_for_FP;
// 
// 		if (n_rem >0){
// 			type_dist_for_FP.add_instance ( OBJECTIVE_TYPE_ORIGINAL, n_originals_for_FP );
// 			type_dist_for_FP.add_instance ( OBJECTIVE_TYPE_RANDOM_FROM_LONG_DIKIN, n_random_for_FP );
// 
// 			if ( 0 ) {
// 				type_dist_for_FP.add_instance ( OBJECTIVE_TYPE_PERTURBED_FROM_LONG_DIKIN, n_rem );
// 			}
// 			else{
// 				type_dist_for_FP.add_instance ( OBJECTIVE_TYPE_RANDOM_FROM_LONG_DIKIN, n_rem);
// 
// 			}
// 		}
//     }
//     else if (nf < MPI::COMM_WORLD.Get_size() - 1){
//         int n_rem = MPI::COMM_WORLD.Get_size() - 1 - nf;
//         type_dist_for_FP.add_instance ( OBJECTIVE_TYPE_RANDOM_FROM_LONG_DIKIN, n_rem);
//     }

    return 0;
}

#ifdef USING_IOPTIMIZE
int tester ( string argfilename ){
    ioptimize_wrapper test_iop;
    string sl= "slave_log.txt";
    string ml= "master_log.txt";

    int n_random_points = 10;

    test_iop.set_seed ( my_default_seed );
    test_iop.setCenProbType ( solve_original_centering_problem );
    test_iop.initialize2 ( argfilename );
    IFSLAVE test_iop.redirectoutput ( sl );
    IFMASTER test_iop.redirectoutput ( ml );
// 	test_iop.set_seed ( my_default_seed );
    test_iop.use_paramOpt();





    if ( type > 3 ) {
        test_iop.set_walk_type ( type-3 );
    } else {
        test_iop.set_walk_type ( type );
    }

    req_stat_vectors r;

    IFMASTER {

//         test_iop.calculate_analytic_center_self();
        for ( int i=1; i < WORLDSIZE; ++i ) {
//             test_iop._send_mpi ( i,TAG_SENDING_ANALYTIC_CENTER, r.v_req, r.v_stat );
        }
// 		test_iop.get_all_IPA();
// 		test_iop.get_all_DPA();

    }

    IFSLAVE {
//         test_iop._recv_mpi ( 0,TAG_SENDING_ANALYTIC_CENTER, r.v_req, r.v_stat );
    }

// 	test_<ioptimize_wrapper>(test_iop, false,false);


    vector<vector<double> >master_random_points;
    vector<vector<double> >slave_random_points;


    IFMASTER {


// 		COUT_TESTER_D(__LINE__, " M")
        test_iop.iop_set_seed ( 2 );
// 		COUT_TESTER_D(__LINE__, " M")
// 		test_iop.print(cout, 2);


        test_iop.init_sampler();
// 		COUT_TESTER_D(__LINE__, " M")
// 		ofstream out("ml.txt");

// 		test_iop.print(out, 1);
        test_iop.get_next_n_random_point ( master_random_points,n_random_points );

// 		test_iop.get_task_pointer();
// 		IOP_CHECK_CALL ( Iop_get( test_iop.get_task_pointer(), IOP_COS_PROB_CEN ) );

// 		for (int i = 0; i < n_random_points;++i){
// 			vector<double> rp;
// 			test_iop.get_next_random_point(rp);
// 			master_random_points.push_back(rp);
// 		}
// 		COUT_TESTER_D(__LINE__, " M")
// 		test_iop.writeParams("mparams.txt");

    }


    IFSLAVE {
        // 		test_iop.set_seed(2);


        test_iop.after_receive_update();
        test_iop.iop_set_seed ( 2 );

        // 		usleep(10000000);

        // 		cout << "init sampler: " <<
        test_iop.init_sampler();
        // 		COUT_TESTER_D(__LINE__, " S")
// 		ofstream out("sl.txt");
        // 		test_iop.print(out, 1);
        // 		test_iop.write_Mps("after_init_sampler");
        // 		cout << " get_next_n_random_point " << endl;
        // 		cout << "--------" <<endl;
// 		for (int i = 0; i < n_random_points;++i){
// 			vector<double> rp;
// 			test_iop.get_next_random_point(rp);
// 			slave_random_points.push_back(rp);
// 		}
        // 		cout <<
        test_iop.get_next_n_random_point ( slave_random_points,n_random_points );
        // 		cout << endl << "--------" <<endl;
        //

        // 		slave_random_points.resize(1);
        // 		test_iop.get_next_random_point(slave_random_points[0]);
        // 		test_iop.write_Mps("after_random_point");

        // 		COUT_TESTER_D(__LINE__, " S")
// 		test_iop.writeParams("sparams.txt");
    }

// 	COUT_TESTER(__LINE__)
// 	usleep(1000000);
// 	COUT_TESTER(__LINE__)

    for ( int i = 0; i < n_random_points; ++i ) {

        IFMASTER {
            cout << "MMM " << i+1<<": ";
// 			vector_line_print<double>(master_random_points[i],cout);
        }
        IFSLAVE {
            cout << "S" << SLAVEID +1 << " " << i+1<<": ";
// 			vector_line_print<double>(slave_random_points[i],cout);
        }
// 		usleep(500000);
        IFMASTER cout << endl;
// 		usleep(500000);

// // 		int k;
// // 		cin >> k;
    }
    char *hostName = getenv ( "HOSTNAME" );

    cout << " process " << SLAVEID+1 << " hostname " << hostName << endl;

    return 0;
}
#endif
/** current version works only for original FP NO RANDOM WALK*/
int serial_FP_run_tester( int argc, char *argv[], string argfilename){
	time_for_FP.start();

	running_serial = true;
	cpfp_executable my_runner ( argfilename );
	unsigned  my_seed = my_seed_generator.lth_seed_lcm ( MACRO_SLAVE_ID + 1 );
	
	
	my_runner.set_seed ( my_seed );
	if ( my_runner.initialize(solver_type,ipopt_solver) == 1 ) {
		objective_improvement_coefficient_FP = 1;
		objective_improvement_coefficient_CP = 1;
        cout << " no problem in initialization" << endl; 
	}
	vector<int> column_types = my_runner.get_column_types();
	for (unsigned i =0; i < column_types.size();++i){
		if (column_types[i]>0) integer_columns_in_slave.push_back(i);
	}
	double sb = my_runner.calculate_pseudo_bound();
	my_runner.set_pseudo_bounds(sb);
    
    int res1 = my_runner.optimize(1);
	//cout << "  with first lp res1: "<<res1 <<endl;
	int res = my_runner.optimize(2);
	//cout << "  with first lp res: "<<res <<endl;

	double LB_for_original_objective = my_runner.getObjValue(1);
    cout << "LB_for_original_objective: "<<LB_for_original_objective <<endl;

	my_runner.update_objective_cutoff_constraint_lb ( LB_for_original_objective );
	int current_stage = 0;
	int FP_total_iteration_counter = 0;

// 	LP_WRITE(my_runner,current_stage,FP_total_iteration_counter, __LINE__)


	My_solution current_lp_optimum;
	current_lp_optimum.set_original_size(my_runner.get_ncols());
// 	double d_frac;
//  	int i_frac = my_runner.get_solution(current_lp_optimum, d_frac);

	My_solution starting_solution;
	My_solution rounded_solution;
	My_solution previous_rounded_solution;
	My_solution best_incumbent_solution_at_slave_FP;
	My_solution best_starting_solution_so_far;



	int stage1_iterlim;
	int stage2_iterlim;
	int stage1_resetlim;
	int stage2_resetlim;

    read_FP_parameters_from_file(my_runner);
	my_runner.update_objective_cutoff_constraint_ub(starting_objective);

	my_runner.get_FP_INT_PARAM(ENM_FP_INT_PARAM_MAXITER1, stage1_iterlim);
	my_runner.get_FP_INT_PARAM(ENM_FP_INT_PARAM_MAXITER2, stage2_iterlim);
	my_runner.get_FP_INT_PARAM(ENM_FP_INT_PARAM_MAXITERW1, stage1_resetlim);
	my_runner.get_FP_INT_PARAM(ENM_FP_INT_PARAM_MAXITERW2, stage2_resetlim);

	int rounding_type = initial_rounding_type;
	//	runner.get_FP_INT_PARAM(ENM_FP_INT_PARAM_ROUNDING, rounding_type);
	// 	runner.get_FP_INT_PARAM(ENM_FP_INT_PARAM_MINFLIP,minflip);

	my_runner.get_FP_DOUBLE_PARAM(ENM_FP_DOUBLE_PARAM_ALPHA_INIT, initial_alpha);
// 	int tere =
	my_runner.get_FP_DOUBLE_PARAM(ENM_FP_DOUBLE_PARAM_ALPHA_QUOD, alpha_reduction);
	my_runner.get_FP_DOUBLE_PARAM(ENM_FP_DOUBLE_PARAM_ALPHA_DIST, alpha_dist);

    my_runner.print_FP_parameters();
    
    //DEBWAIT("stage1_iterlim: ", stage1_iterlim)

	int ncols = my_runner.get_ncols();

	int nBinInt = my_runner.get_number_of_binaries();
	int nGenInt = my_runner.get_number_of_general_integers();
	int nCont   = ncols - nBinInt - nGenInt;
	my_runner.init_FP_stage2();

	cout << "seed: " << my_runner.get_seed() << endl;


	int missedDecr = 0;
	int maxMissedDecr = stage1_resetlim;
	int minFracIt = 0;
	double minFracDbl = std::numeric_limits< double >::max();

// 	int minFracInt = ncols;

	double best_incumbent_objective_at_serial_FP = std::numeric_limits< double >::max();

	starting_solution.set_original_size(ncols);
	rounded_solution.set_original_size(ncols);
	previous_rounded_solution.set_original_size(ncols);
	best_incumbent_solution_at_slave_FP.set_original_size(ncols);
	best_starting_solution_so_far.set_original_size(ncols);


	int big_iteration_counter = 0;
	int external_iteration_counter =0;

	bool fp_continue_flag = true;

	int stage_1_restarts = 0;
	int stage_2_restarts = 0;
	int total_s1_restarts = 0;
	int total_s2_restarts = 0;
	bool original_continues = false;



	bool found = false;



	while (fp_continue_flag){
		external_iteration_counter++;
		big_iteration_counter++;
		double curFrac_d = std::numeric_limits< double >::max();
		int curFrac_i = ncols;

		if (external_iteration_counter == 1 && !original_continues){
			/** create initial starting solution*/
			static int ki = 0;
			ki++;
			//cout << "\t\t\t\tXXXXXXXXXXXXXXXXXXXX in line : " << __LINE__ << " for the " << ki <<"th time " << endl;

			int opt_res = my_runner.optimize();
			curFrac_i = my_runner.get_solution(starting_solution,curFrac_d);
            //cout << "curFrac_i:" << curFrac_i << endl;
			//best_starting_solution_so_far.clear();
// 			cout << "--" << endl;
// 			LP_WRITE(my_runner,current_stage,FP_total_iteration_counter, __LINE__)

// 			starting_solution.shortlineprint_selected_indices(integer_columns_in_slave,false, cout);

// 			cout << "--XX" << endl;
			if (opt_res<=0) {
				cout << "problem in solving lp line: "<< __LINE__ << " opt_res: "<< opt_res << " external_iteration_counter: " << external_iteration_counter << " big_iteration_counter: "<< big_iteration_counter << endl;
				break;
			}
			stage_1_restarts =0;
			stage_2_restarts =0;
			minFracDbl = curFrac_d;
// 			minFracInt = curFrac_i;
			minFracIt = FP_total_iteration_counter;
			if(curFrac_i ==0){ /** int sol*/
				double new_obj;
				if (my_runner.calculate_original_objective_value_of_solution(new_obj,starting_solution) != 0) {
					cout << "problem in line: " << __LINE__ << endl;
				}
				cout << "Serial FP finds a solution in line "<<__LINE__ << " with obj: " << starting_solution.get_original_objective() << endl;
				if (current_stage ==2){
					my_runner.unset_bounds_FP_stage_2();
				}

				if (nCont >0){
					My_solution polished_sol;
					double polished_obj;
					int pol_res  = my_runner.polish_integer_solution(starting_solution,polished_sol, polished_obj);
					if (pol_res ==1){
						starting_solution = polished_sol;
					}
				}
				best_starting_solution_so_far = starting_solution;
				cout << "Serial FP  finds a solution with obj: " << starting_solution.get_original_objective() << endl;
				if(best_incumbent_objective_at_serial_FP> starting_solution.get_original_objective()){
					best_incumbent_objective_at_serial_FP= starting_solution.get_original_objective();
					best_incumbent_solution_at_slave_FP = starting_solution;
					found = true;
				}
			}
			else{
				previous_rounded_solution =  rounded_solution;
				rounding_type = rounding_decider(stage_1_restarts,current_stage);
// 				cout << " line " << __LINE__ <<endl;

// 				cout << "--" << endl;
 				starting_solution.shortlineprint_selected_indices(integer_columns_in_slave, false, cout);


// 				cout << "--" << endl;
				my_runner.round_sol(starting_solution,rounded_solution,rounding_type,current_stage);
// 				cout << " line " << __LINE__ <<endl;

			}
			DEB_WAIT(__LINE__, ++_ux,FP_total_iteration_counter )

		}
		else{
			static int kix = 0;
			kix++;
			//cout << "\t\t\t\tYYYYYYYYYYYYYYYYYYYYYYY in line : " << __LINE__ << " for the " << kix <<"th time " << endl;
			history_clear_version1();
			previous_rounded_solution = rounded_solution;
// 			rounded_solution.clear();
			DEB_WAIT(__LINE__, ++_ux,FP_total_iteration_counter )
            alpha = initial_alpha;
		
			
		}

		/** starting solution and rounded solution are generated*/
		/** STAGE 1*/
		// 		cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << endl;
		DEBOFF (found, nBinInt,stage1_iterlim, fp_continue_flag)
		if(!found && (nBinInt>0) && stage1_iterlim >0 && fp_continue_flag){
			missedDecr =0;
			minFracDbl = std::numeric_limits< double >::max();
// 			minFracInt = ncols;
			minFracIt = FP_total_iteration_counter;
            DEBOFF(11)
            
			int local_iteration = 0;
			int local_iterlim = stage1_iterlim;
			maxMissedDecr = stage1_resetlim;
			current_stage = 1;
			stage_1_restarts =0;

// 			rounding_type = initial_rounding_type;
			// 			rounded_solution = best_starting_solution_so_far;

			curFrac_d = std::numeric_limits< double >::max();
			curFrac_i = ncols;
			while ((local_iteration < local_iterlim) && (fp_continue_flag)){
// 				cout << "iteration "<< FP_total_iteration_counter << " " <<  curFrac_i << endl;

				FP_total_iteration_counter++;
				local_iteration++;
				REDUCE_ALPHA
				My_objective_function ob = my_runner.f_type_stage_1_create_obj_for_FP(rounded_solution,alpha);
				my_runner.set_FP_objective_function(ob);
				my_runner.unset_bounds_FP_stage_2();
				stringstream z;
				z<< "stage_"<<current_stage<< "_iter_" << FP_total_iteration_counter;
// 				cout << "FILE_WRITE exporting " << z.str() <<endl;

// 				my_runner.puke_lp_file(z.str().c_str());

				DEB_WAIT(__LINE__, ++_ux,FP_total_iteration_counter )

				int res2 = my_runner.optimize(1);
				if (res2<=0){
					/** to be done*/
					/** LP not optimal */
					cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << " LP not optimal " << endl;

					break;
				}
				My_solution pre_start = starting_solution;
				curFrac_i = my_runner.get_solution(starting_solution, curFrac_d);
				DEB_WAIT(__LINE__, ++_ux,FP_total_iteration_counter )

// 				if (pre_start == starting_solution) {
//
// 					cout << " starting solution did not change line: " <<__LINE__ << " iter: " << FP_total_iteration_counter <<endl;
//
// 				}
				curFrac_i = my_runner.get_binary_fractionality_of_solution(starting_solution,curFrac_d);
// 				if (pre_start == starting_solution) {
// 					cout << " starting solution did not change line: " <<__LINE__  <<  " iter: " << FP_total_iteration_counter << " curFrac_i: "<<curFrac_i  <<endl;
//
// 				}

				previous_rounded_solution = rounded_solution;
				if(curFrac_d < minFracDbl){
					if (curFrac_d/minFracDbl < 0.9)
						missedDecr = 0;
					minFracDbl = curFrac_d;
// 					minFracInt = curFrac_i;
					minFracIt = FP_total_iteration_counter;

					best_starting_solution_so_far = starting_solution;
				}
				else
					missedDecr++;

				rounding_type = rounding_decider(stage_1_restarts,current_stage);
				int n_changed = my_runner.round_sol(starting_solution, rounded_solution, rounding_type ,current_stage);

// 				cout << "line " << __LINE__ << " "<< rounding_type <<endl;

				DEB_WAIT(__LINE__, ++_ux,FP_total_iteration_counter )


				if(curFrac_i == 0) {
// 					cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << " binary variables did not change at iteration: " << FP_total_iteration_counter << endl;
					break; /** binary variables did not change*/
				}
				if (missedDecr > maxMissedDecr && stage2_iterlim>0) {
// 					cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << " too many itarations without 10% fractionality improvement stage1 " <<  endl;
					break; /** too many itarations withour 10% fractionality improvement */
				}
				if (stage_1_restarts > 100) {
					// 					total_s1_restarts+= stage_1_restarts;
					// 					stage_1_restarts = 0;
// 					cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << " too many stage_1_restarts" <<endl;
					break;
				}

// 				cout << " iteration " << FP_total_iteration_counter << " n_changed: " << n_changed <<endl;

				if (n_changed ==0){

					n_changed += my_runner.getNextIntegerPoint(starting_solution,rounded_solution,current_stage);

				}
				DEB_WAIT(__LINE__, ++_ux,FP_total_iteration_counter )

				My_solution pre;
				// 				pre.set_original_size(ncols);
				bool insert = history_insert_version1(rounded_solution,alpha,pre) ;
				while(!insert && stage_1_restarts < 10*local_iterlim){
					stage_1_restarts++;

					int changed= my_runner.f_type_restart_1_or_2(starting_solution,rounded_solution,pre,current_stage,FP_total_iteration_counter);
					DEB_WAIT(__LINE__, ++_ux,FP_total_iteration_counter )

					if (changed ==0){
						continue;
					}
					else insert =   history_insert_version1(rounded_solution,alpha,pre) ;

				}

				fp_continue_flag = FP_iterate_decider_slave( FP_total_iteration_counter,time_for_FP.stop());
// 				cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ << " fp_continue_flag " << fp_continue_flag<<  endl;


			}
        }
        int res_ext = my_runner.optimize_extra_mip(starting_solution, current_stage);
		if (res_ext == 1){
			My_solution extra_sol = my_runner.get_solution_of_extra_mip(current_stage);
			starting_solution = extra_sol;

		}
		curFrac_i = my_runner.get_fractionality_of_solution(starting_solution,curFrac_d);
		if (curFrac_i == 0) {
			double new_obj;
			if (my_runner.calculate_original_objective_value_of_solution(new_obj,starting_solution) != 0) {
				cout << "problem in line: " << __LINE__ << endl;
			}
			if (current_stage ==2){
				my_runner.unset_bounds_FP_stage_2();
			}


			if (nCont > 0){
				My_solution polished_sol;
				double polished_obj;
				int pol_res  = my_runner.polish_integer_solution(starting_solution,polished_sol, polished_obj);
				if (pol_res ==1){
					starting_solution = polished_sol;
				}
			}
			best_starting_solution_so_far = starting_solution;

			/** found a feasible solution*/
			if(best_incumbent_objective_at_serial_FP > starting_solution.get_original_objective()){


				best_incumbent_objective_at_serial_FP = starting_solution.get_original_objective();
				best_incumbent_solution_at_slave_FP = starting_solution;
				found = true;

			}


		}
		history_clear_version1();

		/** STAGE 2*/
		if (!found && stage2_iterlim>0 && fp_continue_flag){


			missedDecr =0;
			minFracDbl = std::numeric_limits< double >::max();
// 			minFracInt = ncols;
			minFracIt = FP_total_iteration_counter;

			int local_iteration = 0;
			int local_iterlim = stage2_iterlim;
			maxMissedDecr = stage2_resetlim;
			current_stage = 2;
			stage_2_restarts =0;
// 			rounding_type = initial_rounding_type;
			// 			runner.get_FP_INT_PARAM(ENM_FP_INT_PARAM_ROUNDING,rounding_type);
			curFrac_d = std::numeric_limits< double >::max();
			curFrac_i = ncols;
			starting_solution = best_starting_solution_so_far;
			cout << "Using best point from iter: " << minFracIt << endl;

			DEB_WAIT(__LINE__, ++_ux,FP_total_iteration_counter )

			while ((local_iteration < local_iterlim) && (fp_continue_flag)){
				FP_total_iteration_counter++;
				local_iteration++;
				REDUCE_ALPHA
				My_objective_function ob = my_runner.f_type_stage_2_create_obj_for_FP(rounded_solution,alpha);
				my_runner.set_FP_objective_function(ob);
				my_runner.unset_bounds_FP_stage_2();
// 				my_runner.puke_lp_file("temp");

				my_runner.f_type_stage_2_set_bounds(rounded_solution);
// 				stringstream z;
// 				z<< "stage_"<<current_stage<< "_iter_" << FP_total_iteration_counter;
// 				cout << "FILE_WRITE exporting " << z.str() << ".lp  alpha:"<< alpha << endl;

// 				my_runner.puke_lp_file(z.str().c_str());
// 				ob.print();
				DEB_WAIT(__LINE__, ++_ux,FP_total_iteration_counter )

				int res2 = my_runner.optimize();
				if (res2<=0){
					/** to be done*/
					/** LP not optimal */
					break;
				}

				curFrac_i = my_runner.get_solution(starting_solution, curFrac_d);
// 				starting_solution.print();
				DEB_WAIT(__LINE__, ++_ux,FP_total_iteration_counter )
				previous_rounded_solution= rounded_solution;
				rounding_type = rounding_decider(stage_2_restarts,current_stage);
				int n_changed = my_runner.round_sol(starting_solution, rounded_solution,rounding_type, current_stage);
				curFrac_i = my_runner.get_fractionality_of_solution(starting_solution,curFrac_d); /** might not be needed */
				// external_frac_tester(starting_solution,temp,__LINE__);

				if(curFrac_d < minFracDbl){
					if (curFrac_d/minFracDbl < 0.9)
						missedDecr = 0;
					minFracDbl = curFrac_d;
// 					minFracInt = curFrac_i;
					minFracIt = FP_total_iteration_counter;

					best_starting_solution_so_far = starting_solution;
				}
				else
					missedDecr++;

				DEB_WAIT(__LINE__, ++_ux,FP_total_iteration_counter )

				if(curFrac_i == 0) { /** we have a feasible solution */
// 					cout << "line: " << __LINE__ << " curFrac_i: " << curFrac_i << endl;
					double new_obj;
					if (my_runner.calculate_original_objective_value_of_solution(new_obj,starting_solution) != 0) {
						cout << "problem in line: " << __LINE__ << endl;
					}
					if (current_stage ==2){
						my_runner.unset_bounds_FP_stage_2();
					}
// 					cout << "SLAVE " << SLAVEID << " finds a solution in line: "<< __LINE__ << " with obj: " << starting_solution.get_original_objective() << endl;

					if (nCont > 0){
						My_solution polished_sol;
						double polished_obj;
// 						cout << "polishing: "<< endl;

						int pol_res  = my_runner.polish_integer_solution(starting_solution,polished_sol, polished_obj);
						// external_frac_tester(starting_solution,temp,__LINE__);

// 						cout << "polished objective: "<< polished_obj << " pol_res: " << pol_res << endl;
						//
						if (pol_res ==1){
							starting_solution = polished_sol;
							// //                     usleep(100000000);
						}
					}
					best_starting_solution_so_far = starting_solution;
// 					cout << "SLAVE " << SLAVEID << " finds a solution in line: "<< __LINE__ << " with obj: " << starting_solution.get_original_objective() << endl;

					/** found a feasible solution*/
					if(best_incumbent_objective_at_serial_FP> starting_solution.get_original_objective()){

						best_incumbent_objective_at_serial_FP = starting_solution.get_original_objective();
						best_incumbent_solution_at_slave_FP = starting_solution;
						found = true;
					}
					// external_frac_tester(starting_solution,temp,__LINE__);

					break;
				}
				if (missedDecr > maxMissedDecr){
					break;
				}
				if (stage_2_restarts>100) {
					// 					cout << "SLAVE: " << SLAVEID << " line: " << __LINE__ <<  " too many stage_2_restarts " <<endl;
					break; /** too many restarts     */
				}


// 				cout << "n_changed: " <<n_changed <<endl;
				if (n_changed ==0){
					n_changed += my_runner.getNextIntegerPoint(starting_solution,rounded_solution,current_stage);

				}
				DEB_WAIT(__LINE__, ++_ux,FP_total_iteration_counter )

				My_solution pre;
				bool insert = history_insert_version1(rounded_solution,alpha,pre) ;
				while(!insert && stage_2_restarts < 100){
					stage_2_restarts++;
// 					int change =
					my_runner.f_type_restart_1_or_2(starting_solution,rounded_solution,pre,current_stage,FP_total_iteration_counter);

					DEB_WAIT(__LINE__, ++_ux,FP_total_iteration_counter )
					insert =   history_insert_version1(rounded_solution,alpha,pre) ;


				}

				fp_continue_flag = FP_iterate_decider_slave(FP_total_iteration_counter,time_for_FP.stop());


			}

			curFrac_i = my_runner.get_fractionality_of_solution(starting_solution,curFrac_d);
			if (curFrac_i == 0) {
				double new_obj;
				if (my_runner.calculate_original_objective_value_of_solution(new_obj,starting_solution) != 0) {
					cout << "problem in line: " << __LINE__ << endl;
				}

				if(best_incumbent_objective_at_serial_FP > starting_solution.get_original_objective()){
					best_incumbent_objective_at_serial_FP = starting_solution.get_original_objective();
					best_incumbent_solution_at_slave_FP = starting_solution;
					found = true;
				}

			}
			history_clear_version1();
		}
		if(found){
			
			UPDATE_RUNNER_FP(best_incumbent_objective_at_serial_FP, my_runner)

			found = false;
		}

		total_s1_restarts += stage_1_restarts;
		total_s2_restarts += stage_2_restarts;


		if (fp_continue_flag) fp_continue_flag = FP_iterate_decider_slave(FP_total_iteration_counter,time_for_FP.stop());
        
        
		cout << "SLAVE " << SLAVEID<< " time: " << time_for_FP.stop() << " it: " << FP_total_iteration_counter << " fp_C: "<< fp_continue_flag << " best_obj: " << best_incumbent_objective_at_serial_FP<<  endl;

	}

	cout << " serial FP finalizes at line "<< __LINE__ << " with obj: " << best_incumbent_objective_at_serial_FP <<endl;


	FP_iterate_decider_slave(-1,-1);
	my_runner.undo_FP_stage2();



	return 0;
}

bool using_AC3=true;

int main ( int argc, char *argv[] ){
	if (using_AC3) ipopt_solver = new OsiIpoptSolverInterface;
	global_break_timer.start();
    global_time.reset();

	CP_FP_decider_int = CP_FP_n_rounds_of_FP_than_m_CP;
	
    FP_iterate_decider_master = bind(time_limit_master_checker, placeholders::_1);
    FP_iterate_decider_slave = bind (FP_iterate_use_time_limit,placeholders::_1, placeholders::_2);
    #ifdef USING_MPI
    MPI::Init_thread ( argc, argv, MPI_THREAD_MULTIPLE );
    #endif

    IFMASTER global_time.start();
    string filename;
	DEBON();
    pre_run ( argc, argv, filename );
	DEBON("pre_run completed");
	if (!file_exists(global_filename)){
		DEBON("file does not exist");
		IFMASTER 
			cout << " file " << global_filename.c_str() << " does not exists \n breaking. " << endl;
		#ifdef USING_MPI
			MPI_Barrier(MPI::COMM_WORLD);
			MPI_Finalize();
		#endif
		return -1;
	}

    IFMASTER {
        type_dist_for_CP.print();
        type_dist_for_FP.print();
		cout << " max_time_limit: " <<max_time_limit <<endl;
    }
    DEBON();
    //DEBON()
    if ( signaler.get_MY_WORLD_SIZE()==1 ) {
        cout << "running serial" << endl;
        serial_FP_run_tester ( argc,argv, filename );
    } else {
		DEBON();
        IFMASTER master_run ( argc, argv,filename );
        IFSLAVE  slave_run ( argc, argv,filename );
    }
    #ifdef USING_MPI
	MPI_Barrier(MPI::COMM_WORLD);
	#endif
    post_run();

    #ifdef USING_MPI
    MPI_Barrier(MPI::COMM_WORLD);
    MPI_Finalize();
    #endif
    return 0;
}


