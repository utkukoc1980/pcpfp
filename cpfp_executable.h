
#include "myserializablerowcut.h"
#include "my_solution.h"
#include "my_objective_function.h"
#include "parallel_cpfp_stuff.h"
#include <queue>
#include <random>

// #define USE_SET

// #define FP_DEFAULT_PARAMETER_SETTINGS \
//     FP_INT_PARAM[ENM_FP_INT_PARAM_MAXITER1]= 10000;\
//     FP_INT_PARAM[ENM_FP_INT_PARAM_MAXITER2] = 2000;\
//     FP_INT_PARAM[ENM_FP_INT_PARAM_MAXITERW1] = 70;\
//     FP_INT_PARAM[ENM_FP_INT_PARAM_MAXITERW2] = 600;\
//     FP_INT_PARAM[ENM_FP_INT_PARAM_MINFLIP]= 20 ;\
//     FP_INT_PARAM[ENM_FP_INT_PARAM_ROUNDING]= 8;\
//     FP_INT_PARAM[ENM_FP_INT_PARAM_OPTIMIZATION_ALGORITHM]  = 1;\
//     FP_INT_PARAM[ENM_FP_INT_PARAM_LP_PRESOLVE] = 0;\
//     FP_INT_PARAM[ENM_FP_INT_PARAM_MIP_PRESOLVE] = 0;\
//     FP_DBL_PARAM[ENM_FP_DOUBLE_PARAM_FLIP_TRESHOLD]=  0.001; \
//     FP_DBL_PARAM[ENM_FP_DOUBLE_PARAM_ALPHA_INIT] = 1.0; \
//     FP_DBL_PARAM[ENM_FP_DOUBLE_PARAM_ALPHA_QUOD] = 0.9; \
//     FP_DBL_PARAM[ENM_FP_DOUBLE_PARAM_ALPHA_DIST] = 0.005;


#define epInt 1e-5
    
    
#define OPT_RESULT(xmip)\
    int retval =-9; \
    if ( xmip->isProvenOptimal() ) retval= 1;\
    if ( xmip->isProvenPrimalInfeasible() ) retval = 0;\
    if ( xmip->isDualObjectiveLimitReached() ) retval = -2;\
    if ( xmip->isIterationLimitReached() )     retval = -3;\
    if ( xmip->isPrimalObjectiveLimitReached() ) retval = -4;\
    if ( xmip->isProvenDualInfeasible() ) retval = -5;\
    if ( xmip->isAbandoned() ) retval = -6;\
    if (retval<=0) {\
        stringstream s;\
        s << "prob_" << " lp_" << retval;\
        xmip->writeLp(s.str().c_str());\
    }\
    return retval;
    
    
    
/** parameter set */
enum FP_INT_PARAM_SET {
    ENM_FP_INT_PARAM_BGN                =    0, /**< begin tag of the parameter setting */
    ENM_FP_INT_PARAM_MAXITER1            =    0,				/** maxit1	   = 0 to disable stage1 or max iterations for stage 1 (10000) */
    ENM_FP_INT_PARAM_MAXITER2            =    1, /**	 maxit2	   = 0 to disable stage2 or max iterations for stage 2 (2000)\n*/
    ENM_FP_INT_PARAM_MAXITERW1            =   2,  /** 	 maxit1wi  = max stage1 iterations without 10% improvement (70)\n" */
    ENM_FP_INT_PARAM_MAXITERW2            =    3, /**   maxit2wi  = max stage2 iterations without 10% improvement (600)\n" */
    ENM_FP_INT_PARAM_MINFLIP             = 4,       /** minimum # of variables to flip per iteration (20)*/
    ENM_FP_INT_PARAM_ROUNDING             = 5,      /** s|r|h|q|a : rounding threshold (standard 1/2, random in 0-1,\n"
				"					 random in .25-.75, random with quadratic distribution, auto) (a)\n" */
    ENM_FP_INT_PARAM_OPTIMIZATION_ALGORITHM = 6,  /** p|d|b|a (primal, dual, barrier, auto*/
    ENM_FP_INT_PARAM_LP_PRESOLVE = 7,               /** use presolving on LPs */
    ENM_FP_INT_PARAM_MIP_PRESOLVE = 8,              /** use presolving on mip prior to FP **/
    ENM_FP_INT_PARAM_END                   =    9, /**< end tag of the parameter setting */
};


enum FP_DOUBLE_PARAM_SET {
    ENM_FP_DOUBLE_PARAM_BGN =0,
    ENM_FP_DOUBLE_PARAM_FLIP_TRESHOLD = 0, /**integrality threshold below which variables are not flipped =0.001 */
    ENM_FP_DOUBLE_PARAM_ALPHA_INIT = 1, /** alphainit = start value of the original objective function fraction (1.0)\n */
    ENM_FP_DOUBLE_PARAM_ALPHA_QUOD = 2, /** alphaquot = redunction factor of alpha (0.9)*/
    ENM_FP_DOUBLE_PARAM_ALPHA_DIST = 3, /** alphadist = difference up to which two alpha values are considered equal (0.005)*/
    ENM_FP_DOUBLE_PARAM_END =4,


};




#ifndef GENERATE_CUT_TYPES_H
#define GENERATE_CUT_TYPES_H
class Generate_cut_type
{
    int generate_cut[CUT_TYPE_END]; /** 0 for do not generate that cut
										 = 1 generate that type of cut
										 > 1 may be used later */

public:
    Generate_cut_type();
    Generate_cut_type ( const Generate_cut_type &other );
    Generate_cut_type &operator= ( const Generate_cut_type &other );
    void generate_all();
    void generate_none();
    void set_generate_type ( CUT_TYPE_ENM index = CUT_TYPE_BGN );
    void unset_generate_type ( CUT_TYPE_ENM index = CUT_TYPE_BGN );
    int get_generate_type ( CUT_TYPE_ENM index = CUT_TYPE_BGN );
#ifdef USING_MPI
    void _send_mpi ( int destination, int tag, vector<MPI_Request> &v_request, vector<MPI_Status> v_status, MPI_SEND_TYPE SEND_TYPE = MPI_ISEND_NONBLOCKING );
    void _recv_mpi ( int source, int tag, vector<MPI_Request> &v_request, vector<MPI_Status> v_status, MPI_RECV_TYPE RECV_TYPE = MPI_IRECV_NONBLOCKING );
#else
        void _send_alt1(int destination, int tag, my_signaler communicator, ALT1_SEND_TYPE SEND_TYPE = ALT1_SEND_NONBLOCKING);
        void _recv_alt1(int source,int tag, my_signaler communicator,  ALT1_RECV_TYPE RECV_TYPE = ALT1_RECV_NONBLOCKING);
#endif
    void print ( std::ostream & out = std::cout );

};

#endif

#ifndef FP_INFORMATION_CLASS_H
#define FP_INFORMATION_CLASS_H

class FP_additional_information
{
public:

    int n_rows_added;
    int n_cols_added;

    vector<int> added_row_indices1;
    vector<int> added_row_indices2;

    vector<int> added_col_indices;

    vector<int> binary_indices;
    vector<int> general_integer_indices;


//  bool org_objective_maximize;
//  void change_rhs(OsiSolverInterface *mip, vector<double> new_solution);

};

class FP_information{
public:
    int iterlim;
    int itercount;
    int resetlim;
    int resetcount;
    int rounding_type;


    int stage2_iterlim;
    int stage2_itercount;
    int stage2_resetlim;
    int stage2_resetcount;
};

#endif

#ifndef FP_HISTORY
#define FP_HISTORY

class FP_history{
public:

    vector<int> int_frac_list;
    vector<double> double_frac_list;
    vector<My_solution> solution_history;
    vector<double> alpha_list;

    bool compare_ids ( unsigned i1, unsigned i2 );


//  bool org_objective_maximize;
//  void change_rhs(OsiSolverInterface *mip, vector<double> new_solution);
    int compare_last_with_previous ( bool use_hist = false );
    int compare_last_with_previous_v1 ( bool use_hist = false );
    void clear ( unsigned keep_latest_n = 0 );


};
#endif

#ifndef CPFP_EXECUTABLE_H
#define CPFP_EXECUTABLE_H

class cpfp_executable{

    /** these are defined for line search*/
    vector<int> LS_indices;
    vector<double> LS_lambda;
    vector<double> LS_one_over_dif;
    vector<int> LS_d_k; 
    
    
    OsiClpSolverInterface *mipCP;
    OsiSolverInterface *mip;
    OsiSolverInterface *extra_mip;

    
    OsiSolverInterface *AC_mip;
    bool initialized_AC_mip;
    My_solution Approximate_AC;
    bool initial_solved_AC;
    
    
    string filename;
    bool initialized;
    bool initialized_extra;
    bool initialized_mipCP;
    
    
    
    bool FP_initialized;

    int ncols;
    unsigned seed;
    My_objective_function original_objective_function;
    My_objective_function auxilary_objective_function;

    My_objective_function FP_objective_function;

    My_solution LP_optimum;
    My_solution best_incumbent;
    vector<int> column_types;
//     vector<int> IntVars;

    int index_of_original_objective_cutoff_constraint;
    double value_of_original_objective_cutoff_constraint_ub;
    double value_of_original_objective_cutoff_constraint_lb;

    bool initial_solved;
    bool initial_solved_extra;
    bool initial_solved_mipCP;
    bool initial_solved_mipAC;
    
// 	function <int(int nm, vector<double> &zikkim, int perturbation_level)> sampler;
    double pseudo_bound;

    FP_additional_information fp_status;
    FP_information fp_inf;
    FP_history fp_hist;
	FP_history hist_v1;

    vector<double> original_lbs;
    vector<double> original_ubs;
    
    
    vector<int> integer_columns_in_slave;
    
	int FP_INT_PARAM[ENM_FP_INT_PARAM_END];
	double FP_DBL_PARAM[ENM_FP_DOUBLE_PARAM_END];
    ofstream rounded_out; 
    ofstream relaxed_out; 

    

//     ofstream rounded_out ("roundedIntVars.txt",ios_base::app);
//     ofstream relaxed_out ("relaxedIntVars.txt",ios_base::app);
//     function < int (double &d_send_recv, bool send_recv, int to_from) > communication_function;

public:
    int original_objective_direction;

    int  initialize(int solver_type = 2); 
                        /** solver type = 1 use CLP, = 2 use Cplex  
                         * returns 1 objective improvent is multiple of 1,
	                    *  i.e., all objective coefficients are integer
                        *        all variables are integer or binary*/
    cpfp_executable();
    cpfp_executable ( string arg_filename );
    cpfp_executable ( string arg_filename, unsigned seedx );
    void set_seed ( unsigned arg );
    unsigned get_seed() const;
    void apply_cuts ( vector<mySerializableRowCut> cut_list,std::ostream &out = std::cout );
    int calculate_original_objective_value ( double &retval, const  double *solution );
    int calculate_original_objective_value_of_solution ( double &retval, My_solution &solution );

    int calculate_obj_of_a_solution(My_objective_function &obj, My_solution &sol);

    int f_type_stage_2_set_bounds(My_solution &rd);
        
    void set_objective_function ( const My_objective_function &new_obj, int which_mip );
    void set_FP_objective_function ( const My_objective_function &new_obj );

    void update_objective_cutoff_constraint_ub ( double arg );
    void update_objective_cutoff_constraint_lb ( double arg );

    double get_objective_cutoff_constraint_ub() const;
    double get_objective_cutoff_constraint_lb() const;

    int optimize(int which_mip = 1); 
                        /** if which_mip = 1 solve mip if 2 solve extra_mip 3 for mipCP 4 for AC_mip
                         * returns 1 if optimal
                                0 if infeasible
                                -2 if DualObjectiveLimitReached
                                -3 if IterationLimitReached
                                -4 if PrimalObjectiveLimitReached
                                -5 if ProvenDualInfeasible
                                -6 if Abandoned*/

    double getObjValue(int which_mip = 1);
    const double * getColSolution();
    vector<double> getColSolution_asVector();
    int  find_obj_value_of_a_solution ( double &retval, const My_solution &sol, const My_objective_function &obj );
    int post_iterate ( vector<double> &starting_solution );

    int get_ncols();

    vector<int> get_column_types() const;
    My_solution get_best_incumbent();


    My_objective_function get_auxilary_objective_function() const ;
    My_objective_function get_original_objective_function() const ;

    double get_best_incumbents_objective();

    double get_integer_infeasibility ( const double * solution );

    int get_solution ( My_solution &msolution, double &integer_infeasibility );
    int get_solution_stage ( My_solution &msolution, double &integer_infeasibility, int stage );


    int get_fractionality_of_solution( My_solution &msolution, double &double_infeasibility );
    int get_binary_fractionality_of_solution( My_solution &msolution, double &double_infeasibility );

    OsiCuts generate_cuts ( Generate_cut_type in_cut_type_to_generate );


    void fp_inf_itercount_plus();
    void reset_fp_inf();
    void reset_fp_lim();

#ifdef USE_SET

    set<mySerializableRowCut> generate_rowcutlist ( Generate_cut_type in_cut_type_to_generate );
#else
    vector<mySerializableRowCut> generate_rowcutlist ( Generate_cut_type in_cut_type_to_generate, vector<mySerializableRowCut> &allreceivedcuts, int &previousnumber_of_cuts );
#endif

    int init_FP_stage2(); /** does not   set objective */
    int init_FP_stage1(); /** does not   set objective */
    int undo_FP_stage2(); /** does not reset objective */
    int undo_FP_stage1(); /** does not reset objective */
//     int round(vector<double> &r, const vector<double> i);

    int pre_FP ( int iterlim, int resetlim, double timelim, int rounding_type,vector<double> &starting_solution, vector<double> &rounded_solution );
    int pre_FP_v1 ( int iterlim, int resetlim, int st2_iterlim, int st2_resetlim,  int rounding_type,vector<double> &starting_solution, vector<double> &rounded_solution );
	int pre_FP_v1_sol ( int iterlim, int resetlim, int st2_iterlim, int st2_resetlim,  int rounding_type, My_solution &starting_solution, My_solution &rounded_solution );

    int iterate_FP ( vector<double> &starting_solution, vector<double> &rounded_solution, int rounding_type ); /** return >0 if a feasible solution is found */

    /** rounding_type   = 0; standart : threshold = 0.5 ;
     *                  = 1; completely random threshold ;
     *                  = 2; random threshold .25-.75 ;
     *                  = 3; random flip;
     *                  = 4; aggressive restart;
     *
     *                  = 9; random method ;
     */


    int intoptimize();
    double r_gen();


    void puke_lp_file ( const char *filename );
    void puke_mps_file ( const char *filename );


    void print_smt_to_somewhere ( std::ostream &out = std::cout );
    //     int run_FP_slave(int iterlim, int resetlim, double timelim);



    int create_starting_solution();




    /** basically round it */
    int round(vector<double> &starting_solution, vector<double> &rounded_solution, int arg_rounding_type, int stage = 2);

    int round_sol(My_solution &starting_solution, My_solution &rounded_solution, int arg_rounding_type, int stage = 2);
//     int round_sol_v2(My_solution &starting_solution, My_solution &rounded_solution, int &arg_rounding_type, int stage = 2);

    int f_type_stage_1_set_obj(vector<double> &rounded_solution);
    int f_type_stage_1_set_obj(My_solution &my_sol_rounded_solution);
    int f_type_stage_2_set_obj(vector<double> &rounded_solution);
    int f_type_stage_2_set_obj(My_solution &my_sol_rounded_solution);

    My_objective_function create_FP_dist_obj_for_solution_stage_1(My_solution &my_sol_rounded_solution);

    My_objective_function create_FP_dist_obj_for_solution_stage_2(My_solution &my_sol_rounded_solution);


    My_objective_function f_type_stage_1_create_obj_for_FP(My_solution &my_sol_rounded_solution, double alpha);


    My_objective_function f_type_stage_2_create_obj_for_FP(My_solution &my_sol_rounded_solution, double alpha);
    
    

    My_objective_function f_type_stage_1_or_2_create_obj_for_FP(My_solution &my_sol_rounded_solution, double alpha, int stage);
    
    


    int f_type_stage_1_check_improvement();
    int f_type_stage_1_iterate();


    int f_type_stage_2_check_improvement();
    int f_type_stage_2_iterate();




    int f_type_restart_1 ( vector<double> &starting_solution, vector<double> &rounded_solution, vector<double> &previous_rounded );
    int f_type_restart_1 ( My_solution &starting_solution, My_solution  &rounded_solution, My_solution  &previous_rounded );


    int f_type_restart_2 ( vector<double> &starting_solution, vector<double> &rounded_solution, vector<double> &previous_rounded, int iter );
    int f_type_restart_2 ( My_solution &starting_solution, My_solution  &rounded_solution, My_solution  &previous_rounded ,int iter);

    int f_type_restart_1_or_2 ( My_solution &starting_solution, My_solution  &rounded_solution, My_solution  &previous_rounded, int stage, int iter );

    int aggressive_restart(My_solution &rounded_solution);
    int aggressive_flip ( My_solution& rounded_solution, double flip_probability  );
    
    
    
    int mip_set_FP_objective_stage_1(My_solution &my_sol_rounded_solution, double alpha);
    int mip_set_FP_objective_stage_2(My_solution &my_sol_rounded_solution, double alpha);


    int set_bounds_FP_stage_2(My_solution &my_sol_rounded_solution);
    int unset_bounds_FP_stage_2();

    int get_number_of_binaries() const;
    int get_number_of_general_integers() const;

    int FP_history_check(My_solution &rounded_solution, double alpha);
    int FP_history_check(My_solution &rounded_solution, double alpha, int stage);

    void clear_FP_history(unsigned int keep_latest_n = 0);
    int add_to_history(int nfrac, double d_frac, My_solution sol);
    int compare_last_with_previous ( bool use_hist);

    int get_history(My_solution &z, int n);

    int stage_decider(int previous_stage);

    void print_solution_wth_column_types(My_solution &sol,ostream &out = std::cout);


//     int add_to_hist_v1(My_solution rounded_sol, double alpha);
//     int  get_history_v1(My_solution &z, int n);
//     void clear_hist_v1(unsigned int keep_latest_n = 0);
    /** assumes FP_bounds are not set ??*/
    int set_pseudo_bounds(double sb);
    double calculate_pseudo_bound();
    void print_bounds();


    /** assume bounds  are unset*/
    /** fix integer variables and minimize original objective value*/

    int polish_integer_solution(My_solution &sol, My_solution &new_sol, double &new_obj);


    int get_FP_INT_PARAM(FP_INT_PARAM_SET param, int &arg) const ;
    int get_FP_DOUBLE_PARAM(FP_DOUBLE_PARAM_SET param, double &arg) const ;

    int  set_FP_INT_PARAM(FP_INT_PARAM_SET param, int arg) ;
    int  set_FP_DOUBLE_PARAM(FP_DOUBLE_PARAM_SET param, double arg);


    int getNextIntegerPoint(My_solution &cont_sol, My_solution &int_sol, int stage);
    
	int getDivergentPoint(My_solution &cont_sol, My_solution &int_sol, int stage, int slave_id);
	
	
    int set_bounds_extra_mip(My_solution &s, int stage);

    int optimize_extra_mip( My_solution& s, int stage );
    
    My_solution get_solution_of_extra_mip(int stage);
    
    void print_FP_parameters(ostream &out = std::cout) const;
    
    double get_lambda_start( My_solution& s_sol, My_solution& t_sol,int index);    
    double get_lambda_min( My_solution& s_sol, My_solution& t_sol, int index, double &abs_one_over_dif);    
	
    int FP_line_search_fill_lambda_dif (My_solution &start_sol, My_solution &end_sol, double alpha_min, double aplha_max);
    
    /**  updates lambda file*/
    int FP_line_search_get_next_solution(My_solution &to_be_updated);

    //     int FP_line_search_get_next_solution(My_solution &to_be_updated);
};

typedef pair<double, int> dipair;   // (sigma, index) pair and comparison function
class dipairCmp {
    public: bool operator() (const dipair& p1, const dipair& p2) const {return p1.first>p2.first;}
};

#endif
// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
