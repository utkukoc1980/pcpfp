#ifndef PARALLEL_CPFP_H
	#define PARALLEL_CPFP_H
	#include "OsiCuts.hpp"

	#include "OsiSolverInterface.hpp"
	#include "OsiSolverParameters.hpp"
	#include "OsiClpSolverInterface.hpp"
	
// #define USE__CPLEX	
	
// #ifdef USE__CPLEX
	#include "OsiCpxSolverInterface.hpp"
// #endif

	#include "CoinPackedVectorBase.hpp"

	#include "CglKnapsackCover.hpp"
	#include "CglSimpleRounding.hpp"
	#include "CglRedSplit.hpp"
	#include "CglGomory.hpp"


	#include "CglAllDifferent.hpp"
	#include "CglPreProcess.hpp"
	#include "CglLandP.hpp"
	#include "../templates.h"
	

#ifdef USING_MPI
	#include "mpi.h"
	
	#define MACRO_SLAVE_ID MPI::COMM_WORLD.Get_rank()
    #define IFMASTER if(MPI::COMM_WORLD.Get_rank() == 0 )
	#define IFSLAVE if(MPI::COMM_WORLD.Get_rank() > 0 )
	#define IFSLAVE_y(y) if(MPI::COMM_WORLD.Get_rank() == y )
	#define NSLAVES MPI::COMM_WORLD.Get_size() - 1
    
    #define CPFP_ERROR(zikkim, out)\
        out << "ERROR: " << zikkim << " (in rank " << MPI::COMM_WORLD.Get_rank() << " of line " << __LINE__ << " of file " << __FILE__ <<endl;

         
	#define COMMON_DATA_SEND_ALL \
	MPI_Send(&size,1,MPI_INT,destination,tag,MPI_COMM_WORLD);\
	MPI_Send(indices2,size,MPI_INT, destination,tag,MPI_COMM_WORLD);\
	MPI_Send(elements2,size,MPI_DOUBLE, destination,tag,MPI_COMM_WORLD);

	#define COMMON_DATA_SEND_FUNCTION \
	MPI_Send(&original_size,1,MPI_INT,destination,tag,MPI_COMM_WORLD);\
	COMMON_DATA_SEND_ALL

	#define COMMON_DATA_SEND_SOLUTION \
	COMMON_DATA_SEND_FUNCTION \
	MPI_Send(&original_objective,1,MPI_DOUBLE,destination,tag,MPI_COMM_WORLD);\
	MPI_Send(&auxiliary_objective,1,MPI_DOUBLE,destination,tag,MPI_COMM_WORLD);

	#define COMMON_DATA_SEND_ROWCUT \
	MPI_Send(&sense,1,MPI_CHAR,destination,tag,MPI_COMM_WORLD);\
	COMMON_DATA_SEND_ALL \
	MPI_Send(&norm,1,MPI_DOUBLE,destination,tag,MPI_COMM_WORLD);\
	MPI_Send(efficiency,CUT_EFFICIENCY_END,MPI_DOUBLE, destination,tag,MPI_COMM_WORLD);\
	MPI_Send(&lb,1,MPI_DOUBLE,destination,tag,MPI_COMM_WORLD);\
	MPI_Send(&ub,1,MPI_DOUBLE,destination,tag,MPI_COMM_WORLD);\

	#define COMMON_DATA_RECV_ALL \
	MPI_Recv(&size,1,MPI_INT,source,tag,MPI_COMM_WORLD,&v_status[1]);\
	indices.resize(size);\
	MPI_Recv(&indices[0],size,MPI_INT, source,tag,MPI_COMM_WORLD,&v_status[2]);\
	elements.resize(size);\
	MPI_Recv(&elements[0],size,MPI_DOUBLE, source,tag,MPI_COMM_WORLD,&v_status[3]);

	#define COMMON_DATA_RECV_FUNCTION\
	MPI_Recv(&original_size,1,MPI_INT,source,tag,MPI_COMM_WORLD,&v_status[0]);\
	COMMON_DATA_RECV_ALL

	#define COMMON_DATA_RECV_SOLUTION\
	COMMON_DATA_RECV_FUNCTION \
	MPI_Recv(&original_objective,1,MPI_DOUBLE,source,tag,MPI_COMM_WORLD,&v_status[4]);\
	MPI_Recv(&auxiliary_objective,1,MPI_DOUBLE,source,tag,MPI_COMM_WORLD,&v_status[5]);

	#define COMMON_DATA_RECV_ROWCUT\
	MPI_Recv(&sense,1,MPI_CHAR,source,tag,MPI_COMM_WORLD,&v_status[0]);\
	COMMON_DATA_RECV_ALL\
	MPI_Recv(&norm, 1, MPI_DOUBLE, source,tag,MPI_COMM_WORLD,&v_status[4]);\
	MPI_Recv(&efficiency[0], CUT_EFFICIENCY_END, MPI_DOUBLE, source,tag,MPI_COMM_WORLD,&v_status[5]);\
	MPI_Recv(&lb, 1, MPI_DOUBLE, source,tag,MPI_COMM_WORLD,&v_status[6]);\
	MPI_Recv(&ub, 1, MPI_DOUBLE, source,tag,MPI_COMM_WORLD,&v_status[7]);

	#define OBJECT_RECV(source, TAG, object)\
		req_stat_vectors rvv;\
		object._recv_mpi(source, TAG,rvv.v_req, rvv.v_stat, MPI_RECV_BLOCKING);

	#define OBJECT_SEND(destination,TAG,object)\
		req_stat_vectors rvv;\
		object._send_mpi(destination, TAG, rvv.v_req, rvv.v_stat, MPI_SEND_BLOCKING);

	#define OBJECT_SEND_ALL_SLAVES(TAG, object)\
		for(int destination = 1; destination < MPI::COMM_WORLD.Get_size();++destination){\
			OBJECT_SEND_MPI(destination, TAG, object)\
		}
    struct req_stat_vectors{
		vector<MPI_Request> v_req;
		vector<MPI_Status> v_stat;
	};
    
    #define PROBE(source, TAG)\
            MPI_Status status;\
            int t = MPI_Iprobe ( source,TAG, MPI::COMM_WORLD, &flag, &status );
            
#define RECV_D (source, TAG)  MPI_Recv ( &r_d, 1, MPI_DOUBLE,source,TAG,MPI_COMM_WORLD, &status);

#define RECV_I (source, TAG) MPI_Recv ( &r_i, 1, MPI_INT,source,TAG,MPI_COMM_WORLD, &status);    
     
#define SEND_I(destination, TAG) MPI_Send ( &s_i, 1, MPI_INT, destination,  TAG, MPI::COMM_WORLD );

#define SEND_D(destination, TAG) MPI_Send ( &s_d, 1, MPI_DOUBLE, destination,  TAG, MPI::COMM_WORLD );
//  

               
#define FINALIZE MPI_Finalize();
	
            
            
#else 
	#define MACRO_SLAVE_ID signaler.get_MY_SLAVE_ID()
    #define IFMASTER if(signaler.get_MY_SLAVE_ID() == 0)
	#define IFSLAVE if(signaler.get_MY_SLAVE_ID()  > 0 )
	#define IFSLAVE_y(y) if(signaler.get_MY_SLAVE_ID()  == y )
	#define NSLAVES (signaler.get_MY_WORLD_SIZE() - 1)
	
	#define CPFP_ERROR(zikkim, out)\
        out << "ERROR: " << zikkim << " (in rank " << signaler.get_MY_SLAVE_ID() << " of line " << __LINE__ << " of file " << __FILE__ <<endl;

	#define COMMON_DATA_SERIALIZE_ALL \
	    stringstream s; \
        s << size << " ";\
        for (unsigned i = 0; i< indices.size(); ++i)  {s << indices[i] << " " ;};\
        for (unsigned i = 0; i< elements.size(); ++i) {s << elements[i] << " ";};
	
    #define COMMON_DATA_UNSERIALIZE_ALL(text) \
        if (text == "") return; \
        std::vector<std::string> s = split_string(text, ' ');\
        size = stoi(s[0]);\
        indices.resize(unsigned(size));\
        elements.resize(unsigned(size));\
        for (unsigned i = 0; i< indices.size(); ++i) {indices[i] = stoi(s[i+1]);};\
        for (unsigned i = 0; i< elements.size(); ++i) {elements[i] = stod(s[size + i+1]);};
        
    #define COMMON_DATA_SERIALIZE_FUNCTION \
        COMMON_DATA_SERIALIZE_ALL\
        s << original_size << " ";
    
    #define COMMON_DATA_UNSERIALIZE_FUNCTION(text) \
        COMMON_DATA_UNSERIALIZE_ALL(text)\
        original_size = stoi(s[2*size +1]);
   
    #define COMMON_DATA_SERIALIZE_SOLUTION \
        COMMON_DATA_SERIALIZE_FUNCTION \
        s << original_objective << " " << auxiliary_objective << " ";
        
    #define COMMON_DATA_UNSERIALIZE_SOLUTION(text) \
        COMMON_DATA_UNSERIALIZE_FUNCTION(text)\
        original_objective = stod(s[2*size +2]);\
        auxiliary_objective = stod(s[2*size +3]);
 
    #define COMMON_DATA_SERIALIZE_ROWCUT \
        COMMON_DATA_SERIALIZE_ALL \
        s << sense << " " << norm << " " << lb << " " << ub << " ";\
        for (unsigned i = 0; i< CUT_EFFICIENCY_END; ++i) {s << efficiency[i] << " ";}; 
        
    #define COMMON_DATA_UNSERIALIZE_ROWCUT(text) \
        COMMON_DATA_UNSERIALIZE_ALL(text)\
        sense = (s[2*size +1]).c_str()[0]; \
        norm = stod(s[2*size +2]);\
        lb =  stod(s[2*size +3]);\
        ub =  stod(s[2*size +4]);\
        for (unsigned i = 0; i< CUT_EFFICIENCY_END; ++i) {efficiency[i] = stod(s[2*size +5+i]);}; 
    
   	#define OBJECT_RECV(SOURCE, TAG, object)\
		object._recv_alt1(SOURCE, TAG,signaler);

	#define OBJECT_SEND(DESTINATION,TAG,object)\
		object._send_alt1(DESTINATION, TAG,signaler);

	#define OBJECT_SEND_ALL_SLAVES(TAG, object)\
		object._send_alt1(DESTINATION_ALL_SLAVES,TAG,signaler);

    #define PROBE(SOURCE, TAG)\
        flag = signaler.check_signal(SOURCE,TAG);
        
    #define RECV_D(SOURCE, TAG)\
        int ___s =1;\
        r_d = signaler.receive_signal_d(SOURCE,TAG,___s);\
        if (___s == TAG_SIGNALER_RETURNS_ERROR) r_d = TAG_SIGNALER_RETURNS_ERROR;
        
        
    #define RECV_I(SOURCE, TAG)\
        int ___s =1;\
        r_i = signaler.receive_signal_i(SOURCE,TAG,___s);\
        if (___s == TAG_SIGNALER_RETURNS_ERROR) r_i = TAG_SIGNALER_RETURNS_ERROR;
        
    
    #define SEND_I(DESTINATION, TAG) \
        signaler.send_signal_i(DESTINATION,TAG,s_i);

    #define SEND_D(DESTINATION, TAG) signaler.send_signal_d(DESTINATION,TAG,s_d);
//  

        
        

    #define FINALIZE 
#endif


	#define COMPARE_LAST_N 3

	#define __CUT_EQUALITY_TOLERANCE 1e-6

	#define SLEEP_TIME_BETWEEN_SIGNALS 0.001
	#define SLEEP_TIME_BETWEEN_OBJECTS 0.001


	/** TAGLIST */
    #define DOUBLE_TAG_LIMIT 1000 //ANY TAG greater than this value corresponds to a double information 
	
    #define TAG_TESTING 99
    

    
    #define DESTINATION_ALL_SLAVES -111
	#define TAG_SENDING_OBJECTIVE_FUNCTION 100
	#define TAG_SENDING_SOLUTION 101
	
    #define TAG_SIGNALER_RETURNS_ERROR -10999
	#define TAG_SIGNALER_RETURNS_NOERROR -10199
	#define TAG_SIGNALER_NO_NEW_SIGNAL_TIME_LIMIT -10201
	#define TAG_SIGNALER_NO_NEW_SIGNAL -10200
	#define TAG_SIGNALER_NEW_SIGNAL -10000
	#define TAG_SIGNALER_PROBE_SUCCESSFUL 1234
	
	
	#define TAG_SENDING_ROWCUT 102
	#define TAG_SENDING_GENERATE_CUT_TYPE 103
	#define TAG_SENDING_OBJECTIVE_VALUE 1104
	#define TAG_SENDING_OBJECTIVE_VALUE_IN_FP 1105 
	#define TAG_SENDING_ANALYTIC_CENTER 106
	#define TAG_ALL_DONE_WITH_TIME_LIMIT 107
	
	#define TAG_SENDING_LOWERBOUND 1108
    #define TAG_SENDING_BESTOBJECTIVE 1109
	
	#define TAG_SIGNAL_COMMAND 200
	#define TAG_CHANGE_STATUS_TO_CP 202
	#define TAG_CHANGE_STATUS_TO_FP 203
	#define TAG_CHANGE_STATUS_DONE_WITH_FP 204
	#define TAG_CHANGE_STATUS_NEXT_FP_ITERATION 205
	#define TAG_SENDING_DONE_MINIMUM_NUMBER_OF_ITERATIONS 206
	#define TAG_ALL_DONE_MINIMUM_NUMBER_OF_ITERATIONS 207
	#define TAG_SENDING_ALL_DONE_FINAL_NUMBER_OF_ITERATIONS 208
	#define TAG_SENDING_DONE_WITH_TIME_LIMIT 209
	
	
    #define FOUND_A_FEASIBLE_SOLUTION_IN_FP 1301

	#define TAG_SENDING_SIGNAL 800
	#define TAG_ABORT 	999

	#define SIGNAL_SENDING_OBJECTIVE_FUNCTION 201
// 	#define SIGNAL_SENDING_OBJECTIVE_FUNCTION 201


	#define   my_default_seed 5489

	#define CUT_EFFICIENCY_K 1.0
	
	
	
	using namespace std;

	enum COMMAND_LIST_ENM{
		COMMAND_LIST_BEGIN = 222000, /* begin tag for command list enm*/
		COMMAND_LIST_DO_OBJ_FP_OR = 222001,
		COMMAND_LIST_DO_OBJ_FP_OR_ALT_ROUNDING = 222002,
		COMMAND_LIST_DO_OBJ_FP_RW = 222003,
		COMMAND_LIST_DO_OBJ_FP_RR = 222004,
		COMMAND_LIST_DO_OBJ_FP_ORIGINAL = 222005,
		COMMAND_LIST_DO_CP = 222006,
		COMMAND_LIST_END = 222007,
		/*
		Command	ALG	START POINT		RESTART 	ROUNDING
		FP	OBJ	LPOPT		back to LP opt new rounding	1
		BASIC	RR opt	perturbation amount	new RR	2
		PR opt		new PR	3
		R interior		new R int	4
		P Interior		new P int	
		LATTICE WALK		new LW	
		LP opt NEIGH			
		
		
		
		CP	GOMORY	LPOPT			
		KNAPSACK	PR opt			
		.	RR opt			
		.				
		
		
		*/
	};
	
	enum CUT_EFFICIENCY_ENM
	{
        CUT_EFFICIENCY_BGN = 0, /**< begin tag for cut efficiency >*/
		CUT_EFFICIENCY_NORM_ALPHA = 1, /** ||a|| */
		CUT_EFFICIENCY_NORM_C_ORG = 2, /** || original_objective|| */
		CUT_EFFICIENCY_NORM_C_AUX = 3, /** || auxiliary_objective|| */
		CUT_EFFICIENCY_VIOLATION = 4, /** v(a,b,x*) = a^T x*  - b */
		CUT_EFFICIENCY_RELATIVE_VIOLATION = 5, /** v(.) / |b| */
		CUT_EFFICIENCY_DISTANCE = 6, /** d() = v(.) / ||a|| */
		CUT_EFFICIENCY_ADJUSTED_DISTANCE = 7, /** ad() = v(.) / (||\hat(a)|| + 1) where \hat(a)_j = {a_j if x*_j >0; 0 ow } */
		CUT_EFFICIENCY_DISTANCE_VARIANT = 8, /** dv() = v(.)^k / (|a_1 . a_2 . ... . a_r|)^(k/r). take k = 1 */

		CUT_EFFICIENCY_OBJECTIVE_PARALLELISM_WRT_ORG = 9, /**o(a) = a^T c / ||a||.||c|| */
		CUT_EFFICIENCY_OBJECTIVE_PARALLELISM_WRT_AUX = 10, /**o(a) = a^T c / ||a||.||c|| */
		CUT_EFFICIENCY_EXPECTED_IMPROVEMENT_WRT_ORG  = 11, /** ||c|| o(a) d() */
		CUT_EFFICIENCY_EXPECTED_IMPROVEMENT_WRT_AUX  = 12, /** ||c|| o(a) d() */
													/** NOTE: need to consider which c to use here */

		CUT_EFFICIENCY_SUPPORT = 13, /** 1-s(a), s(a) = |\bar(N)| / |N|
										* where \bar(N) = {j \in N: |a_j| > 0} */

		CUT_EFFICIENCY_INTEGRAL_SUPPORT = 14, /** i(a) = |\bar(N)_I| / |\bar(N)|
												* where
												* \bar(N)_I = \bar(N)I intersection N_I */
		// 		CUT_EFFICIENCY_ROTATED_DISTANCE = 6, /** v(.) / ||\bar(a)||  where \bar(a) = a -D^T(DD^T)^(-1)D a */
		// 		CUT_EFFICIENCY_DISTANCE_WITH_BOUNDS = 7, /** min(||x-x*||:a^T x<= b, l <= x <= u) */
		// 		CUT_EFFICIENCY_ROTATED_DISTANCE_WITH_BOUNDS = 8, /** min(||x-x*||:\bar(a)^T x<= \bar(b), l <= x <= u) */


        CUT_EFFICIENCY_END = 15, /**< end tag for cut efficiency >*/
    };


	enum CUT_TYPE_ENM
	{
		CUT_TYPE_BGN                   =    0, /**< begin tag of the cut type */
		CUT_TYPE_GOMORY 				=    0, /**< use gomory cuts */
		CUT_TYPE_KNAPSACK				= 	 1, /**< use knapsack cuts */
		CUT_TYPE_REDSPLIT				= 	 2, /**< use reducsed split cuts */
		CUT_TYPE_SIMPLE_ROUNDING		= 	 3, /**<use simple rounding cuts */
        CUT_TYPE_ALL_DIFFERENT 			= 	 4,
		CUT_TYPE_L_AND_P				= 	 5,
		CUT_TYPE_LIFT_AND_PROJECT 		= 	 6, 
		CUT_TYPE_END					=    7, /**<end tag for the cut type */
	};


	enum OBJECTIVE_TYPE_ENUM
	{
		OBJECTIVE_TYPE_BGN 								= 0,  /**< begin tag of objective type*/
		OBJECTIVE_TYPE_ORIGINAL 						= 0,  /**< use original objective function*/
		OBJECTIVE_TYPE_RANDOM_FROM_SHORT_DIKIN			= 1,  /**< use a random function generated from the analytic center to a random point (random point by short Dikin Walks )*/
		OBJECTIVE_TYPE_RANDOM_FROM_LONG_DIKIN			= 2, /**< use a random function generated from the analytic center to a random point (random point by long Dikin Walks )*/
		OBJECTIVE_TYPE_RANDOM_FROM_HIT_RUN				= 3, /**< use a random function generated from the analytic center to a random point (random point by shit and run )*/
		OBJECTIVE_TYPE_PERTURBED_FROM_SHORT_DIKIN		= 4, /**< use a poerturbed function generated from the analytic center to a random point (random point by short Dikin Walks )*/
		OBJECTIVE_TYPE_PERTURBED_FROM_LONG_DIKIN		= 5,/**< use a poerturbed function generated from the analytic center to a random point (random point by long Dikin Walks )*/
		OBJECTIVE_TYPE_PERTURBED_FROM_HIT_RUN			= 6,/**< use a poerturbed function generated from the analytic center to a random point (random point by hit and run )*/
		
		OBJECTIVE_TYPE_RANDOM_WALK_SHORT_DIKIN          = 7,  /**< use a random point generated from the analytic center to a random point (random point by short Dikin Walks )*/
		OBJECTIVE_TYPE_RANDOM_WALK_LONG_DIKIN           = 8, /**< use a random point generated from the analytic center to a random point (random point by long Dikin Walks )*/
		OBJECTIVE_TYPE_RANDOM_WALK_HIT_RUN              = 9, /**< use a random point generated from the analytic center to a random point (random point by shit and run )*/
		
		
		
		
		OBJECTIVE_TYPE_END								= 10, /**<end tag for objective type*/
	};



	enum OBJECTIVE_PERTURBATION_ENUM{ /** percent value */
		OBJECTIVE_PERTURBATION_BGN 			=0,

		OBJECTIVE_PERTURBATION_ZERO 		=0,
		OBJECTIVE_PERTURBATION_ONE 		=1,
		OBJECTIVE_PERTURBATION_FIFTY		=50,
		OBJECTIVE_PERTURBATION_NINENINE		=99,
		OBJECTIVE_PERTURBATION_HUNDRED		=100,

		OBJECTIVE_PERTURBATION_END 		=101,

	};

	enum OBJECTIVE_SENSE {MAXIMIZE, MINIMIZE};

	enum MPI_SEND_TYPE{MPI_SEND_BLOCKING =1, MPI_ISEND_NONBLOCKING = 2,};
    enum MPI_RECV_TYPE{MPI_RECV_BLOCKING =1, MPI_IRECV_NONBLOCKING = 2,};

    enum ALT1_SEND_TYPE{ALT1_SEND_BLOCKING =1, ALT1_SEND_NONBLOCKING = 2,};
    enum ALT1_RECV_TYPE{ALT1_RECV_BLOCKING =1, ALT1_RECV_NONBLOCKING = 2,};



    
    

#endif

	/* This counts the number of args */
	#define NARGS_SEQ(_1,_2,_3,_4,_5,_6,_7,_8,N,...) N
	#define NARGS(...) NARGS_SEQ(__VA_ARGS__, 8, 7, 6, 5, 4, 3, 2, 1)
	
	/* This will let macros expand before concating them */
	#define PRIMITIVE_CAT(x, y) x ## y
	#define CAT(x, y) PRIMITIVE_CAT(x, y)
	
	/* This will pop the last argument off */
	#define POP_LAST(...) CAT(POP_LAST_, NARGS(__VA_ARGS__))(__VA_ARGS__)
	#define POP_LAST_1(x1)
	#define POP_LAST_2(x1, x2) x1
	#define POP_LAST_3(x1, x2, x3) x1, x2
	#define POP_LAST_4(x1, x2, x3, x4) x1, x2, x3
	#define POP_LAST_5(x1, x2, x3, x4, x5) x1, x2, x3, x4
	#define POP_LAST_6(x1, x2, x3, x4, x5, x6) x1, x2, x3, x4, x5
	#define POP_LAST_7(x1, x2, x3, x4, x5, x6, x7) x1, x2, x3, x4, x5, x6
	#define POP_LAST_8(x1, x2, x3, x4, x5, x6, x7, x8) x1, x2, x3, x4, x5, x6, x7
	
	/* This will return the last argument */
	#define LAST(...) CAT(LAST_, NARGS(__VA_ARGS__))(__VA_ARGS__)
	#define LAST_1(x1) x1
	#define LAST_2(x1, x2) x2
	#define LAST_3(x1, x2, x3) x3
	#define LAST_4(x1, x2, x3, x4) x4
	#define LAST_5(x1, x2, x3, x4, x5) x5
	#define LAST_6(x1, x2, x3, x4, x5, x6) x6
	#define LAST_7(x1, x2, x3, x4, x5, x6, x7) x7
	#define LAST_8(x1, x2, x3, x4, x5, x6, x7, x8) x8
	
	#define MPI_SEND_BLOCK_SWITCH(...) { \
		if(SEND_TYPE == MPI_SEND_BLOCKING) {MPI_Send(POP_LAST(__VA_ARGS__));} \
		if(SEND_TYPE == MPI_ISEND_NONBLOCKING) {MPI_Isend(__VA_ARGS__);} \
			}
			

    #define TAG_CONTROL(TAG1)\
    if( tag != TAG1 ) { \
		if( tag != TAG_TESTING ) { \
			cout << "ERROR: tag: " << tag << " BUT should be TAG1" << TAG1 << ": " <<  TAG1 << endl;}}
	
// 	#define IFBLOCKING if(MPI_SEND == MPI_SEND_BLOCKING)
// 	#define IFNONBLOCKING  if (MPI_SEND == MPI_ISEND_NONBLOCKING)
    


#ifndef MPI_TEST_FUNCTIONS
#ifdef USING_MPI
	#define MPI_TEST_FUNCTIONS


	template <class Test_class>
	void test_(Test_class &test_object, bool send_operator_equal = true, bool send_a_copy= true){
		int my_id = MPI::COMM_WORLD.Get_rank();
		int world_size = MPI::COMM_WORLD.Get_size();
		req_stat_vectors r;

		if(my_id == 0){
			cout << "MASTER SENDS                : ";
			test_object.print();
			for(int i =1; i < world_size;++i ){
// 				cout << "MASTER: sending " << i << endl;

				test_object._send_mpi(i,TAG_TESTING, r.v_req, r.v_stat);
// 				cout << "MASTER: sending  " << i << " done."<< endl;

			}
// 			cout << "__LINE__  " << __LINE__  << endl;

			vector<Test_class> asd;
			if (send_operator_equal) asd.resize(world_size);
			vector<Test_class> asd2;
			if (send_a_copy) asd2.resize(world_size);

			for(int i =1; i < world_size;++i ){
// 				cout << "MASTER: receiving " << i << endl;
				if (send_operator_equal) {
					asd[i]._recv_mpi(i,TAG_TESTING, r.v_req, r.v_stat);
				}
				if (send_a_copy){ 
					asd2[i]._recv_mpi(i,TAG_TESTING, r.v_req, r.v_stat);
				}
// 				cout << "MASTER: receiving " << i << " done."<< endl;

			}
			if (send_operator_equal){
				for(int i =1; i < world_size;++i ){
					cout << "MASTER RECEIVES =(slave " << i << ")  : ";
					asd[i].print();
				}
				
			}
			if (send_a_copy){
				for(int i =1; i < world_size;++i ){
					cout << "MASTER RECEIVES CP(slave " << i << ") : ";
					asd[i].print();
				}
			}
		}
		else{
// 			cout << "SLAVE " << my_id <<  " active" << endl;
			test_object._recv_mpi(0,TAG_TESTING, r.v_req, r.v_stat);
// 			COUT_TESTER(__LINE__)
			
			cout << "SLAVE "<< my_id <<  " RECEIVES            : ";
			test_object.print();
// 			COUT_TESTER(__LINE__)
			
			if (send_operator_equal){
				Test_class t = test_object;
				/**
				* DO SOMETING HERE
				* */
				cout << "SLAVE "<< my_id <<  " SENDS AN EQUAL      : ";
				t.print();

				t._send_mpi(0,TAG_TESTING, r.v_req, r.v_stat);
			}
			if (send_a_copy){
				Test_class t_copy(test_object);
				cout << "SLAVE "<< my_id <<  " SENDS A COPY        : ";
				t_copy.print();
				t_copy._send_mpi(0,TAG_TESTING, r.v_req, r.v_stat);
			}
		}
		return;
	}

	template <class Test_class>
	void test_default(bool send_operator_equal = true, bool send_a_copy= true){
		Test_class my_test_object;
		test_<Test_class>(my_test_object, send_operator_equal, send_a_copy);

	}
#endif
#endif

#ifdef USING_MPI
#ifndef SEND_AND_RECV_SIGNAL_AS_TAG_AND_NUMBER
	#define SEND_AND_RECV_SIGNAL_AS_TAG_AND_NUMBER

	class EASY_SEND_AND_RECEIVE{
        int send_recv_tag;
		int number_tag;
		double double_tag;


		vector<MPI_Request> v_request;
		vector<MPI_Status> v_status;


        stringstream signal_sender_filename;
        stringstream signal_receiver_filename;
        void EASY_Recv();
    public:
		int receive_blocking_signal_as_tag_and_number(int source);
		int receive_blocking_signal_as_tag_and_double(int source);
		int send_blocking_signal_as_tag_and_number(int tag, int number, int destination);
		int send_blocking_signal_as_tag_and_double(int tag, double d_number, int destination);
		int get_tag() const;
		int get_number() const;
		double get_double() const;
		void print(std::ostream &out = std::cout) const;


	};


#endif
#endif
    
#ifndef MY_SIGNALER_V1
#define MY_SIGNALER_V1
    
    class signal;
    class object;
    // the structure of the input_output file system
    class my_signaler{
        Timer ts; 
		Timer to; 
        vector<signal> signal_list;
    	vector<object> object_list;
        
    	unsigned latest_signal_line;
    	unsigned latest_object_line;

        int MY_SLAVE_ID; // uniq slave id, 0 'if master'
        int MY_WORLD_SIZE; // number of slaves + 1 (master) similar to MPI_COMM_WORLD.Get_Size

        int send_recv_tag; // latest send or received tag (signal or data)
		int number_tag; // The latest received / sent integer signal 
		double double_tag; // the latest received / sent double signal 

        stringstream latest_received_objectstringstream; // the latest received / sent object -- serialized 

        // object is serialized stringstream of cut, objective function, or solution
        
        ofstream sender; // the filestream to write data (object)
        ifstream receiver; // the filestream to read data (object)
        
        string communication_world_uniq_file_prefix; // uniq prefix for the run (“in case multiple runs of different files happen at the same time)

        stringstream object_sender_filename;  // the name of the file to write (send) data (object )
        stringstream object_receiver_filename; // the name of the file to read (recv) data (object )
        
        ofstream signal_sender; // the filestream to write signal (int or double )
        ifstream signal_receiver; // the filestream to read signal (int or double )

        stringstream signal_sender_filename;   // the name of the file to write (send) signal (int or double)
        stringstream signal_receiver_filename; // the name of the file to read (recv) signal (int or double)
        
        int signal_check(double ms = SLEEP_TIME_BETWEEN_SIGNALS);
		int object_check(double ms = SLEEP_TIME_BETWEEN_OBJECTS);
        
        /* 
        
        filename convention =
     	slave reads (recv) signals from prefix_slave_MYSLAVEID_signal.in  // 
		slave writes (send) signals to prefix_slave_MYSLAVEID_signal.out 

     	slave reads (recv) objects (data) from prefix_slave_MYSLAVEID_data.in 
		slave writes (send) objects (data) to prefix_slave_MYSLAVEID_data.out 
		
		file systems 		

		*_signal.in
		source TAG value \n  
		source is the slave_ID of the received signal (0 for master), TAG is the signal info and value is either integer or double (based on TAG)

		*_signal.out

		destination TAG value \n  
		destination is the destination slave_ID for signal (0 for master), TAG is the signal info and value is either integer or double (based on TAG)
		if destination is 'A' in prefix_MASTER_signal.out then send to all slaves, 	
		
		*_data.in

		source TAG data \n 
		source is the slave_ID of the received data (0 for master), TAG is the datatype info and data stringstream for object solution, function,AC or cut (based on TAG)

		*_data.out

		destination TAG data \n
		destination is the destination slave_ID for data (0 for master), TAG is the datatype info and data stringstream for object solution, function,AC or cut (based on TAG)
		if destination is 'A' in prefix_MASTER_data.out then send to all slaves, 	

        */
        
    public:
	
    	my_signaler();
		my_signaler(string s, int ID, int WS);
		my_signaler(const my_signaler &other);
        my_signaler &operator= (const my_signaler& other);
    	int receive_signal_as_tag_and_number(int source, int &flag,ALT1_RECV_TYPE RECV_TYPE = ALT1_RECV_NONBLOCKING ); //-1 for any source
		double receive_signal_as_tag_and_double(int source, int &flag, ALT1_RECV_TYPE RECV_TYPE = ALT1_RECV_NONBLOCKING ); //-1 for any source

        void clean_files(bool sig = true, bool obj=true);
		
		string recv_object (int source, int tag, ALT1_RECV_TYPE RECV_TYPE = ALT1_RECV_NONBLOCKING ); //-1 for any source

		void send_signal_d(int destination, int tag, double d_number){int retval = send_signal_as_tag_and_double(destination, tag, d_number); return;};
		void send_signal_i(int destination, int tag, int number){int retval = send_signal_as_tag_and_number(destination, tag, number);return;};

		int  send_signal_as_tag_and_number(int destination, int tag, int number); 
		int  send_signal_as_tag_and_double(int destination, int tag, double d_number);
        
		int send_object ( int destination, int tag, string objs, ALT1_SEND_TYPE SEND_TYPE = ALT1_SEND_NONBLOCKING);
		
		int get_tag() const;
		int get_number() const;
		double get_double() const;

        void set_tag(int tag);
		//int get_latest_received_objectstring() const;
		
		void print(std::ostream &out = std::cout) const;
        int get_MY_SLAVE_ID()const; 
        int get_MY_WORLD_SIZE() const;
        
		
		int check_signal(int source, int tag); //returns 1 if there is a new signal tag -1 for any signal 
        int receive_signal_i(int source, int tag, int &flag,ALT1_RECV_TYPE RECV_TYPE = ALT1_RECV_NONBLOCKING );
        double receive_signal_d(int source, int tag, int &flag,ALT1_RECV_TYPE RECV_TYPE = ALT1_RECV_NONBLOCKING );
        
        string generate_slaveX_signal_in(int sl_id);
        string generate_slaveX_object_in(int sl_id);
        
        
    };
    
      class signal{
        int source; 
        int tag;
        int intval;
        double doubleval;
        bool used; 
        bool intmi;
        friend class my_signaler;
      public:
        signal();
        signal(int s, int t, int i);
        signal(int s, int t, double d);
        void print();
    };
    class object{
        
        int source; 
        int tag;
        string object_string;
        bool used; 
        friend class my_signaler;
    public:
        object();
        object(int s, int t, string ss);
    };
    
#endif
    

