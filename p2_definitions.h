#include "myserializablerowcut.h"
#include "my_solution.h"
#include "my_objective_function.h"
#include "cpfp_executable.h"
#include "others.h"
#include <random>

#define WAITDEBUGGING

#define SLEEP_TIME_BETWEEN_FP_PROBES 1000
#define SLEEP_TIME_BETWEEN_AC_PROBES 10000
#define FILENAME_FOR_SLAVES "slave_log_"

#define LOG_SET 1
//1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18

#define COUT_SET
//2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18


#ifdef USING_MPI
    #define SIGNAL_ALL_SLAVES_INT(TAG,NUMBER)\
    for(int i=1; i < MPI::COMM_WORLD.Get_size();++i) {\
        master_mpi_signaler.send_blocking_signal_as_tag_and_number(TAG,NUMBER,i);\
    }

    #define SIGNAL_ALL_SLAVES_DOUBLE(TAG,DOUBLE)\
    for(int i=1; i < MPI::COMM_WORLD.Get_size();++i) {\
        master_mpi_signaler.send_blocking_signal_as_tag_and_double(TAG,DOUBLE,i);\
    }
    
    #define SLAVEID MPI::COMM_WORLD.Get_rank() - 1
    #define N_SLAVES MPI::COMM_WORLD.Get_size() - 1
    #define WORLDSIZE MPI::COMM_WORLD.Get_size() 

  
#else 

    #define SIGNAL_ALL_SLAVES_INT(TAG,NUMBER)\
        signaler.send_signal_i(DESTINATION_ALL_SLAVES, TAG, NUMBER); 

    #define SIGNAL_ALL_SLAVES_DOUBLE(TAG,DOUBLE) \
		signaler.send_signal_d(DESTINATION_ALL_SLAVES, TAG, DOUBLE); 
    
    #define SLAVEID signaler.get_MY_SLAVE_ID() 
    #define N_SLAVES (signaler.get_MY_WORLD_SIZE() - 1)
    #define WORLDSIZE signaler.get_MY_WORLD_SIZE()


#endif



// The multiple macros that you would need anyway [as per: Crazy Eddie]

#ifndef DEB

#define DEB_0()\
    cout << " DEBUG call at Slave " << SLAVEID << " in file " << __FILE__ << " in function "  << __func__ << " at line " << __LINE__ << endl;

#define DEB_1(A)  \
    cout << " DEBUG call at Slave " << SLAVEID << " in file " << __FILE__ << " in function "  << __func__ << " at line " << __LINE__ << " with "<<  A  << endl;
    
#define DEB_2(A,B)\
    cout << " DEBUG call at Slave " << SLAVEID << " in file " << __FILE__ << " in function "  << __func__ << " at line " << __LINE__ << " with "<<  A  << " " << B << endl;
#define DEB_3(A,B,C)\
    cout << " DEBUG call at Slave " << SLAVEID<< " in file " << __FILE__ << " in function "  << __func__ << " at line " << __LINE__ << " with "<<  A  << " " << B << " " << C <<endl;
#define DEB_4(A,B,C,D)\
    cout << " DEBUG call at Slave " << SLAVEID << " in file " << __FILE__ << " in function "  << __func__ << " at line " << __LINE__ << " with "<<  A  << " " << B << " " << C << " " << D << endl;

// The interim macro that simply strips the excess and ends up with the required macro
#define DEB_X(x,A,B,C,D,FUNC, ...)  FUNC  

// The macro that the programmer uses 
#define DEBON(...)                    DEB_X(,##__VA_ARGS__,\
                                          DEB_4(__VA_ARGS__),\
                                          DEB_3(__VA_ARGS__),\
                                          DEB_2(__VA_ARGS__),\
                                          DEB_1(__VA_ARGS__),\
                                          DEB_0(__VA_ARGS__)\
                                         ) 

  #define DEBWAIT(...)\
        int __debwait; \
        DEBON (__VA_ARGS__)\
        cin >> __debwait;
    #define DEBOFF(...)
    
#endif
    
    


// #define DEBUGGING

    
#ifdef DEBUGGING
# define __MSS(INTVALUE) cout << "MASTR SEND SIGN " << INTVALUE << " __LINE__: " << __LINE__ << endl;
# define __MSD(INTVALUE) cout << "MASTR SEND DATA " << INTVALUE << " __LINE__: " << __LINE__ << endl;
# define __MRS(INTVALUE) cout << "MASTR RECV SIGN " << INTVALUE << " __LINE__: " << __LINE__ << endl;
# define __MRD(INTVALUE) cout << "MASTR RECV DATA " << INTVALUE << " __LINE__: " << __LINE__ << endl;

# define __SSS(INTVALUE) if(create_cout()) cout << "SLAVE SEND SIGN " << INTVALUE << " __LINE__: " << __LINE__ << endl;
# define __SSD(INTVALUE) if(create_cout()) cout << "SLAVE SEND DATA " << INTVALUE << " __LINE__: " << __LINE__ << endl;
# define __SRS(INTVALUE) if(create_cout()) cout << "SLAVE RECV SIGN " << INTVALUE << " __LINE__: " << __LINE__ << endl;
# define __SRD(INTVALUE) if(create_cout()) cout << "SLAVE RECV DATA " << INTVALUE << " __LINE__: " << __LINE__ << endl;



#else


# define __MSS(INTVALUE)
# define __MSD(INTVALUE)
# define __MRS(INTVALUE)
# define __MRD(INTVALUE)

# define __SSS(INTVALUE)
# define __SSD(INTVALUE)
# define __SRS(INTVALUE)
# define __SRD(INTVALUE)

#endif

#ifdef WAITDEBUGGING

# define BWAIT \
if(blockingWait) {\
	int kz; \
	cout << "blockingWait ";\
	if (SLAVEID == 0) cout << "Master ";\
	else cout << "Slave " << SLAVEID;\
	cout << "\t line " << __LINE__ << " of file " << __FILE__ << endl;\
	cin >> kz;\
  }

#else

# define BWAIT

#endif

//#define USING_IOPTIMIZE

#ifdef USING_IOPTIMIZE
# include "start_point_generator6.h"
#endif

#define USE_POST_CP_FP
#ifdef USE_POST_CP_FP
# define USE_POST_CP_FP_iterlim 100
#endif

#undef USE_POST_CP_FP



#define DOUT \
//
