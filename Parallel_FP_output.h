
#ifndef Parallel_FP_output_H
    #define Parallel_FP_output_H
	#include "parallel_cpfp_stuff.h"
    #include "my_solution.h" 
    using namespace std;
    
    class Parallel_FP_output{
        
       
        string lpfile;
        char OVPR;
        My_solution sol;
        double obj;
        
        
    public:
        Parallel_FP_output(bool found2,double obj2, My_solution sol2){
            obj = obj2;
            sol = sol2;
            
        }
        string get_FP_lpfile() const {return lpfile;}
        char get_OVPR() const {return OVPR;}
        void set_OVPR(char S){OVPR = S; return;}
        
    };
    
    

#endif
