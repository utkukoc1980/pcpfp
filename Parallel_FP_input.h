
#ifndef Parallel_FP_input_H
    #define Parallel_FP_input_H
	#include "parallel_cpfp_stuff.h"
    using namespace std;
    
    class Parallel_FP_input{
        
        int maxIter1;
        int maxIter2;
        int maxit1wi;
        int maxit2wi;
        int minChange;
        double trsmin;
        int trsrnd;
        bool primal;
        bool barrier;
        bool autoalg;
        int s3me;
        bool presolve;
        bool mippresolve;
        double imp;
        double alpha_start;
        double alpha_quot;
        double alphadist; 
        string lpfile;
        char OVPR;
        

        
    public:
        int get_FP_maxIter1() const {return maxIter1;}
        int get_FP_maxIter2() const {return maxIter2;}
        int get_FP_maxit1wi() const {return maxit1wi;}
        int get_FP_maxit2wi() const {return maxit2wi;}
        int get_FP_minChange()const {return minChange;}
        double get_FP_trsmin() const {return trsmin;}
        int get_FP_trsrnd() const {return trsrnd;}
        bool get_FP_primal() const {return primal;}
        bool get_FP_barrier() const {return barrier;}
        bool get_FP_autoalg() const {return autoalg;}
        int get_FP_s3me() const {return s3me;}
        bool get_FP_presolve() const {return presolve;}
        bool get_FP_mippresolve() const {return mippresolve;}
        double get_FP_imp() const {return imp;}
        double get_FP_alpha_start() const {return alpha_start;}
        double get_FP_alpha_quot() const {return alpha_quot;}
        double get_FP_alphadist() const {return alphadist;}
        string get_FP_lpfile() const {return lpfile;}
        char get_FP_OVPR() const {return OVPR;}
        
    };
    
    

#endif
