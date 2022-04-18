#include "parallel_cpfp_stuff.h"

#ifndef MY_OBJECTIVE_FUNCTION_H
#define MY_OBJECTIVE_FUNCTION_H
class My_objective_function{
		int original_size;
		int size;
		vector<int> indices;
		vector<double> elements;
	public:
		My_objective_function();
		My_objective_function(vector<double> arg);
		My_objective_function(const My_objective_function &other);

        vector<double>  get_elements() const;
        double  get_element(int whichindex) const;
		void set_elements(vector<double>  arg);

        int get_index(int whichindex) const;

        vector<int> get_indices() const;
		void set_indices(vector<int>  arg);

		int get_size() const;
		void set_size(int arg);

		int get_original_size() const;
		void set_original_size(int arg);

		int get_objective_function(vector<double> &retval, ostream &asd3 = std::cout)const;
		int get_objective_function(double *retval, ostream &asd3= std::cout)const;

        double calculate_norm();
        
		My_objective_function operator=(const My_objective_function &other);
        
        My_objective_function operator+=(const My_objective_function &other);
                
        
        
        

		void normalize(std::ostream &out= std::cout);


		void divide_coefficients(double x);

		int calculate_value(double &retval, const double *x);
 		int calculate_value_of_vector(double &retval, vector<double> &sol);
#ifdef USING_MPI
		void _send_mpi(int destination, int tag, vector<MPI_Request> &v_request, vector<MPI_Status> v_status, MPI_SEND_TYPE SEND_TYPE = MPI_ISEND_NONBLOCKING);
        void _recv_mpi(int source, int tag, vector<MPI_Request> &v_request, vector<MPI_Status> v_status, MPI_RECV_TYPE RECV_TYPE = MPI_IRECV_NONBLOCKING);
#else
        void _send_alt1(int destination, int tag, my_signaler comm, ALT1_SEND_TYPE SEND_TYPE = ALT1_SEND_NONBLOCKING);
        void _recv_alt1(int source,int tag, my_signaler comm,  ALT1_RECV_TYPE RECV_TYPE = ALT1_RECV_NONBLOCKING);
#endif
		void print(std::ostream & out = std::cout);
		void shortlineprint(std::ostream & out = std::cout);
		
        bool is_empty();
        
        
        
	};
    
    My_objective_function My_objective_function_multiply_sum(double l1, const My_objective_function &lhs, double r, const My_objective_function &rhs);
    
    
#endif
// kate: indent-mode cstyle; indent-width 1; replace-tabs on;
