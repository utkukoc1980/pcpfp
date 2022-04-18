#include "parallel_cpfp_stuff.h"

#ifndef MY_SOLUTION_H
#define MY_SOLUTION_H


	class My_solution{
		int original_size;
		int size;
		vector<int> indices;
		vector<double> elements;
		double original_objective;
		double auxiliary_objective;
	public:
		vector<double>  get_elements() const;
		void set_elements(vector<double>  arg);
		vector<int> get_indices() const;
		void set_indices(vector<int>  arg);
		int get_size() const;
		void set_size(int arg);
		int get_original_size() const;
		void set_original_size(int arg);
		double get_original_objective() const;
		void set_original_objective(double  arg);
		double get_auxiliary_objective() const;
		void set_auxiliary_objective(double arg);
		
		double get_element(int index,int &status ) const;
        int set_element(int index, double  arg);
        
        int element_pp(int index);
        int element_pe(int index,double arg);
        int element_mm(int index);
        int element_me(int index,double arg);
        
		
		My_solution();
		My_solution(vector<double> arg);
		My_solution(int ncols, const double *arg);

		My_solution(const My_solution &rhs);

		My_solution operator=(const My_solution &other);
		void normalize(std::ostream &out= std::cout);

        vector<double>  get_solution_vector() const;


		void divide_coefficients(double x);
#ifdef USING_MPI
		void _send_mpi(int destination, int tag, vector<MPI_Request> &v_request, vector<MPI_Status> v_status, MPI_SEND_TYPE SEND_TYPE = MPI_ISEND_NONBLOCKING);
		void _recv_mpi(int source, int tag, vector<MPI_Request> &v_request, vector<MPI_Status> v_status, MPI_RECV_TYPE RECV_TYPE = MPI_IRECV_NONBLOCKING);
#else
		void _send_alt1(int destination, int tag, my_signaler comm, ALT1_SEND_TYPE = ALT1_SEND_NONBLOCKING);
		void _recv_alt1(int source, int tag, my_signaler comm, ALT1_RECV_TYPE RECV_TYPE = ALT1_RECV_NONBLOCKING);
#endif
        
        
		void print(std::ostream & out = std::cout);
        void shortlineprint( std::ostream & out = std::cout);
        void shortlineprint_selected_indices(vector<int> selected_indices, bool pos ,  std::ostream & out = std::cout);
        
        
//         _selected_indices
        
        bool operator==(const My_solution &rhs) const;
        
        bool equality_on_selected_indices(const My_solution &rhs, vector<int> selected_indices) const;
        
        void print_selected_indices(vector<int> selected_indices, std::ostream & out = std::cout) const;
        void print_selected_nonzero_indices(vector<int> selected_indices, std::ostream & out = std::cout) const;
        
        bool is_empty();
        bool operator<(const My_solution &rhs) const;

		void clear(); /** clears the content except original size*/
	};

#endif
// kate: indent-mode cstyle; indent-width 1; replace-tabs on;
