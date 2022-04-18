#include "parallel_cpfp_stuff.h"
#include "my_objective_function.h"

#ifndef MYSERIALIZABLEROWCUT_H
#define MYSERIALIZABLEROWCUT_H
	class mySerializableRowCut{

		double lb;
		double ub;
		char sense;
// 		int original_size;

		int size;
		vector<int> indices;
		vector<double> elements;
		double norm;
        double efficiency[CUT_EFFICIENCY_END];

		CoinPackedVector packedvector;
		OsiRowCut rowcut;


	public:


		void internal_sort();

		mySerializableRowCut();
		mySerializableRowCut(const mySerializableRowCut &other);

		mySerializableRowCut(vector<int> ind, vector<double> ele, double llb, double uub);

		mySerializableRowCut(OsiRowCut* cut2);

		mySerializableRowCut(OsiRowCut cut);



		void print(std::ostream &out = std::cout, bool eff=true);

		CoinPackedVector generate_CoinPackedVector(std::ostream &out= std::cout);

		OsiRowCut* generate_OsiRowCutPointer();

		OsiRowCut generate_OsiRowCut(std::ostream &out= std::cout);

		void normalize(std::ostream &out= std::cout);
		void divide_coefficients(double x);
        void multiply_coefficients(double x);
        
		bool operator<(const mySerializableRowCut &rhs) const;
        bool operator==(const mySerializableRowCut &rhs) const;
        mySerializableRowCut &operator=(const mySerializableRowCut &other);

		mySerializableRowCut(double arg_lb,
							 double arg_ub,
							 char arg_sense,
					   int arg_size,
					   vector<double> arg_elements,
					   vector<int> arg_indices	);

		double get_lb() const ;
		void set_lb(double arg) ;

		double get_ub() const;
		void set_ub(double arg);

		double get_norm() const;
		void set_norm(double arg) ;

		char get_sense() const ;
		void set_sense(char arg) ;


		int get_size() const;
		void set_size(int arg) ;
	/*	int get_original_size() const;
		void set_original_size(int arg) ;
	*/
		vector<double>  get_elements() const;
		void set_elements(vector<double>  arg);

		vector<int> get_indices() const;
		void set_indices(vector<int>  arg);

		CoinPackedVector get_packedvector() const;
		void set_packedvector(CoinPackedVector arg);
		OsiRowCut get_rowcut() const;
		void set_rowcut(OsiRowCut arg);


		void set_efficiency(double arg[CUT_EFFICIENCY_END]);
        void set_efficiency(CUT_EFFICIENCY_ENM enm , double arg);

		void reset_efficiency();
		double* get_efficiency();
		double get_efficiency(CUT_EFFICIENCY_ENM enm) const ;
#ifdef USING_MPI
		void _send_mpi(int destination, int tag, vector<MPI_Request> &v_request, vector<MPI_Status> v_status, MPI_SEND_TYPE SEND_TYPE = MPI_ISEND_NONBLOCKING);
		void _recv_mpi(int source, int tag, vector<MPI_Request> &v_request, vector<MPI_Status> v_status, MPI_RECV_TYPE RECV_TYPE = MPI_IRECV_NONBLOCKING);
#else
        void _send_alt1(int destination, int tag, my_signaler comm, ALT1_SEND_TYPE SEND_TYPE = ALT1_SEND_NONBLOCKING);
        void _recv_alt1(int source,int tag, my_signaler comm,  ALT1_RECV_TYPE RECV_TYPE = ALT1_RECV_NONBLOCKING);
#endif
//         void calculate_efficiency(const vector<double> &sol_elements, const vector<int> &sol_indices);

        void calculate_efficiency(const double *x); /** assumes size x = original_size*/
        void calculate_efficiency_v(const vector<double> x); /** assumes size x = original_size*/

        void calculate_efficiency_wrt_objective(My_objective_function &obj,bool is_org); /** assumes objective_indices are sorted*/


        void update_efficiency_integral_support(const vector<int> &column_types, int n_binary_plus_n_int);
	};

#endif
// kate: indent-mode cstyle; indent-width 1; replace-tabs on;
