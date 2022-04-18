#ifndef PARALLEL_CP
	#define PARALLEL_CP

// 	#include <boost/archive/binary_oarchive.hpp>
// 	#include <boost/archive/binary_iarchive.hpp>
// 	#include <boost/archive/text_oarchive.hpp>
// 	#include <boost/archive/text_iarchive.hpp>
// 	#include <boost/serialization/string.hpp>
// 	#include <boost/serialization/vector.hpp>
// 	#include <boost/serialization/list.hpp>
	#include "OsiCuts.hpp"

	#include "OsiSolverInterface.hpp"
	#include "OsiSolverParameters.hpp"
	#include "OsiClpSolverInterface.hpp"


// 	extern std::vector< char > OVPR;	//

	#include "CoinPackedVectorBase.hpp"

	#include "CglKnapsackCover.hpp"
	#include "CglSimpleRounding.hpp"
	#include "CglRedSplit.hpp"
	#include "CglGomory.hpp"


	#include "CglAllDifferent.hpp"
	#include "CglPreProcess.hpp"
	#include "CglLandP.hpp"


	#include "../templates.h"

	enum CUT_TYPE_ENM
	{
		CUT_TYPE_BGN                   =    0, /**< begin tag of the cut type */
		CUT_TYPE_GOMORY 				=    0, /**< use gomory cuts */
		CUT_TYPE_KNAPSACK				= 	 1, /**< use knapsack cuts */
		CUT_TYPE_REDSPLIT				= 	 2, /**< use reducsed split cuts */
		CUT_TYPE_SIMPLE_ROUNDING		= 	 3, /**<use simple rounding cuts */
		CUT_TYPE_END					=    4, /**<end tag for the cut type */
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
		OBJECTIVE_TYPE_END								= 7, /**<end tag for objective type*/
	};

	enum OBJECTIVE_PERTURBATION_ENUM{

	};

#define	_VERSION_1_FOR_CHANGING_OBJECTIVE
// #define USING_IOPTIMIZE

// 	#ifdef USING_IOPTIMIZE
// 		#define USE_NEW_IOPTIMIZE
// 		#include "../start_point_generator5.h"
// 	#endif

#define __CUT_EQUALITY_TOLERANCE 1e-6

#define USE_SET
	using namespace std;

	enum OBJECTIVE_SENSE {MAXIMIZE, MINIMIZE};

	const unsigned  my_default_seed = 5489;

	int sample100100(int n, vector<double> &vec, unsigned seed = my_default_seed){
		// initialize PRNG at first call
		static uniform_real_distribution<double> uni(-100,100);
		static mt19937 prng(seed); // Mersenne twister MT19937

		for (int i = 0; i < n; i++) vec[i] = uni(prng);
		return n;
	}

	int sampleAB(int n, vector<double>& vec, double l, double u, unsigned seed = my_default_seed){
		// initialize PRNG at first call
		vec.resize(n);
		static uniform_real_distribution<double> uni(l,u);
		static mt19937 prng(seed); // Mersenne twister MT19937

		for (int i = 0; i < n; i++) vec[i] = uni(prng);
		return n;
	}

	int sampleNN(int n, vector<double>& vec, int NN, unsigned seed = my_default_seed){
		return sampleAB(n,vec, -NN, NN, seed);

		// initialize PRNG at first call
		vec.resize(n);
		static uniform_real_distribution<double> uni(-NN,NN);
		static mt19937 prng(seed); // Mersenne twister MT19937

		for (int i = 0; i < n; i++) vec[i] = uni(prng);
		return n;
	}



	class my_objective_function{
		int size;
		vector<int> indices;

		vector<double> elements;

		public:
			friend class boost::serialization::access;
// 			template<class Archive>
//
// 			void serialize(Archive &ar, const unsigned int version){
//
// 					ofstream sout("load_temp2");
//
// 					sout <<__LINE__ << endl;
//
// 					ar & size;
// 					sout << "size " << size << endl;
//
// 					sout <<__LINE__ << endl;
//
// 					ar & indices;
// 					sout <<__LINE__ << endl;
// 					ar & elements;
// 					sout <<__LINE__ << endl;
//
//
// 				}
			template<class Archive>

			void load(Archive &ar, const unsigned int version){

				ofstream sout("load_temp2___");

				sout <<__LINE__ << endl;

				ar & size;
				sout << "size " << size << endl;

				sout <<__LINE__ << endl;

				ar & indices;
				sout <<__LINE__ << endl;
				ar & elements;
				sout <<__LINE__ << endl;


			}
			template<class Archive>
			void save(Archive &ar, const unsigned int version)const{
				ofstream sout("save_temp2");

				sout <<__LINE__ << endl;

				ar & size;
				sout <<__LINE__ << endl;
				sout << "size " << size << endl;

				for (int i = 0; i <size; ++i){
					sout << "x_" << indices[i] << ": " << elements[i] << endl;
				}
				ar & indices;
				sout <<__LINE__ << endl;
				ar & elements;
				sout <<__LINE__ << endl;


			}
			BOOST_SERIALIZATION_SPLIT_MEMBER();

			my_objective_function(){}
			my_objective_function(vector<double> arg){

				for (int i = 0; i < arg.size(); ++i){
					if (!double_equality(arg[i],0)){
						indices.push_back(i);
						elements.push_back(arg[i]);
					}
				}
				normalize();
				size = elements.size();
			}
			my_objective_function(const my_objective_function &other){
				size= other.size;
				elements = other.elements;
				indices = other.indices;

			}

// 			void set_function()
			vector<double>  get_elements() const { return elements;}
			void set_elements(vector<double>  arg) { elements = arg; return; }

			vector<int> get_indices() const { return indices;}
			void set_indices(vector<int>  arg) { indices = arg; return; }

			int get_size() const { return size;}
			void set_size(int arg) { size = arg; return; }

			int get_objective_function(vector<double> &retval, int size2, ostream &asd3 = std::cout)const{

				retval.resize(size2);
				asd3 << __LINE__ << endl;
				asd3 << "size "<< size << endl;
				asd3 << "size2 "<< size2 << endl;
				asd3 << "indices.size() "<< indices.size() << endl;
				for (int i = 0 ;i<size;++i){
					asd3 << __LINE__ << endl;
					asd3 << "indices["<<i<<"]: " << indices[i] << endl;

					if (indices[i]>=size2) return -1;
					asd3 << __LINE__ << endl;
					retval[indices[i]] = elements[i];
					asd3 << __LINE__ << endl;

				}
				asd3 << __LINE__ << endl;

				return 0;
			}
			my_objective_function operator=(const my_objective_function &other){
				size= other.size;
				elements = other.elements;
				indices = other.indices;
				return *this;

			}
			void normalize(std::ostream &out= std::cout){
				double norm =0;
				if(norm <  1.0e-20){
					norm =0;
					for (int i =0;i < elements.size();++i){
						norm += elements[i]*elements[i];
					}
				}
				divide_coefficients(sqrt(norm));

				/*double lhs = 0, rhs= 0, divisor = DBL_MAX;
				 *			if (lb > -DBL_MAX ) {
				 *				lhs = fabs(lb);
				 *				divisor = lhs;
				 }
				 if (ub < DBL_MAX){
					 rhs = fabs(ub);
					 if (rhs < divisor) divisor = rhs;
				 }



				 if (divisor < 1.0e-20) {
					 divisor = norm;
				 }

				 if(fabs (divisor - 1) > 1.0e-20) divide_coefficients(divisor);
				 */
			}

			void divide_coefficients(double x){
				if (x<0){
					x*=-1;
				}

				for (int i = 0; i < elements.size();++i){
					elements[i]/=x;
				}

				return;
			}


	};



	class mySerializableRowCut{
		double lb;
		double ub;
		char sense;
		int size;
		vector<double> elements;
		vector<int> indices;
		double norm;

		CoinPackedVector packedvector;
		OsiRowCut rowcut;


	public:
		friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive &ar, const unsigned int version){
			ar & lb;
			ar & ub;
			ar & sense;
			ar & size;
			ar & norm;

			ar & elements;
			ar & indices;
			//print();
		}

		void internal_sort(){
			vector<pair<int,double>> v;
			for (int i = 0; i < size; ++i){
				pair<int,double> asd(indices[i],elements[i]);
				v.push_back(asd);
			}
			sort(v.begin(),v.end());

// 			indices.assign();
			for (int i = 0; i < size; ++i){
				indices[i] = v[i].first;
				elements[i] = v[i].second;
			}
			generate_OsiRowCut();
			return;
		}

		mySerializableRowCut(){};


		mySerializableRowCut(const mySerializableRowCut &other){
			lb = other.lb;
			ub = other.ub;
			sense = other.sense;
			size = other.size;
			elements = other.elements;
			indices = other.indices;
			norm = other.norm;
		}

		mySerializableRowCut(vector<int> ind, vector<double> ele, double llb, double uub){
			lb = llb;
			ub = uub;
			size = ind.size();
			elements = ele ;
			indices = ind;
			sense = 'U';
			generate_OsiRowCut();
		}

		mySerializableRowCut(OsiRowCut* cut2){

			lb = cut2->lb();
			ub = cut2->ub();
			sense = cut2->sense();

			packedvector = CoinPackedVector(cut2->mutableRow());
// 			rowcut = OsiRowCut(cut2);

// 			OsiRowCut cut(cut2);
			//comment out these two
			// 			OsiRowCut rowcut;
			// 			CoinPackedVector packedvector;
// 			OsiRowCut dummy = new OsiRowCut(cut2);
// 			cout << cut2 << endl;
// 			rowcut = dummy;
// 			rowcut = new OsiRowCut(cut2);



// 			lb = rowcut.lb();
// 			ub = rowcut.ub();
// 			sense = rowcut.sense();
//
//
// 			packedvector = rowcut.row();

			size = packedvector.getNumElements();

			int *indices2; // = new int[size];
			double *elements2;// = new double[size];

			indices2  = packedvector.getIndices();
			elements2 = packedvector.getElements();

			elements.resize(size);
			indices.resize(size);

			for (int i = 0; i<size; ++i){
				elements[i] = elements2[i];
				indices[i] = indices2[i];
			}
			norm = 0;
			internal_sort();
			normalize();


			// 			delete[] indices2;
			// 			delete[] elements2;
// 			rowcut.setRow(packedvector);
// 			rowcut.setLb(lb);
// 			rowcut.setUb(ub);

		}

		mySerializableRowCut(OsiRowCut cut){
			//comment out these two
// 			OsiRowCut rowcut;
// 			CoinPackedVector packedvector;

			rowcut = OsiRowCut(cut);
			lb = rowcut.lb();
			ub = rowcut.ub();
			sense = rowcut.sense();


			packedvector = rowcut.mutableRow();

			size = packedvector.getNumElements();

			int *indices2; // = new int[size];
			double *elements2;// = new double[size];

			indices2  = packedvector.getIndices();
			elements2 = packedvector.getElements();

			elements.resize(size);
			indices.resize(size);

			for (int i = 0; i<size; ++i){
				elements[i] = elements2[i];
				indices[i] = indices2[i];
			}
			norm = 0;
			internal_sort();
			normalize();
// 			delete[] indices2;
// 			delete[] elements2;

		}



		void print(std::ostream &out = std::cout){
// 		    out << "printing mySerializableRowCut"  << endl;
			if (size < 1) {
				out << "size < 1 nothing to print " << endl;

				out << "elements.size() " << elements.size() << endl;
				out << "indices.size() " << indices.size() << endl;
				out << "lb "<< lb << endl;
				out << "ub "<< ub << endl;
				return;
			}

			if (lb > -DBL_MAX)
				out << lb << " <= ";

			for (int i = 0; i<size; ++i){
				if (elements[i] > 0 )
					out << "+" << elements[i] << " * x_" << indices[i] << " " ;
				else
					out << elements[i] << " * x_" << indices[i] << " " ;
			}

			if (ub < DBL_MAX)
				out << " <= " << ub << " Eff: " << rowcut.effectiveness() << "\t NORM : "<< norm << endl;
			else
				out << " Eff: " << rowcut.effectiveness() << "\t NORM : "<< norm << endl;
			return;


			if (size < 1) return;
			out << lb << " <= ";
			for (int i = 0; i<size; ++i){
				if (elements[i] > 0 )
					out << "+" << elements[i] << " * x_" << indices[i] << " " ;
				else
					out << elements[i] << " * x_" << indices[i] << " " ;


			}
			out << " <= " << ub  << endl;
			return;
		}

		CoinPackedVector generate_CoinPackedVector(std::ostream &out= std::cout){

			//comment out these two
// 			OsiRowCut rowcut;
// 			CoinPackedVector packedvector;


			int *indices2 = new int[size];
			double *elements2 = new double[size];
			for (int i = 0; i<size; ++i){
				elements2[i] = elements[i];
				indices2[i] = indices[i];
			}
// 			out << "before calling CoinPackedVector(size,indices2,elements2);" << endl;
			packedvector = CoinPackedVector(size,indices2,elements2);
// 			out << "after calling CoinPackedVector(size,indices2,elements2);" << endl;

			delete[] indices2;
			delete[] elements2;

// 			out << "deleted indices2 and elements2  returning " << endl;
			return packedvector;

		}

		OsiRowCut* generate_OsiRowCutPointer(){
			//comment out these two
			// 			OsiRowCut rowcut;
			// 			CoinPackedVector packedvector;

			rowcut.setRow(generate_CoinPackedVector());

			/* make this more complicated wrt sense */
			rowcut.setLb(lb);
			rowcut.setUb(ub);
			return &rowcut;
		}

		OsiRowCut generate_OsiRowCut(std::ostream &out= std::cout){
			//comment out these two
// 			OsiRowCut rowcut;
// 			CoinPackedVector packedvector;
			out << "setRow "<< endl;
			rowcut.setRow(generate_CoinPackedVector(out));
			out << "setLB "<< endl;

			/* make this more complicated wrt sense */
			rowcut.setLb(lb);
			out << "setUB "<< endl;

			rowcut.setUb(ub);

			out << "return "<< endl;
			out << "rowcut.print() "<< endl;

			rowcut.print();
			out << "rowcut.print() end"<< endl;

            out << "this->print(out)" << endl;
            this->print();
			return rowcut;
		}

		void normalize(std::ostream &out= std::cout){
			if(norm <  1.0e-20){
				norm =0;
				for (int i =0;i < elements.size();++i){
					norm += elements[i]*elements[i];
				}
			}
			divide_coefficients(sqrt(norm));

			/*double lhs = 0, rhs= 0, divisor = DBL_MAX;
			if (lb > -DBL_MAX ) {
				lhs = fabs(lb);
				divisor = lhs;
			}
			if (ub < DBL_MAX){
				rhs = fabs(ub);
				if (rhs < divisor) divisor = rhs;
			}



			if (divisor < 1.0e-20) {
				divisor = norm;
			}

			if(fabs (divisor - 1) > 1.0e-20) divide_coefficients(divisor);
			*/
		}

		void divide_coefficients(double x){
			if (x<0){
				x*=-1;
			}

			for (int i = 0; i < elements.size();++i){
				elements[i]/=x;
			}

			if (lb > -DBL_MAX ) {
				lb/=x;
			}

			if (ub < DBL_MAX){
				ub/=x;
			}
			norm /= (x*x);

			return;
		}

		bool operator<(const mySerializableRowCut &rhs) const {


			if (size != rhs.size) return (size < rhs.size);


// 			if (sense != rhs.sense) return (sense < rhs.sense);


			if (fabs(lb - rhs.lb) > 1e-10) return (lb < rhs.lb);
			if (fabs(ub - rhs.ub) > 1e-10) return (ub < rhs.ub);
			for (int i = 0; i < size; ++i){
				if (indices[i] != rhs.indices[i]) return (indices[i] < rhs.indices[i]);
				if (fabs(elements[i] - rhs.elements[i]) > 1e-10) return (elements[i] < rhs.elements[i]);
			}

			if (fabs(norm - rhs.norm) < 1e-20) return (norm<rhs.norm);

			return false;

		}

			// NEED TO CHECK IF LBs or UBs = inf and sense is not an issue;
		bool operator==(const mySerializableRowCut &rhs) const{

			if (size != rhs.size) return (false);
			if (fabs (norm - rhs.norm) > __CUT_EQUALITY_TOLERANCE) return (false);
// 			if (sense != rhs.sense) return (false);
			if (fabs(lb - rhs.lb) > __CUT_EQUALITY_TOLERANCE) return (false);
			if (fabs(ub - rhs.ub) > __CUT_EQUALITY_TOLERANCE) return (false);
			for (int i = 0; i < size; ++i){
				if (indices[i] != rhs.indices[i]) return (false);
				if (fabs(elements[i] - rhs.elements[i]) > __CUT_EQUALITY_TOLERANCE) return (false);
			}

			return true;
		}

		mySerializableRowCut(double arg_lb,
							 double arg_ub,
							 char arg_sense,
					   int arg_size,
					   vector<double> arg_elements,
					   vector<int> arg_indices	):
						lb(arg_lb),
						ub(arg_ub),
						sense(arg_sense),
						size(arg_size),
						elements(arg_elements),
						indices(arg_indices){}

		double get_lb() const { return lb;}
		void set_lb(double arg) { lb = arg; return; }

		double get_ub() const { return ub;}
		void set_ub(double arg) { ub = arg; return; }

		double get_norm() const { return norm;}
		void set_norm(double arg) { norm = arg; return; }

		char get_sense() const {return sense;}
		void set_sense(char arg) {sense = arg; return;}

		int get_size() const { return size;}
		void set_size(int arg) { size = arg; return; }

		vector<double>  get_elements() const { return elements;}
		void set_elements(vector<double>  arg) { elements = arg; return; }

		vector<int> get_indices() const { return indices;}
		void set_indices(vector<int>  arg) { indices = arg; return; }

		CoinPackedVector get_packedvector() const { return packedvector;}
		void set_packedvector(CoinPackedVector arg) { packedvector = arg; return; }

		OsiRowCut get_rowcut() const { return rowcut;}
		void set_rowcut(OsiRowCut arg) { rowcut = arg; return; }



	};


	vector<mySerializableRowCut> rowCutSelectionProcess(vector<mySerializableRowCut> initial_list, int size = 1){
		vector<mySerializableRowCut> final_list;

		return initial_list;
	}

	class mySerializableSolution{
		double original_objective;
		double auxilary_objective;
		vector<double> solution;

	public:
		friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive &ar, const unsigned int version){
			ar & original_objective;
			ar & auxilary_objective;
			ar & solution;

		}

		bool operator==(const mySerializableSolution &rhs){
			if(solution.size() != rhs.solution.size()) return false;
			if (fabs(original_objective - rhs.original_objective) > 1e-10) return false;
			for (int i = 0; i<solution.size(); ++i){
				if (fabs(solution[i] - rhs.solution[i]) > 1e-10) return false;
			}
			return true;
		}
		bool operator<(const mySerializableSolution &rhs){
			if(solution.size() != rhs.solution.size()) return false;
			for (int i = 0; i<solution.size(); ++i){
				if (solution[i] < (rhs.solution[i] - 1e-10)) return true;
				if (solution[i] > (rhs.solution[i] + 1e-10)) return false;
			}
			if (original_objective < (rhs.original_objective - 1e-10 )) return false;
			if (original_objective > (rhs.original_objective + 1e-10 )) return false;

			return true;
		}

		mySerializableSolution(){}

		mySerializableSolution(int n, const double * sol){
			solution.resize(n);
			for (int i = 0;i <n ;++i) solution[i] = sol[i];
		}

		mySerializableSolution(double arg_original_objective,
							   double arg_auxilary_objective,
						 vector<double> arg_solution):original_objective(arg_original_objective),
						 auxilary_objective(arg_auxilary_objective),
						 solution(arg_solution){}




		mySerializableSolution(const mySerializableSolution &other){
			original_objective = other.original_objective;
			auxilary_objective = other.auxilary_objective;
			solution = other.solution;
		}
		mySerializableSolution operator=(const mySerializableSolution &other){
			original_objective = other.original_objective;
			auxilary_objective = other.auxilary_objective;
			solution = other.solution;
			return *this;
		}

		double get_original_objective() const { return original_objective;}
		void set_original_objective(double arg) { original_objective = arg; return; }

		double get_auxilary_objective() const { return auxilary_objective;}
		void set_auxilary_objective(double arg) { auxilary_objective = arg; return; }

		vector<double> get_solution() const { return solution;}
		void set_solution(vector<double> arg) { solution = arg; return; }

		void print(std::ostream & out = std::cout){
			out << "org_obj : "<< original_objective << endl;
			out << "aux obj : " << auxilary_objective << endl;
			for (int i = 0;i <solution.size() ;++i)
				out << solution[i] << "\t";

			out << endl;
		}


	};

	class Pcp_output{
		bool is_cut;
		bool is_incumbent;
		bool is_feasible;
		vector<mySerializableRowCut> rowcutlist;
		vector<mySerializableSolution> solutionlist;

		int serialization_size_rowcutlist;
		int serialization_size_solutionlist;


	public:
		friend class boost::serialization::access;
// 		template<class Archive>
// 		void serialize(Archive &ar, const unsigned int version){
// // 			ar & input_type;
// 			ar & is_cut;
// 			ar & is_incumbent;
// 			ar & serialization_size_rowcutlist;
// 			for (int i = 0; i < serialization_size_rowcutlist; ++i){
// 				ar & rowcutlist[i];
// 			}
// 			ar & serialization_size_solutionlist;
// 			for (int i = 0; i < serialization_size_solutionlist; ++i){
// 				ar & solutionlist[i];
// 			}
// 		}

		template<class Archive>
		void save(Archive &ar, const unsigned int version) const{
			ar & is_cut;
			ar & is_incumbent;
			ar & is_feasible;
			if(is_cut){
			ar & serialization_size_rowcutlist;
				for (int i = 0; i < serialization_size_rowcutlist; ++i){
					ar & rowcutlist[i];
				}
			}
			if (is_incumbent){
				ar & serialization_size_solutionlist;
				for (int i = 0; i < serialization_size_solutionlist; ++i){
					ar & solutionlist[i];
				}
			}
		}
		template<class Archive>
		void load(Archive &ar, const unsigned int version){
			ar & is_cut;
			ar & is_incumbent;
			ar & is_feasible;
			if(is_cut){
			ar & serialization_size_rowcutlist;
				rowcutlist.resize(serialization_size_rowcutlist);
				for (int i = 0; i < serialization_size_rowcutlist; ++i){
					ar & rowcutlist[i];
				}
			}
			if (is_incumbent){
				ar & serialization_size_solutionlist;
				solutionlist.resize(serialization_size_solutionlist);
				for (int i = 0; i < serialization_size_solutionlist; ++i){
					ar & solutionlist[i];
				}
			}
		}
		BOOST_SERIALIZATION_SPLIT_MEMBER();

		bool operator< (const Pcp_output &rhs) const{
			return false;
		}

		Pcp_output(){}
		Pcp_output(int i):is_cut(false),is_incumbent(false),is_feasible(true){}

// 		Pcp_output(int arg_input_type):input_type(arg_input_type){}

		Pcp_output(vector<mySerializableSolution> arg_solutionlist):is_cut(false),
			is_incumbent(true),
			is_feasible(true),
			solutionlist(arg_solutionlist),
			serialization_size_rowcutlist(0),
			serialization_size_solutionlist(0){}

		Pcp_output(vector<mySerializableRowCut> arg_rowcutlist):is_cut(true),
			is_incumbent(false),
			is_feasible(true),
			rowcutlist(arg_rowcutlist){}

		Pcp_output(bool arg_is_cut,
		bool arg_is_incumbent,
		vector<mySerializableRowCut> arg_rowcutlist,
		vector<mySerializableSolution> arg_solutionlist):is_cut(arg_is_cut),
		is_incumbent(arg_is_incumbent),
		is_feasible(true),
		rowcutlist(arg_rowcutlist),
		solutionlist(arg_solutionlist),
		serialization_size_rowcutlist(0),
		serialization_size_solutionlist(0){}

		Pcp_output(const Pcp_output &other){
			is_cut = other.is_cut;
			is_incumbent = other.is_incumbent;
			is_feasible = other.is_feasible;

			for (int i = 0; i< other.rowcutlist.size();++i){
				rowcutlist.push_back(other.rowcutlist[i]);
			}

			for (int i = 0; i< other.solutionlist.size();++i){
				solutionlist.push_back(other.solutionlist[i]);
			}

			serialization_size_rowcutlist = other.serialization_size_rowcutlist;
			serialization_size_solutionlist = other.serialization_size_solutionlist;
// 			rowcutlist = other.rowcutlist;
// 			solutionlist = other.solutionlist;

		}

		Pcp_output operator=(const Pcp_output &other){
			is_cut = other.is_cut;
			is_incumbent = other.is_incumbent;
			is_feasible = other.is_feasible;


			for (int i = 0; i< other.rowcutlist.size();++i){
				rowcutlist.push_back(other.rowcutlist[i]);
			}

			for (int i = 0; i< other.solutionlist.size();++i){
				solutionlist.push_back(other.solutionlist[i]);
			}
			serialization_size_rowcutlist = other.serialization_size_rowcutlist;
			serialization_size_solutionlist = other.serialization_size_solutionlist;

			return *this;
		}

		void add(mySerializableRowCut c ){
			is_cut = true;
			is_feasible = true;
			rowcutlist.push_back(c);
			serialization_size_rowcutlist = rowcutlist.size();
			return;
		}

		void add(mySerializableSolution s){
			is_incumbent = true;
			is_feasible = true;

			solutionlist.push_back(s);
			serialization_size_solutionlist = solutionlist.size();
			return;
		}

		void add(vector<mySerializableRowCut> c){
			is_cut = true;
			is_feasible = true;

			rowcutlist.insert(rowcutlist.end(),c.begin(),c.end());
			serialization_size_rowcutlist = rowcutlist.size();
			return;
		}

		void add(set<mySerializableRowCut> c){
			is_cut = true;
			rowcutlist.insert(rowcutlist.end(),c.begin(),c.end());
			serialization_size_rowcutlist = rowcutlist.size();
			return;
		}

		void add(vector<mySerializableSolution> s){
			is_incumbent = true;
			is_feasible = true;

			solutionlist.insert(solutionlist.end(),s.begin(),s.end());
			serialization_size_solutionlist = solutionlist.size();

			return;
		}



		bool get_is_cut() const { return is_cut;}
		void set_is_cut(bool arg) { is_cut = arg; return; }

		bool get_is_incumbent() const { return is_incumbent;}
		void set_is_incumbent(bool arg) { is_incumbent = arg; return; }

		bool get_is_feasible() const { return is_feasible;}
		void set_is_feasible(bool arg) { is_feasible = arg; return; }

		vector<mySerializableRowCut> get_rowcutlist() const { return rowcutlist;}
		void set_rowcutlist(vector<mySerializableRowCut> arg) { rowcutlist = arg; return; }

		vector<mySerializableSolution> get_solutionlist() const { return solutionlist;}
		void set_solutionlist(vector<mySerializableSolution> arg) { solutionlist = arg; return; }

		void print (std::ostream &out = std::cout){
			out << "is_feasible: " << is_feasible<< endl;
			out <<  "is_cut: "<< is_cut;
			if (is_cut){
				out << " serialization_size_rowcutlist: " << serialization_size_rowcutlist<< endl;
				for (int i = 0; i < serialization_size_rowcutlist; ++i){
					rowcutlist[i].print(out);
				}

			}
			else
				out << endl;

			out << "is_incumbent: " << is_incumbent;
// 			vector<mySerializableRowCut> rowcutlist;
// 			vector<mySerializableSolution> solutionlist;

			if (is_incumbent) {
				out << " serialization_size_solutionlist: " <<  serialization_size_solutionlist<< endl;
				for (int i = 0; i < serialization_size_solutionlist; ++i){
					solutionlist[i].print(out);
				}
			}
			else
				out << endl;
			return;

		}

		int get_serialization_size_rowcutlist() const { return serialization_size_rowcutlist;}
		void set_serialization_size_rowcutlist(int arg) { serialization_size_rowcutlist = arg; return;}

		int get_serialization_size_solutionlist() const { return serialization_size_solutionlist;}
		void set_serialization_size_solutionlist(int arg) { serialization_size_solutionlist = arg; return;}


	};

	class Cut_type{
		bool apply[CUT_TYPE_END];
	public:
		friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive &ar, const unsigned int version){
			ar & apply;
		}

		Cut_type(){}

		Cut_type(bool arg_apply[CUT_TYPE_END]){
			for (int i = CUT_TYPE_BGN; i < CUT_TYPE_END; ++i)
				apply[i] = arg_apply[i];
		}

		Cut_type(const Cut_type &other){
			for (int i = CUT_TYPE_BGN; i < CUT_TYPE_END; ++i)
				apply[i] = other.apply[i];
		}


		Cut_type &operator=(const Cut_type &other){
			for (int i = CUT_TYPE_BGN; i < CUT_TYPE_END; ++i)
				apply[i] = other.apply[i];
			return *this;
		}

		Cut_type(int i){
			for(int i = CUT_TYPE_BGN; i < CUT_TYPE_END;++i){
				apply[i] = false;
			}
			switch (i){
				case 0 : break;
				case 1 :
					apply[CUT_TYPE_GOMORY] = true;
					apply[CUT_TYPE_KNAPSACK] = true;
					apply[CUT_TYPE_REDSPLIT] = true;
					apply[CUT_TYPE_SIMPLE_ROUNDING] = true;
					break;
				default:
					apply[CUT_TYPE_GOMORY] = true;
					break;
			}

		}

		bool get_apply(int i) const{
			if ((i < CUT_TYPE_BGN) || (i>=CUT_TYPE_END) )
				return false;

			return apply[i];
		}

		void set_apply(int i, bool arg){
			if ((i < CUT_TYPE_BGN) || (i>=CUT_TYPE_END) )
				return;
			apply[i] = arg;
			return;
		}

		bool get_apply_gomory() const { return apply[CUT_TYPE_GOMORY];}
		void set_apply_gomory(bool arg) { apply[CUT_TYPE_GOMORY] = arg; return; }
		bool get_apply_knapsack() const { return apply[CUT_TYPE_KNAPSACK];}
		void set_apply_knapsack(bool arg) { apply[CUT_TYPE_KNAPSACK]= arg; return; }
		bool get_apply_redsplit() const { return apply[CUT_TYPE_REDSPLIT];}
		void set_apply_redsplit(bool arg) { apply[CUT_TYPE_REDSPLIT]= arg; return; }
		bool get_apply_simple_rouding() const { return apply[CUT_TYPE_SIMPLE_ROUNDING];}
		void set_apply_simple_rouding(bool arg) { apply[CUT_TYPE_SIMPLE_ROUNDING]= arg; return; }


	};


	class Pcp_input{
		OBJECTIVE_TYPE_ENUM obj_type;
		char OVPR;
		int input_type;
		double new_rhs_for_problem;
		double time_limit;
		vector<mySerializableRowCut> rowcutlist;

		Cut_type cut_type;

		int change_obj;
		int perturbation_type;
		int perturbation_level;
		double perturbation_percentage;
		bool apply_cuts;

		int serialization_size;
// 		int objective_function_size;

// 		vector<double> objective_function;
		my_objective_function func;
		bool dummy;

	public:

		friend class boost::serialization::access;

// 		template<class Archive>
// 		void serialize(Archive &ar, const unsigned int version){
// 			print();
// 			ar & OVPR;
// 			ar & input_type;
// 			ar & new_rhs_for_problem;
// 			ar & time_limit;
// 			ar & serialization_size;
// 			ar & cut_type;
// 			ar & change_obj;
// 			ar & perturbation_type;
// 			ar & perturbation_level;
// 			ar & apply_cuts;
//
// 			ar & objective_function_size;
// 			ar & perturbation_percentage;
//
// 			if (objective_function_size>0){
// 				ar & objective_function;
// 			}
// 			if (serialization_size>0){
// 				ar & rowcutlist;
// 			}
//
// 		}
// 		template<class Archive>
// 		void serialize(Archive &ar, const unsigned int version){
// 			ar & input_type;
// 			ar & new_rhs_for_problem;
// 			ar & time_limit;
// 			ar & serialization_size;
// // 			cout << "serialization_size CDASDSDF GSDFFS DSDF SDF SD SF: "<< serialization_size << endl;
// 			for (int i =0; i < serialization_size; ++i)
// 				ar & rowcutlist[i];
//
// 			ar & cut_type;
// 			ar & change_obj;
// 			ar & perturbation_type;
// 			ar & perturbation_level;
// 			ar & apply_cuts;
// 		}
		template<class Archive>
		void save(Archive &ar, const unsigned int version)const {
//  			print();
			ar & func;

			ar & obj_type;
			ar & OVPR;
			ar & input_type;
			ar & new_rhs_for_problem;
			ar & time_limit;
// 			serialization_size = rowcutlist.size();
			ar & serialization_size;
// 			cout <<__LINE__ << endl;
// 			cout << "serialization_size " << serialization_size<< endl;
// 			cout <<__LINE__ << endl;
			ar & cut_type;
// 			cout <<__LINE__ << endl;
			ar & change_obj;

// 			cout <<__LINE__ << endl;
// 			ar & objective_function_size;
// 			cout << "objective_function_size: " << objective_function_size <<endl;
// 			cout <<__LINE__ << endl;
			ar & perturbation_type;
// 			cout <<__LINE__ << endl;
			ar & perturbation_level;
// 			cout <<__LINE__ << endl;

			ar & perturbation_percentage;
// 			cout <<__LINE__ << endl;
			ar & apply_cuts;
// 			cout <<__LINE__ << endl;
//  			objective_function_size = objective_function.size();


			if (serialization_size > 0){
				for (int i =0; i < serialization_size; ++i){
					ar & rowcutlist[i];
				}
			}
// 			cout <<__LINE__ << endl;


// 			if(obj_type != OBJECTIVE_TYPE_ORIGINAL){
// 				for(int i = 0; i < objective_function_size; ++i ){
// 						cout <<"i= " << i << " line"<< __LINE__ << endl;
// 						cout << "objective_function[i]" << objective_function[i]  << endl;
// 						ar & objective_function[i];
//
// 					}
// 			}
// 			ar & objective_function;
// 			cout << "objective_function_size:  "<< objective_function_size  <<endl;
// 			if(objective_function_size > 0){
// 				cout <<__LINE__ << endl;
//
// 				for(int i = 0; i < objective_function_size; ++i ){
// 					cout <<"i= " << i << " line"<< __LINE__ << endl;
// 					cout << "objective_function[i]" << objective_function[i]  << endl;
// 					ar & objective_function[i];
//
// 				}
// 			}
// 			cout << "objective_function_size:  "<< objective_function_size  <<endl;
// 			cout <<__LINE__ << endl;

// 			ar & dummy;
// 			ar & dummy;
// 			ar & func;
			cout <<__LINE__ << endl;

		}


		template<class Archive>
		void load(Archive &ar, const unsigned int version){

// 			ofstream sout("load_temp");
// 			print(sout);
// 			sout << __LINE__ <<endl;

			ar & func;


			ar & obj_type;
			ar & OVPR;
// 			sout << __LINE__ <<endl;
			ar & input_type;
// 			sout << __LINE__ <<endl;
			ar & new_rhs_for_problem;
// 			sout << __LINE__ <<endl;
			ar & time_limit;
// 			sout << __LINE__ <<endl;
// 			sout << "serialization_size 1:" << serialization_size<< endl;

			ar & serialization_size;
// 			sout << "serialization_size 2: " << serialization_size<< endl;
// 			sout << __LINE__ <<endl;
			ar & cut_type;
// 			sout << __LINE__ <<endl;
			ar & change_obj;
// 			sout << __LINE__ <<endl;
// 			ar & objective_function_size;
// 			sout << "objective_function_size 2: " << objective_function_size <<endl;
// 			sout.flush();
// 			sout << __LINE__ <<endl;
			ar & perturbation_type;
// 			sout << __LINE__ <<endl;
			ar & perturbation_level;
// 			sout << __LINE__ <<endl;
			ar & perturbation_percentage;
// 			sout << __LINE__ <<endl;
			ar & apply_cuts;
// 			sout << __LINE__ <<endl;
// 			sout << "objective_function_size 1: " << objective_function_size <<endl;
// 			sout.flush();
// 			sout << __LINE__ <<endl;

// 			sout <<__LINE__ << endl;

// 			sout <<__LINE__ << endl;
// 			rowcutlist.clear();
// 			sout <<__LINE__ << endl;

			if (serialization_size > 0)	{
// 				sout <<__LINE__ << endl;
				rowcutlist.resize(serialization_size);
// 				sout <<__LINE__ << endl;
				for (int i =0; i < serialization_size; ++i){
// 					sout <<__LINE__ << endl;
					ar & rowcutlist[i];
				}
// 				sout <<__LINE__ << endl;

			}

// 			sout <<__LINE__ << endl;

// 				objective_function = func.get_elements();
// 				sout <<__LINE__ << endl;

// 				objective_function_size = objective_function.size();
// 				sout <<__LINE__ << endl;
// 				sout << "obj_type: " << obj_type << endl;;
// 				sout << "OBJECTIVE_TYPE_ORIGINAL: " << OBJECTIVE_TYPE_ORIGINAL << endl;;
// 			if((obj_type != OBJECTIVE_TYPE_ORIGINAL)){
// 					for(int i = 0; i < objective_function_size; ++i ){
// 						sout <<"i= " << i << " line"<< __LINE__ << endl;
// 						sout << "objective_function[i]" << objective_function[i]  << endl;
// 						ar & objective_function[i];
// 						sout <<"i= " << i << " line"<< __LINE__ << endl;
// 						sout << "objective_function[i]" << objective_function[i]  << endl;
//
// 					}
// 				}
//  			objective_function.clear();
// 			ar & objective_function;
// 			sout << "objective_function_size:  "<< objective_function_size  <<endl;
// 			sout <<__LINE__ << endl;
// 			if(objective_function_size > 0){
// 				sout <<__LINE__ << endl;
//
// 				objective_function.resize(objective_function_size);
// 				sout <<__LINE__ << endl;
// 				sout << "objective_function_size:  "<< objective_function_size  <<endl;
// 				sout << "objective_function.size:  "<< objective_function.size()  <<endl;
//
// 				for(int i = 0; i < objective_function_size; ++i ){
// 					sout <<"i= " << i << " line"<< __LINE__ << endl;
// 					double x =1;
// 					sout << "x " << x<< endl;
// 					sout << "objective_function[i]" << objective_function[i]<< endl;
//
// 					ar & x  ;
// 					sout << "x " << x<< endl;
//
// 					objective_function[i] = x;
// 					sout << "objective_function[i]" << objective_function[i]<< endl;
// 					sout << "x" << x<< endl;
//
// 				}
// 				sout <<__LINE__ << endl;
//
// 			}
// 			sout << "objective_function_size:  "<< objective_function_size  <<endl;
//
// 			sout <<__LINE__ << endl;
//
// 			ar & dummy;
// 			ar & dummy;
// 				ar & func;

// 				sout <<__LINE__ << endl;

		}


		BOOST_SERIALIZATION_SPLIT_MEMBER();

		Pcp_input(){
			obj_type = OBJECTIVE_TYPE_ORIGINAL;
			perturbation_type= 0;
// 			objective_function_size = objective_function.size();
			serialization_size = rowcutlist.size();

		}

        Pcp_input(int type):
			obj_type(OBJECTIVE_TYPE_ORIGINAL),
			OVPR('X'),input_type(type),new_rhs_for_problem(1e40),
			time_limit(0),
			change_obj(0),
			perturbation_type(0),
			perturbation_level(0),
			perturbation_percentage(0),
			apply_cuts(false){
			rowcutlist.clear();
			cut_type = Cut_type(0);
			serialization_size = rowcutlist.size();
// 			objective_function.clear();
// 			objective_function_size = objective_function.size();
			}

		Pcp_input(char arg_OVPR,
				  int arg_input_type,
				  double arg_new_rhs_for_problem,
			double arg_time_limit,
			vector<mySerializableRowCut> arg_rowcutlist,
				Cut_type arg_cut_type,

			int arg_change_obj,
			int arg_perturbation_type,
			int arg_perturbation_level,
			double arg_perturbation_percentage,
			bool arg_apply_cuts	):obj_type(OBJECTIVE_TYPE_ORIGINAL),OVPR(arg_OVPR),input_type(arg_input_type),
			new_rhs_for_problem(arg_new_rhs_for_problem),
			time_limit(arg_time_limit),
			rowcutlist(arg_rowcutlist),
			cut_type(arg_cut_type),
			change_obj(arg_change_obj),
			perturbation_type(arg_perturbation_type),
			perturbation_level(arg_perturbation_level),
			perturbation_percentage(arg_perturbation_percentage),
			apply_cuts(arg_apply_cuts){

				serialization_size = rowcutlist.size();
// 				objective_function.clear();
// 				objective_function_size = objective_function.size();
			}

			Pcp_input(const Pcp_input &other){
				obj_type = other.obj_type;

				OVPR = other.OVPR;
				input_type = other.input_type;
				new_rhs_for_problem = other.new_rhs_for_problem;
				time_limit = other.time_limit;
				rowcutlist = other.rowcutlist;
				cut_type = other.cut_type;
				change_obj = other.change_obj;
				perturbation_type = other.perturbation_type;
				perturbation_level = other.perturbation_level;
				perturbation_percentage = other.perturbation_percentage;
				apply_cuts = other.apply_cuts;
				serialization_size = other.serialization_size;
				func = other.func;
// 				objective_function = other.objective_function;
// 				objective_function_size = other.objective_function_size;

			}

			Pcp_input &operator=(const Pcp_input &other){
				obj_type = other.obj_type;
				OVPR = other.OVPR;
				input_type = other.input_type;
				new_rhs_for_problem = other.new_rhs_for_problem;
				time_limit = other.time_limit;
				rowcutlist = other.rowcutlist;
				cut_type = other.cut_type;
				change_obj = other.change_obj;
				perturbation_type = other.perturbation_type;
				perturbation_level = other.perturbation_level;
				perturbation_percentage = other.perturbation_percentage;
				func = other.func;

				apply_cuts = other.apply_cuts;
				serialization_size = other.serialization_size;


// 				objective_function = other.objective_function;
// 				objective_function_size = other.objective_function_size;
				func = other.func;

				return *this;
			}

		OBJECTIVE_TYPE_ENUM get_obj_type() const { return obj_type;}
		void set_obj_type(OBJECTIVE_TYPE_ENUM arg) { obj_type = arg; return; }

		char get_OVPR() const { return OVPR;}
		void set_OVPR(char arg) { OVPR = arg; return; }

		int get_input_type() const { return input_type;}
		void set_input_type(int arg) { input_type = arg; return; }

		double get_new_rhs_for_problem() const { return new_rhs_for_problem;}
		void set_new_rhs_for_problem(double arg) { new_rhs_for_problem = arg; return; }

		double get_time_limit() const { return time_limit;}
		void set_time_limit(double arg) { time_limit = arg; return; }

		vector<mySerializableRowCut> get_rowcutlist() const { return rowcutlist;}
		void set_rowcutlist(vector<mySerializableRowCut> arg) { serialization_size = arg.size();
		rowcutlist = arg; return; }


		Cut_type get_cut_type() const { return cut_type;}
		void set_cut_type(Cut_type arg) { cut_type = arg; return; }

		int get_change_obj() const { return change_obj;}
		void set_change_obj(int arg) { change_obj = arg; return; }

		int get_perturbation_type() const { return perturbation_type;}
		void set_perturbation_type(int arg) { perturbation_type = arg; return; }

		int get_perturbation_level() const { return perturbation_level;}
		void set_perturbation_level(int arg) { perturbation_level = arg; return; }

		double get_perturbation_percentage() const { return perturbation_percentage;}
		void set_perturbation_percentage(double arg) { perturbation_percentage = arg; return; }

		bool get_apply_cuts() const { return apply_cuts;}
		void set_apply_cuts(bool arg) { apply_cuts = arg; return; }

		int get_serialization_size() const { return serialization_size;}
		void set_serialization_size(int arg) { serialization_size = arg; return; }
		my_objective_function get_func() const {return func;}



// 		vector<double> get_objective_function() const { return objective_function;}
// 		void set_objective_function(vector<double> arg) { objective_function = arg; objective_function_size = objective_function.size();
// 		func = my_objective_function(objective_function);
// 		return;}

		void set_objective_function(vector<double> arg) {
			func = my_objective_function(arg);
			return;
		}

		vector<double> get_objective_function(int size ) const {
			vector<double> retval;
			if (func.get_objective_function(retval, size)==0)
				return retval;


		}


// 		int get_objective_function_size() const { return objective_function_size;}
// 		void set_objective_function_size(int arg) { objective_function_size = arg; return; }



		void OVPR_to_OBJECTIVE(){
			if (OVPR =='O') set_obj_type(OBJECTIVE_TYPE_ORIGINAL);
			if (OVPR =='H') set_obj_type(OBJECTIVE_TYPE_RANDOM_FROM_HIT_RUN);
			if (OVPR =='S') set_obj_type(OBJECTIVE_TYPE_RANDOM_FROM_SHORT_DIKIN);
			if (OVPR =='L') set_obj_type(OBJECTIVE_TYPE_RANDOM_FROM_LONG_DIKIN);

			return;
		}






		void print(std::ostream &out= std::cout){
			out << OVPR << endl;
			out << input_type << endl;
			out << new_rhs_for_problem << endl;
			out << time_limit << endl;
			out << serialization_size << endl;

			for( int i = 0 ; i < serialization_size; ++i){
				out << "rowcut " << i << endl;
				rowcutlist[i].print(out);
				out << "---------" <<endl;
			}
			//out << cut_type << endl;
			out << change_obj << endl;
			out << perturbation_type << endl;
			out << perturbation_level << endl;
			out << apply_cuts << endl;
			return;
			}

		void add(mySerializableRowCut c ){
			rowcutlist.push_back(c);
			serialization_size++;
			return;
		}
		void add(vector<mySerializableRowCut> c){
			rowcutlist.insert(rowcutlist.end(),c.begin(),c.end());
			serialization_size += c.size();
			return;
		}

	};

	class Pcp_executable{
		OsiSolverInterface *mip;
		string filename;
		bool initialized;
		int ncols;
		unsigned seed;
		vector<double> original_objective_coefficients;

		int index_of_original_objective_cutoff_constraint;

// 		#ifdef USING_IOPTIMIZE
// 		ioptimize_wrapper iop_caller;
// 		#endif

		function <int(int nm, vector<double> &zikkim, int perturbation_level)> sampler;

		void initialize(int solver_id = 0){
			mip = new OsiClpSolverInterface;

			const char *f_name_lp = filename.c_str();
			sampler = bind(sampleNN, placeholders::_1, placeholders::_2, placeholders::_3, seed);
			if(strcmp(&(f_name_lp[strlen(f_name_lp)-3]), ".lp") == 0) {
				mip->readLp(f_name_lp);
			}
			else {
				if(strcmp(&(f_name_lp[strlen(f_name_lp)-4]), ".mps") == 0) {
					mip->readMps(f_name_lp);
				}
				else {
					printf("### ERROR: unrecognized file type\n");
					exit(1);
				}
			}

			ncols = mip->getNumCols();
			const double * obj = mip->getObjCoefficients();



			vector<int> indices;
			indices.resize(ncols);
			for (int i = 0;i <ncols;++i){
				original_objective_coefficients.push_back(obj[i]);
				indices[i] = i;
			}

			index_of_original_objective_cutoff_constraint = mip->getNumRows();
			mip->initialSolve();
			mySerializableRowCut tasd(indices,original_objective_coefficients, -1 * mip->getInfinity() , 1e30);

			CoinPackedVector asd = tasd.generate_CoinPackedVector();
// 			mip->writeLp("few");
			mip->addRow(asd,tasd.get_lb(),tasd.get_ub());
// 			mip->writeLp("few2");

// 			#ifdef USING_IOPTIMIZE
// 			iop_caller.set_seed(seed);
// 			iop_caller.initialize(filename);
// 			iop_caller.use_paramOpt();
//
// 			#endif

			initialized = true;

		}


// 		void set_walk_type(Pcp_input &in, std::ostream &out = std::cout){
// 			#ifdef USING_IOPTIMIZE
//
// 			if (in.get_OVPR() == 'S') {
// 				#ifdef DEBUG_my_FP
// 				out << " line " << __LINE__ << endl;
// 				#endif
// 				iop_caller.set_walk_type(1);
// 				#ifdef DEBUG_my_FP
// 				out << " line " << __LINE__ << endl;
// 				#endif
//
// 			}
// 			#ifdef DEBUG_my_FP
// 			out << " line " << __LINE__ << endl;
// 			#endif
//
// 			if (in.get_OVPR() == 'L'){
// 				#ifdef DEBUG_my_FP
// 				out << " line " << __LINE__ << endl;
// 				#endif
// 				iop_caller.set_walk_type(2);
// 				#ifdef DEBUG_my_FP
// 				out << " line " << __LINE__ << endl;
// 				#endif
//
// 			}
// 			if (in.get_OVPR() == 'H'){
// 				#ifdef DEBUG_my_FP
// 				out << " line " << __LINE__ << endl;
// 				#endif
// 				iop_caller.set_walk_type(3);
// 				#ifdef DEBUG_my_FP
// 				out << " line " << __LINE__ << endl;
// 				#endif
//
// 			}
//
// 			#endif
// 			return;
//
// 		}

		OsiCuts get_cuts(Cut_type in_cut_type){
			OsiCuts retval;

			if(in_cut_type.get_apply_gomory()){
				CglGomory cut_generator;
				cut_generator.generateCuts(*mip,retval);
			}
			if(in_cut_type.get_apply_knapsack()){
				CglKnapsackCover cut_generator;
				cut_generator.generateCuts(*mip,retval);
			}
			if(in_cut_type.get_apply_redsplit()){
				CglRedSplit cut_generator;
				cut_generator.generateCuts(*mip,retval);
			}
			if(in_cut_type.get_apply_simple_rouding()){
				CglSimpleRounding cut_generator;
				cut_generator.generateCuts(*mip,retval);
			}
			return retval;
		}

		void vectorSetObjective(vector<double> asd){
			double * asd2 = new double[ncols];
			for(int i =0; i < ncols; ++i){
				asd2[i] = asd[i];
			}
			mip->setObjective(asd2);
			delete[] asd2;
			return;
		}

		int change_obj_v1(){
			vector<double> asd;

			int b = sampler(ncols, asd,100);
			vectorSetObjective(asd);
			return 0;
		}

		int perturb_objective(int perturbation_level){
// 			double* asd = new double[ncols];
			vector<double> perturbation;

			int b = sampler(ncols, perturbation, perturbation_level);

			for (int i = 0; i < ncols; i++){
				perturbation[i] += original_objective_coefficients[i];
			}
			vectorSetObjective(perturbation);
			return 0;
		}

		double get_integer_infeasibility(const double * solution){

 			ofstream asd("zikkim");
			double retval = 0;
			int intretval = 0;
			asd << "ncols " << ncols << endl;
			for (int i = 0; i < ncols; ++i) {
				asd << "type[" << i << "] = " << mip->isContinuous(i) << " solution[" << i << "]: " << solution[i] << endl;

// 				if (type[i] == '0'){
// // 					continue;
// 					asd << "in type = 0" << endl;
// 				}

				if(mip->isContinuous(i)) continue;

// 				if (type[i] == "1" ){
// 					if ((solution[i] < 1e-10) || (solution[i] < 1e-10))
// 						continue;
// 					else
// 						retval += solution[i];
// 				}
// 				if(type[i] == "2"){

				else {
					double roundsol = floor(solution[i] + 0.5);
					double frac = fabs(roundsol - solution[i]);

					if (frac < mip->getIntegerTolerance())
						continue;
					retval+=frac;
					intretval++;
 					asd << roundsol << " " << frac << " " << retval << " " << intretval << endl;

				}

			}
			return retval;
		}

		//increase coefficients of nonbasic, decrease coefficients of basics;
		int alternative_perturbation(int perturbation_level){
			const double * asd = mip->getReducedCost();
			return -1;
		}

		int  calculate_original_objective_value(double &retval, vector<double> solution){
			if (solution.size() != ncols) return -1;
			double obj = 0;
			for (int i = 0 ;i <ncols; ++i){
				obj+= (original_objective_coefficients[i] * solution[i]);
			}
			retval = obj;
			return 0;
		}

		int  calculate_original_objective_value(double &retval, const double * solution){
// 			if (solution.size() != ncols) return -1;
			double obj = 0;
			for (int i = 0 ;i <ncols; ++i){
				obj+= (original_objective_coefficients[i] * solution[i]);
			}
			retval = obj;
			return 0;
		}



		void apply_cuts(Pcp_input in, std::ostream & out = std::cout){
			int ncuts = in.get_rowcutlist().size();
			OsiRowCut* cutlist = new OsiRowCut[ncuts];
			OsiRowCut *asd;


			for (int i =0; i< in.get_rowcutlist().size();++i){
//  				if (i == 0) continue;
// 				out << "cut " << i << endl ;
				OsiRowCut asd2 = in.get_rowcutlist()[i].generate_OsiRowCut(out);
// 				out << "printing cut " << i << endl ;
				in.get_rowcutlist()[i].print(out);
// 				out << "cut printed "<< endl;
// 				out.flush();

				cutlist[i] = in.get_rowcutlist()[i].generate_OsiRowCut(out);
// 				asd = &cutlist[i];
// 				mip->applyRowCuts(1,asd);
// 				out << "cut " << i << endl ;
// 				out.flush();
// 				asd = &asd2;
// 				mip->applyRowCuts(1,asd);
// 				mip->applyRowCuts(1,in.get_rowcutlist()[i].generate_OsiRowCutPointer());
			}
// 			out.flush();
// 			out << "apply_cut "<< endl;
// 			stringstream out1;
// 			out1 << "before.lp";
// 			mip->writeLp(out1.str().c_str());

			mip->applyRowCuts(ncuts,cutlist);
// 			stringstream out2;
// 			out2 << "after.lp";
// 			out << "write LP to "<< out2 << endl;

//  			mip->writeLp(out2.str().c_str());
			delete[] cutlist;

// 			CglPreProcess preprocessor;
//
// 			mip = preprocessor.preProcess(mip);

			return;
		}




	public:
		friend class boost::serialization::access;
		template<class Archive>
		void serialize(Archive &ar, const unsigned int version){
			ar & filename;
			ar & initialized;
			ar & seed;
		}

		Pcp_executable(string arg_filename):filename(arg_filename),initialized(false), seed(my_default_seed){}

		Pcp_executable(string arg_filename, unsigned seedx):filename(arg_filename),initialized(false), seed(seedx){}

		Pcp_executable(){}

		void print (std::ostream &out = std::cout){
			return;
		}

		void set_seed(unsigned arg) { seed = arg; return;}



		Pcp_output execute(Pcp_input in, int function_id){

			Pcp_output output_of_execute(1);

			OsiCuts cuts;

			stringstream s,s2;
			s << "ASDgrr"<< seed ;
 			s2 << "qwe" << seed;
			ofstream asd3(s.str().c_str());


			if (!initialized){
				initialize();

			}
			#ifdef USING_IOPTIMIZE
				set_walk_type(in);
			#endif

			if (mip->getRowUpper()[index_of_original_objective_cutoff_constraint] > in.get_new_rhs_for_problem() ){
				mip->setRowUpper(index_of_original_objective_cutoff_constraint,in.get_new_rhs_for_problem());
			}

// 			mip->writeLp("few3");

			asd3 << "in.get_apply_cuts " << in.get_apply_cuts() << endl;
			if(in.get_apply_cuts()){

				int ncuts = in.get_rowcutlist().size();
				asd3 << "grr ncuts " << ncuts << endl;
				if(ncuts>0){
					apply_cuts(in);
				}
			}
			asd3 << "writing LP" << endl;
			mip->writeLp(s2.str().c_str());

			#ifdef _VERSION_1_FOR_CHANGING_OBJECTIVE

			if (in.get_change_obj()>0){
				asd3 << __LINE__ << endl;

				change_obj_v1();
				asd3 << __LINE__ << endl;

			}
// 			if(in.get_perturbation_type()>0){
// 				if (in.get_perturbation_level() > 0){
// 					perturb_objective(in.get_perturbation_level());
// 				}
// 			}

			if (in.get_func().get_size() > 0){
				asd3 << __LINE__ << endl;
				// 				vectorSetObjective(in.get_objective_function());
				vector<double> x;
				asd3 << __LINE__ << endl;
				asd3 << in.get_func().get_objective_function(x,ncols,asd3) << endl;
				asd3 << __LINE__ << endl;
				vectorSetObjective(x);
				asd3 << __LINE__ << endl;
				asd3 << __LINE__ << endl;
			}
			asd3 << __LINE__ << endl;


			if(initialized) mip->initialSolve();
			else mip->resolve();
			asd3 << __LINE__ << endl;

			if (mip->isProvenDualInfeasible() || mip->isProvenPrimalInfeasible()){
				asd3 << __LINE__ << endl;
				output_of_execute.set_is_feasible(false);
// 				ofstream as3(s2.str().c_str());
// 				as3 << " infeasible " << output_of_execute.get_is_cut() << output_of_execute.get_is_incumbent()<< endl;
				mip->writeLp("few3");
				asd3 << __LINE__ << endl;

				return output_of_execute;
			}
			asd3 << __LINE__ << endl;

			#else





			#endif

// 			ofstream asd2("psd");
// 			char* type =  new char[mip->getNumCols()];
//
// 			asd2 << type[0] << endl;
// 			for (int i =0; i < ncols; ++i){ asd2 << "mip->getColType()[" << i << "]: " <<  mip->getColType()[i]  << endl;}
// 			asd3 << "get get_integer_infeasibility: "<< get_integer_infeasibility(mip->getColSolution()) << endl;

			if(get_integer_infeasibility(mip->getColSolution()) < mip->getIntegerTolerance()){
				double originalObj;

				if (calculate_original_objective_value(originalObj, mip->getColSolution())==0){
					mySerializableSolution s(ncols, mip->getColSolution());
					s.set_original_objective(originalObj);
					s.set_auxilary_objective(mip->getObjValue());
					output_of_execute.add(s);
				}
			}

			else{
				cuts = get_cuts(in.get_cut_type());

				int number_of_cuts = cuts.sizeCuts() ;

				#ifdef USE_SET

				set<mySerializableRowCut> my_cut_set;

				for (int i = 0; i < number_of_cuts; ++i){
					mySerializableRowCut cut3(cuts.rowCut(i));
					my_cut_set.insert(cut3);
				}

				output_of_execute.add(my_cut_set);

				#endif

				#ifdef USE_VEC

				vector<mySerializableRowCut> my_cut_list;

				for (int i = 0; i < number_of_cuts; ++i){
					mySerializableRowCut cut3(cuts.rowCut(i));
					my_cut_list.push_back(cut3);
				}

// 				sort(my_cut_list.begin(),my_cut_list.end());
// 				std::vector<mySerializableRowCut>::iterator it;

// 				it = unique(my_cut_list.begin(),my_cut_list.end());


// 				my_cut_list.resize( std::distance(my_cut_list.begin(),it) );

				output_of_execute.add(my_cut_list);

				#endif


			}

			return output_of_execute;
		}


	};






#endif
// kate: indent-mode cstyle; indent-width 1; replace-tabs on; 
