#include "cpfp_executable.h"
#include <OsiSolverInterface.hpp>
#include <../../Osi/src/Osi/OsiSolverInterface.hpp>


// #define VIRTUAL_BOUND 	1073741824
#define VIRTUAL_BOUND 	1e9
// #define DEBUG_ON_THE_FLY


// #define read_from_file

using namespace std;

int FP_history::compare_last_with_previous ( bool use_hist )
{


    if ( use_hist ) {
        bool last_two_the_same = true;
        bool the_one_before_last_the_same = true;

        if ( solution_history.size() <2 ) {
            return 31000;    /** compares latest two */
        }
        vector<double> previous_rounded = solution_history[solution_history.size() - 2].get_solution_vector();
        vector<double> latest_rounded = solution_history[solution_history.size() - 1].get_solution_vector();
        for ( unsigned i = 0; i < previous_rounded.size(); ++i ) {
            if ( double_inequality ( previous_rounded[i],latest_rounded[i] ) ) {
                last_two_the_same = false;
                break;
            }
        }
        if ( last_two_the_same ) {
            return 1;
        }
        if ( solution_history.size() <3 ) {
            return 30000;    /** compares latest two */
        }
        vector<double> two_previous_rounded = solution_history[solution_history.size() - 3].get_solution_vector();
        for ( unsigned i = 0; i < two_previous_rounded.size(); ++i ) {
            if ( double_inequality ( two_previous_rounded[i],latest_rounded[i] ) ) {
                the_one_before_last_the_same = false;
                break;
            }
        }
        if ( the_one_before_last_the_same ) {
            return 2;
        } else {
            return 29000;
        }



    }
    else {
        if ( int_frac_list.size() <= 1 ) {
            return 32767;
        }
        int  last_id = int_frac_list.size() - 1;
        for ( int current_id = last_id - 1 ; current_id >= 0; current_id-- ) {
            //cout << "last id : "<< last_id << " current_id : " << current_id << endl;
            bool ids = compare_ids ( last_id,current_id );
            if ( ids ) {
                return ( last_id - current_id );
            } else {
                //                 cout << int_frac_list[last_id] << " "<< int_frac_list[current_id] << " " <<
                //                         double_frac_list[last_id] << " "<< double_frac_list[current_id] << endl;

            }
        }
        return 32766;
    }

}


int FP_history::compare_last_with_previous_v1 ( bool use_hist )
{
	if ( use_hist ) {

		if (solution_history.size() <= 1){
			return 29000;
		}
		unsigned last_id = solution_history.size() - 1;
		unsigned size = solution_history.size() - 1;
		for (unsigned i = 2; i < solution_history.size();++i){
			if (solution_history[last_id] == solution_history[size - i]){
				return i;
			}

		}
		return 30000;
	}
	else{
		if ( int_frac_list.size() <= 1 ) {
			return 32767;
		}
		int  last_id = int_frac_list.size() - 1;
		for ( int current_id = last_id - 1 ; current_id >= 0; current_id-- ) {
			//cout << "last id : "<< last_id << " current_id : " << current_id << endl;
			bool ids = compare_ids ( last_id,current_id );
			if ( ids ) {
				return ( last_id - current_id );
			}
		}
		return 32766;


	}
	return -1;

}
bool FP_history::compare_ids ( unsigned i1, unsigned i2 )
{
    if ( int_frac_list[i1]!=int_frac_list[i2] ) {
        return false;
    }
    return double_equality ( double_frac_list[i1], double_frac_list[i2] );

}


void FP_history::clear ( unsigned int keep_latest_n )
{
    if ( keep_latest_n>0 ) {
        if ( keep_latest_n < int_frac_list.size() ) {
            int_frac_list.erase ( int_frac_list.begin(), int_frac_list.end() - keep_latest_n );
        }
        if ( keep_latest_n < double_frac_list.size() ) {
            double_frac_list.erase ( double_frac_list.begin(), double_frac_list.end() - keep_latest_n );
        }
        if ( keep_latest_n < solution_history.size() ) {

            solution_history.erase ( solution_history.begin(), solution_history.end() - keep_latest_n );
        }
    } else {
        int_frac_list.clear();
        double_frac_list.clear();
        solution_history.clear();
    }


}


int  cpfp_executable::find_obj_value_of_a_solution ( double &retval, const My_solution &sol, const My_objective_function &obj )
{
    /** */




	if ( sol.get_original_size() != obj.get_original_size() ) {
		    cout <<" +++++ " << endl;

        return 1;
    }

    vector<int> sol_ind = sol.get_indices();

    vector<int> obj_ind = obj.get_indices();

// 	vector_print<int> (obj_ind, cout);
// 	vector_print<int> (sol_ind, cout);


    vector<double> sol_ele = sol.get_elements();
    vector<double> obj_ele = obj.get_elements();

// 	vector_print<double> (obj_ele, cout);
// 	vector_print<double> (sol_ele, cout);


    retval = 0;
// 	cout << " line: " << __LINE__ << " sol_ind.size(): "  << sol_ind.size() << endl;
// 	cout << " line: " << __LINE__ << " obj_ind.size(): "  << obj_ind.size() << endl;

    for ( unsigned i = 0, j = 0; ( i < sol_ind.size() ) && ( j < obj_ind.size() ); ) {
// 		cout << " i: " << i << " j: " << j << " retval: " << retval <<endl;

// 		cout << "sol["<< i << "]: " << sol_ind[i] << "  " << "obj["<< j << "]: " << obj_ind[j] << endl;
// 		usleep(1000000);
        if ( sol_ind[i] < obj_ind[j] ) {
            ++i;
            continue;
        }
//         else {
			if ( sol_ind[i] > obj_ind[j] ) {
            	++j;
            	continue;
        	}
//         	else{
// 			cout << " line: " << __LINE__ << " retval: " << retval<< endl;

				retval += ( sol_ele[i] * obj_ele[j] );
				++i;
				++j;
// 			}
// 		}
    }
//     cout << " line: " << __LINE__ << endl;
    return 0;
}


Generate_cut_type::Generate_cut_type()
{
    for ( int i= 0; i<CUT_TYPE_END; ++i ) {
        generate_cut[i]= 0;
    }
}

Generate_cut_type::Generate_cut_type ( const  Generate_cut_type &other )
{
    for ( int i= 0; i<CUT_TYPE_END; ++i ) {
        generate_cut[i]= other.generate_cut[i];
    }
}

Generate_cut_type &Generate_cut_type::operator= ( const Generate_cut_type &other )
{
    for ( int i= 0; i<CUT_TYPE_END; ++i ) {
        generate_cut[i]= other.generate_cut[i];
    }

    return *this;
}

void Generate_cut_type::generate_all()
{
    for ( int i= 0; i<CUT_TYPE_END; ++i ) {
        generate_cut[i]= 1;
    }
    return;
}

void Generate_cut_type::generate_none()
{
    for ( int i= 0; i<CUT_TYPE_END; ++i ) {
        generate_cut[i]= 0;
    }
    return;
}

void Generate_cut_type::set_generate_type ( CUT_TYPE_ENM  index )
{
    generate_cut[index] = 1;
    return;
}

void Generate_cut_type::unset_generate_type ( CUT_TYPE_ENM  index )
{
    generate_cut[index] = 0;
    return;
}

int Generate_cut_type::get_generate_type ( CUT_TYPE_ENM  index )
{
    return generate_cut[index];
}
#ifdef USING_MPI
void Generate_cut_type::_send_mpi ( int destination, int tag, vector<MPI_Request> &v_request, vector<MPI_Status> v_status, MPI_SEND_TYPE SEND_TYPE )
{

    if ( tag != TAG_SENDING_GENERATE_CUT_TYPE ) {
        if ( tag != TAG_TESTING ) {
            cout << "ERROR rank: " <<MPI::COMM_WORLD.Get_rank() << " : tag: " << tag << " BUT should be TAG_SENDING_GENERATE_CUT_TYPE: " <<  TAG_SENDING_GENERATE_CUT_TYPE <<" FILE: " << __FILE__<< " LINE: " << __LINE__ << endl ;
        }
    }
    v_request.resize ( 1,MPI_REQUEST_NULL );
    v_status.resize ( 1 );

    MPI_Send ( generate_cut,CUT_TYPE_END, MPI_INT, destination,tag, MPI_COMM_WORLD );
    return;

}

void Generate_cut_type::_recv_mpi ( int source, int tag, vector<MPI_Request> &v_request, vector<MPI_Status> v_status, MPI_RECV_TYPE RECV_TYPE )
{
    if ( tag != TAG_SENDING_GENERATE_CUT_TYPE ) {
        if ( tag != TAG_TESTING ) {
            cout << "ERROR: tag: " << tag << " BUT should be TAG_SENDING_GENERATE_CUT_TYPE: " <<  TAG_SENDING_GENERATE_CUT_TYPE << " FILE: " << __FILE__<< " LINE: " << __LINE__ << endl ;
        }
    }
    v_request.resize ( 1,MPI_REQUEST_NULL );
    v_status.resize ( 1 );
    MPI_Recv ( &generate_cut[0],CUT_TYPE_END,MPI_INT, source, tag,MPI_COMM_WORLD, &v_status[0] );
    return;
}


#else
void Generate_cut_type::_send_alt1(int destination, int tag, my_signaler comm, ALT1_SEND_TYPE SEND_TYPE ){
	TAG_CONTROL(TAG_SENDING_GENERATE_CUT_TYPE)
    stringstream s; \
    s << CUT_TYPE_END << " ";
    for (unsigned i = 0; i< CUT_TYPE_END; ++i)  {s << generate_cut[i] << " " ;};
	
	//we have s as stringstream
	comm.send_object(destination, tag, s.str(),SEND_TYPE);
	return;

}
void Generate_cut_type::_recv_alt1(int source, int tag, my_signaler comm, ALT1_RECV_TYPE RECV_TYPE ){
    TAG_CONTROL(TAG_SENDING_GENERATE_CUT_TYPE)

	int rtag  = tag;
	string sr = comm.recv_object(source, rtag, RECV_TYPE);	
    std::vector<std::string> s = split_string(sr, ' ');
    //assert::assert(CUT_TYPE_END = stoi(s[0]);

    for (unsigned i = 0; i< CUT_TYPE_END; ++i) {generate_cut[i] = stoi(s[i+1]);};

    return;
}

#endif


void Generate_cut_type::print ( std::ostream & out )
{
    for ( int i= 0; i<CUT_TYPE_END; ++i ) {
        out << generate_cut[i] <<  " " ;
    }
    out << endl;
    return;
}

cpfp_executable::cpfp_executable() :initialized ( false ),initialized_extra(false), initialized_mipCP(false),initialized_AC_mip(false), seed ( my_default_seed ) {
// 	FP_DEFAULT_PARAMETER_SETTINGS

// 	cout << "line: " << __LINE__ << " alpha: " << FP_DBL_PARAM[ENM_FP_DOUBLE_PARAM_ALPHA_INIT] << " alpha_red: " << FP_DBL_PARAM[ENM_FP_DOUBLE_PARAM_ALPHA_QUOD]<<  endl ;

}

cpfp_executable::cpfp_executable ( string arg_filename ) : filename ( arg_filename ), initialized ( false ), initialized_extra(false), initialized_mipCP(false), initialized_AC_mip(false), seed ( my_default_seed ) {
// 	FP_DEFAULT_PARAMETER_SETTINGS
// 	cout << "line: " << __LINE__ << " alpha: " << FP_DBL_PARAM[ENM_FP_DOUBLE_PARAM_ALPHA_INIT] << " alpha_red: " << FP_DBL_PARAM[ENM_FP_DOUBLE_PARAM_ALPHA_QUOD]<<  endl ;

}

cpfp_executable::cpfp_executable ( string arg_filename, unsigned int seedx ) :filename ( arg_filename ), initialized ( false ),initialized_extra(false), initialized_mipCP(false), initialized_AC_mip(false), seed ( seedx ) {

}




int  cpfp_executable::initialize(int solver_type, OsiSolverInterface* mips)
{


 		mip 		= new OsiClpSolverInterface;
 		extra_mip 	= new OsiClpSolverInterface;

		mipCP 		= new OsiClpSolverInterface;
        
        AC_mip      = new OsiClpSolverInterface;
		


// #ifdef USE__CPLEX
		if (solver_type == 2){
			mip 		= new OsiCpxSolverInterface;
			extra_mip 	= new OsiCpxSolverInterface;
		}
// #else
		if (solver_type == 1){
			mip 		= new OsiClpSolverInterface;
			extra_mip 	= new OsiClpSolverInterface;
		}
		
// #endif
	CPXENVptr envasd;
	
	switch (solver_type){
		case 1:
			mip 		= new OsiClpSolverInterface;
			extra_mip 	= new OsiClpSolverInterface;
            AC_mip      = new OsiClpSolverInterface;
			break;
		case 2:
			mip 		= new OsiCpxSolverInterface;
			extra_mip 	= new OsiCpxSolverInterface;
            AC_mip 	= new OsiCpxSolverInterface;
			break;
		default:
			mip 		= new OsiCpxSolverInterface;
			extra_mip 	= new OsiCpxSolverInterface;
            AC_mip 	= new OsiCpxSolverInterface;
			
	
	}

	

	mip->messageHandler()->setLogLevel ( 0 );
	extra_mip->messageHandler()->setLogLevel ( 0 );
	mipCP->messageHandler()->setLogLevel ( 0 );
    AC_mip->messageHandler()->setLogLevel ( 0 );
	
	mip->setHintParam ( OsiDoDualInResolve );
	extra_mip->setHintParam ( OsiDoDualInResolve );
	mipCP->setHintParam ( OsiDoDualInResolve );
	AC_mip->setHintParam ( OsiDoDualInResolve );
	
// 	mip->setHintParam(OsiDoPrimalInResolve);
	if(mips != nullptr){
		cout << "mips not nullptr" << endl; 
		//mips = new OsiIpoptSolverInterface;
	}
	else {
		cout <<  "mips IS nullptr" << endl; 
	}
    const char *f_name_lp = filename.c_str();
// 		sampler = bind(sampleNN, placeholders::_1, placeholders::_2, placeholders::_3, seed);
	
    if ( strcmp ( & ( f_name_lp[strlen ( f_name_lp )-3] ), ".lp" ) == 0 ) {
		mip->readLp ( f_name_lp );
		extra_mip->readLp ( f_name_lp );
		mipCP->readLp ( f_name_lp );
        AC_mip->readLp ( f_name_lp );

    } else {
        if ( strcmp ( & ( f_name_lp[strlen ( f_name_lp )-4] ), ".mps" ) == 0 ) {
			mip->readMps ( f_name_lp );
			extra_mip->readMps ( f_name_lp );
			mipCP->readMps ( f_name_lp );

            AC_mip->readMps ( f_name_lp );

        } else {
            printf ( "### ERROR: unrecognized file type\n" );
            exit ( 1 );
        }
    }


    ncols = mip->getNumCols();
    //cout << "ncols "  << ncols << endl;
    const double * obj = mip->getObjCoefficients();
    
    original_objective_direction = mip->getObjSense();
// 	const char* col_types = mip->getColType();
// // 	cout << "col_types " << col_types<< endl;
// 	cout << "mip->getColType() " << mip->getColType()<< endl;
// 	cout << "col_types1 " << col_types[1]<< endl;
// 	cout << "col_types2 " << col_types[2]<< endl;

    column_types.resize ( ncols,0 );
    for ( int i = 0; i <ncols; ++i ) {
        
        if ( mip->isBinary ( i ) ) {
            column_types[i] =1;
            fp_status.binary_indices.push_back ( i );
			integer_columns_in_slave.push_back(i);
        } else if ( mip->isIntegerNonBinary ( i ) ) {
            column_types[i] = 2;
            fp_status.general_integer_indices.push_back ( i );
			integer_columns_in_slave.push_back(i);

        }
        //cout << "column_types[i]: " << column_types[i]<< endl;
        //cout << "obj[i]: " << obj[i]<< endl;
    }
    
// 	mip->addCol();

    vector<double> original_objective_coefficients;

    int retval = 1;

    vector<int> indices;
    indices.resize ( ncols );
    for ( int i = 0; i <ncols; ++i ) {
        original_objective_coefficients.push_back ( obj[i] );
        indices[i] = i;
        /** calculating if all objective fucntions are integer or not */
        if ( retval == 1 ) {
            if ( ( ( column_types[i]==0 ) && ( double_inequality ( original_objective_coefficients[i], 0, 1e-6 ) ) ) ) {
                /** continious variable with nonzero objective */
                retval = 0;
            } else if ( !is_integer ( original_objective_coefficients[i] ) ) {
                retval = 0;
            }
        }
    }

    const double *lbs = mip->getColLower();
	const double *ubs = mip->getColUpper();

	original_lbs.resize(ncols);
	original_ubs.resize(ncols);

	for (int i =0;i < ncols;++i){
		original_lbs[i] = lbs[i];
		original_ubs[i] = ubs[i];
	}


    original_objective_function = My_objective_function ( original_objective_coefficients );
    auxilary_objective_function = My_objective_function ( original_objective_coefficients );

    index_of_original_objective_cutoff_constraint = mip->getNumRows();
	value_of_original_objective_cutoff_constraint_ub = mip->getInfinity();
    value_of_original_objective_cutoff_constraint_lb = -1 * mip->getInfinity();
    mySerializableRowCut tasd ( indices,original_objective_coefficients,  value_of_original_objective_cutoff_constraint_lb, value_of_original_objective_cutoff_constraint_ub );

    CoinPackedVector asd = tasd.generate_CoinPackedVector();

	mip->addRow ( asd,tasd.get_lb(),tasd.get_ub() );
	mipCP->addRow ( asd,tasd.get_lb(),tasd.get_ub() );
	extra_mip->addRow ( asd,tasd.get_lb(),tasd.get_ub() );

	AC_mip->addRow ( asd,tasd.get_lb(),tasd.get_ub() );
	double *asd2 = new double[AC_mip->getNumCols()];
	for ( int i = 0; i< AC_mip->getNumCols(); ++i ) {
		asd2[i] =0;	
	}
	AC_mip->setObjective(asd2);
	
	if(mips != nullptr){
		mips = AC_mip;

		mips->initialSolve();
		const double *s = mips->getColSolution();
		vector<double> retval;
		retval.resize ( ncols );
		cout << " centerL :" ;
		for ( int i = 0; i <ncols; ++i ) {
			retval[i]=s[i];
			cout << s[i] << ", ";
		}
		cout << endl ;
		free(asd2);		
		
	}


    initialized = true;
	initialized_extra = true;
	initialized_mipCP = true;
    initialized_AC_mip = true;

	initial_solved = false;
	initial_solved_extra = false;
	initial_solved_mipCP = false;
    initial_solved_AC = false;


// 	rounded_out.open ("roundedIntVars.txt", ios_base::app);

//  	relaxed_out.open("relaxedIntVars.txt",ios_base::app);
    return retval;
}

My_objective_function cpfp_executable::get_auxilary_objective_function() const
{
    return auxilary_objective_function;
}

My_objective_function cpfp_executable::get_original_objective_function() const{
	return original_objective_function;
}

unsigned int cpfp_executable::get_seed() const
{
    return seed;
}

void cpfp_executable::set_seed ( unsigned int arg )
{
    seed = arg;
    return;
}

void cpfp_executable::set_objective_function ( const My_objective_function& new_obj, int which_mip =1 )
{

// 	const double *obj = mip->getObjCoefficients();
// 	cout << "objs : " ;
// 	for ( int i = 0; i< mip->getNumCols(); ++i ) {
// 		cout << " " << obj[i];
// 	}
// 	cout << endl;
// 	cout << endl;


// 	if (ncols != mip->getNumCols())
// 		cout << " ncols: "<< ncols << " mip->getNumCols(): "<<mip->getNumCols() << endl;
// 		cout << " line: " << __LINE__ << endl;

	auxilary_objective_function = new_obj;
// 	cout << " line: " << __LINE__ << endl;

	double *asd2 = new double[mip->getNumCols()];
	for ( int i = ncols; i< mip->getNumCols(); ++i ) {
		asd2[i] =0;
	}

// 	cout << " line: " << __LINE__ << endl;
    if ( new_obj.get_objective_function ( asd2 ) == -1 ) {
        cout << "-----------------" <<endl;
        cout << "-----------------" <<endl;
        cout << "-----------------" <<endl;
        cout << "-----------------" <<endl;
    }
//     cout << " line: " << __LINE__ << endl;

    switch (which_mip){
		case 1:
			mip->setObjective ( asd2 );
			break;
		case 2:
			extra_mip->setObjective ( asd2 );
			break;
		case 3:
			mipCP->setObjective ( asd2 );
			break;
        case 4:
            AC_mip->setObjective ( asd2 );
			break;
            
    }

//     if (which_mip == 1 ) mip->setObjective ( asd2 );
// 	if (which_mip == 2 ) extra_mip->setObjective ( asd2 );
// 	if (which_mip == 3 ) mipCP->setObjective ( asd2 );



// 	mip->setObjCoeffSet();


// 	obj = mip->getObjCoefficients();
// 	cout << "objs : " ;
// 	for ( int i = 0; i< mip->getNumCols(); ++i ) {
// 		cout << " " << obj[i];
// 	}
// 	cout << endl;
// 	cout << endl;



    delete[] asd2;
    return;
}


void cpfp_executable::set_FP_objective_function ( const My_objective_function& new_obj )
{

    FP_objective_function = new_obj;
	double *asd2 = new double[mip->getNumCols()];
	for ( int i = ncols; i< mip->getNumCols(); ++i ) {
		asd2[i] =0;
	}

    if ( new_obj.get_objective_function ( asd2 ) == -1 ) {

    }
    mip->setObjective ( asd2 );

    delete[] asd2;
    return;

}





int cpfp_executable::calculate_original_objective_value ( double& retval, const double* solution )
{
    return original_objective_function.calculate_value ( retval,solution );
}

int cpfp_executable::calculate_original_objective_value_of_solution(double &retval, My_solution & sol){
    int ret =  find_obj_value_of_a_solution(retval,sol,original_objective_function);

    if (ret == 0){
        sol.set_original_objective(retval);

    }
    return ret;
}


void cpfp_executable::print_smt_to_somewhere ( ostream& out )
{
	const double *obj = mip->getObjCoefficients();
	cout << "objs : " ;
	for ( int i = 0; i< mip->getNumCols(); ++i ) {
		cout << " " << obj[i] << " X_" << i;
	}
	cout << endl;
}


OsiCuts cpfp_executable::generate_cuts ( Generate_cut_type in_cut_type_to_generate )
{
// 	IFSLAVE_y(1) cout << " runner generate_cuts line  " <<__LINE__<<  endl;


// 	bool a;
// 	mipCP = mip->clone(a);



    OsiCuts retval;
    if ( in_cut_type_to_generate.get_generate_type ( CUT_TYPE_GOMORY ) > 0 ) {
// 		IFSLAVE_y(1) cout << " runner generate_cuts line  " <<__LINE__<<  endl;

        CglGomory cut_generator;
        cut_generator.generateCuts ( *mipCP,retval );
//         IFSLAVE_y ( 3 ) cout << 1<< endl;
    }
    if ( in_cut_type_to_generate.get_generate_type ( CUT_TYPE_KNAPSACK ) > 0 ) {
// 		IFSLAVE_y(1) cout << " runner generate_cuts line  " <<__LINE__<<  endl;

        CglKnapsackCover cut_generator;
        cut_generator.generateCuts ( *mipCP,retval );
//         IFSLAVE_y ( 3 ) cout << 2<< endl;
    }
    if ( in_cut_type_to_generate.get_generate_type ( CUT_TYPE_REDSPLIT ) > 0 ) {

// 		IFSLAVE_y(1) cout << " runner generate_cuts line  " <<__LINE__<<  endl;

        if ( mipCP->basisIsAvailable() ) {
            CglRedSplit cut_generator;
            cut_generator.generateCuts ( *mipCP,retval );
        }
//         IFSLAVE_y ( 3 ) cout << 3<< endl;
    }
    if ( in_cut_type_to_generate.get_generate_type ( CUT_TYPE_SIMPLE_ROUNDING ) > 0 ) {

// 		IFSLAVE_y(1) cout << " runner generate_cuts line  " <<__LINE__<<  endl;

        CglSimpleRounding cut_generator;
// 		IFSLAVE_y(1) cout << " runner generate_cuts line  " <<__LINE__<<  endl;

        cut_generator.generateCuts ( *mipCP,retval );
// 		IFSLAVE_y(1) cout << " runner generate_cuts line  " <<__LINE__<<  endl;

//         IFSLAVE_y ( 3 ) cout << 4<< endl;
    }

    if ( in_cut_type_to_generate.get_generate_type ( CUT_TYPE_ALL_DIFFERENT ) > 0 ) {
        CglAllDifferent cut_generator;
// 		cout << " \n\n\n WORKING \n\n\n" << endl;
        cut_generator.generateCuts ( *mipCP,retval );

    }

    if ( in_cut_type_to_generate.get_generate_type ( CUT_TYPE_L_AND_P ) > 0 ) {
        CglLandP cut_generator;

// 		cut_generator;
        // 		cout << " \n\n\n WORKING \n\n\n" << endl;
        cut_generator.generateCuts ( *mipCP,retval );
    }

    if ( in_cut_type_to_generate.get_generate_type ( CUT_TYPE_LIFT_AND_PROJECT ) > 0 ) {
        CglLandP cut_generator;

        // 		cut_generator;
        // 		cout << " \n\n\n WORKING \n\n\n" << endl;
        cut_generator.generateCuts ( *mipCP,retval );
    }

// 	if (in_cut_type_to_generate.get_generate_type(CUT_TYPE_MIR_ROUNDING) > 0){
//
//
// 		CglMixedIntegerRounding cut_generator;
//
// 		cut_generator.generateCuts ( *mip,retval );
//     }

    return retval;
}
#ifdef USE_SET
set< mySerializableRowCut> cpfp_executable::generate_rowcutlist ( Generate_cut_type in_cut_type_to_generate )
{
    OsiCuts cuts = generate_cuts ( in_cut_type_to_generate );
    int number_of_cuts = cuts.sizeCuts() ;
// 	#ifdef USE_SET

    set<mySerializableRowCut> my_cut_set;

    for ( int i = 0; i < number_of_cuts; ++i ) {
        mySerializableRowCut cut3 ( cuts.rowCut ( i ) );
        my_cut_set.insert ( cut3 );
    }
// 	#endif
    return my_cut_set;


}
#else
vector< mySerializableRowCut> cpfp_executable::generate_rowcutlist ( Generate_cut_type in_cut_type_to_generate, vector<mySerializableRowCut> &allreceivedcuts, int &previousnumber_of_cuts )
{
//     IFSLAVE_y(1) cout << " runner generate_rowcutlist line  " <<__LINE__<<  endl;
    OsiCuts cuts = generate_cuts ( in_cut_type_to_generate );
//     IFSLAVE_y(1) cout << " runner generate_rowcutlist line  " <<__LINE__<<  endl;

    int number_of_cuts = cuts.sizeCuts() ;

    vector<mySerializableRowCut> my_cut_set;
//     IFSLAVE_y(1) cout << " runner generate_rowcutlist line  " <<__LINE__<<  endl;

    for ( int i = previousnumber_of_cuts; i < number_of_cuts; ++i ) {
        mySerializableRowCut cut3 ( cuts.rowCut ( i ) );

        cut3.calculate_efficiency ( mip->getColSolution() );
        cut3.calculate_efficiency_wrt_objective ( original_objective_function,true );
        cut3.calculate_efficiency_wrt_objective ( auxilary_objective_function,false );
        cut3.update_efficiency_integral_support ( column_types, fp_status.binary_indices.size() + fp_status.general_integer_indices.size() );

        if ( cut3.get_efficiency ( CUT_EFFICIENCY_VIOLATION ) > 0 ) {
            my_cut_set.push_back ( cut3 );
        }
    }

    previousnumber_of_cuts = number_of_cuts;
    return my_cut_set;


}


#endif

void cpfp_executable::apply_cuts ( vector< mySerializableRowCut > cut_list,std::ostream &out )
{
    int ncuts = cut_list.size();
    OsiRowCut* cutlist = new OsiRowCut[ncuts];
// 	OsiRowCut *asd;


    for ( int i =0; i< ncuts; ++i ) {
//         IFSLAVE_y(3) cout << __LINE__ << endl;
// 		double temp = cut_list[i].get_norm();
// 		cut_list[i].multiply_coefficients(sqrt(temp));
// 		cut_list[i].set_norm(temp);
        cutlist[i] = cut_list[i].generate_OsiRowCut ( out );
//         IFSLAVE_y(3) cutlist[i].print();
    }

    mip->applyRowCuts ( ncuts,cutlist );
	extra_mip->applyRowCuts ( ncuts,cutlist );
	mipCP->applyRowCuts ( ncuts,cutlist );
    delete[] cutlist;
    return;
}



void cpfp_executable::update_objective_cutoff_constraint_ub ( double arg )
{
    value_of_original_objective_cutoff_constraint_ub = arg;
	mip->setRowUpper ( index_of_original_objective_cutoff_constraint,arg );
	extra_mip->setRowUpper ( index_of_original_objective_cutoff_constraint,arg );
	mipCP->setRowUpper ( index_of_original_objective_cutoff_constraint,arg );
    AC_mip->setRowUpper ( index_of_original_objective_cutoff_constraint,arg );
	

    return;
}

double cpfp_executable::get_objective_cutoff_constraint_ub () const
{
    return value_of_original_objective_cutoff_constraint_ub;
}
double cpfp_executable::get_objective_cutoff_constraint_lb () const
{
    return value_of_original_objective_cutoff_constraint_lb;
}

void cpfp_executable::update_objective_cutoff_constraint_lb ( double arg )
{
    value_of_original_objective_cutoff_constraint_lb = arg;

	mip->setRowLower ( index_of_original_objective_cutoff_constraint,arg );
	extra_mip->setRowLower ( index_of_original_objective_cutoff_constraint,arg );
	mipCP->setRowLower ( index_of_original_objective_cutoff_constraint,arg );
    AC_mip->setRowLower ( index_of_original_objective_cutoff_constraint,arg );
    return;
}

int cpfp_executable::optimize(int which_mip)
{
//#define DEBUG_ON_THE_FLY
	#ifdef DEBUG_ON_THE_FLY
	static int xxx;

	cout << "SLAVE  line " << __LINE__ << "  information: " << xxx << endl;
	xxx++;
	#endif

	bool current_initial_solved = initial_solved;
	switch(which_mip){
		case 1:
			if (!initial_solved){
				mip->initialSolve();
				initial_solved = true;
			}
			else{
				mip->resolve();
			}
			break;
		case 2:
			if (!initial_solved_extra){
				extra_mip->initialSolve();
				initial_solved_extra = true;
			}
			else{
				extra_mip->resolve();
			}
			break;
		case 3:
			if (!initial_solved_mipCP){
				mipCP->initialSolve();
				initial_solved_mipCP = true;
			}
			else{
				mipCP->resolve();
			}
			break;
		case 4:
			if (!initial_solved_AC){
				AC_mip->initialSolve();
				initial_solved_AC = true;
			}
			else{
				AC_mip->resolve();
			}
			break;
	}
	


    #ifdef DEBUG_ON_THE_FLY
    cout << "SLAVE  line " << __LINE__ << "  information: " << xxx << endl;
	xxx++;
	#endif

//     if(mip->isDualObjectiveLimitReached())



    //cout << "which_mip ==" << which_mip << endl;
	if (which_mip == 1) {OPT_RESULT(mip);}
	if (which_mip == 2) {OPT_RESULT(extra_mip);}
	if (which_mip == 3) {OPT_RESULT(mipCP);}
	if (which_mip == 4) {OPT_RESULT(AC_mip);}
	

	cout <<" THIS IS NOT WRITING " <<endl;

    if ( mip->isProvenOptimal() ) {
        return 1;
    }
    if ( mip->isProvenPrimalInfeasible() ) {
		mip->writeLp("infeas");
        return 0;
    }
    if ( mip->isDualObjectiveLimitReached() ) {
        return -2;
    }
    if ( mip->isIterationLimitReached() ) {
        return -3;
    }
    if ( mip->isPrimalObjectiveLimitReached() ) {
        return -4;
    }
    if ( mip->isProvenDualInfeasible() ) {
		mip->writeLp("isProvenDualInfeasible");

        return -5;
    }


    if ( mip->isAbandoned() ) {
		
        return -6;
    }



//     mip->is
// 	mip->solveFromHotStart();
    return -1;

}


double cpfp_executable::getObjValue(int which_mip)
{
	if (which_mip == 1) return mip->getObjValue();
	if (which_mip == 2) return extra_mip->getObjValue();
 	if (which_mip == 3) return mipCP->getObjValue();
    if (which_mip == 4) return AC_mip->getObjValue();
    
	return -1e99;
}

const double * cpfp_executable::getColSolution()
{
    return mip->getColSolution();
}

vector<double> cpfp_executable::getColSolution_asVector()
{
    const double *s = mip->getColSolution();
    vector<double> retval;
    retval.resize ( ncols );
    for ( int i = 0; i <ncols; ++i ) {
        retval[i]=s[i];
    }
    return retval;
}

int cpfp_executable::get_ncols()
{
    if ( !initialized ) {
        initialize();
    }
    return ncols;
}


double cpfp_executable::get_integer_infeasibility ( const double *solution )
{

    double retval = 0;
    int intretval = 0;
// 		asd << "ncols " << ncols << endl;
    for ( int i = 0; i < ncols; ++i ) {
        if ( mip->isContinuous ( i ) ) {
            continue;
        } else {
            double roundsol = floor ( solution[i] + 0.5 );
            double frac = fabs ( roundsol - solution[i] );

            if ( frac < mip->getIntegerTolerance() ) {
                continue;
            }
            retval+=frac;
            intretval++;
        }

    }
    return retval;
}


int cpfp_executable::get_solution ( My_solution& msolution, double &integer_infeasibility ){


// 	cout << " line: " << __LINE__ << " mip->isProvenOptimal(): " << mip->isProvenOptimal() << endl;
    const double *solution = mip->getColSolution();
    integer_infeasibility= 0;
    int intretval = 0;
// 	cout << "ncols " << ncols << " getNumcols: "<<mip->getNumCols() << endl;
//     IFSLAVE_y(3) cout << "SLAVE " << MPI::COMM_WORLD.Get_rank() << " ";


// 	cout <<" --------------" << endl;

    for ( int i = 0; i < ncols; ++i ) {
//         cout << solution[i] << " ";

//         IFSLAVE_y(3) cout << " " << solution[i] ;
        if ( mip->isContinuous ( i ) ) {
            continue;
        } else {
            double roundsol = floor ( solution[i] + 0.5 );
            double frac = fabs ( roundsol - solution[i] );

            if ( frac < mip->getIntegerTolerance() ) {
                continue;
            }
            integer_infeasibility+=frac;
            intretval++;
        }
//         cout << solution[i] << " ";

    }
//     cout << endl;
//     IFSLAVE_y(3) cout << endl;
//     if ( intretval ==0 ) {
        msolution = My_solution ( ncols,solution );
        double original_objective_value;

// 		msolution.print();
        calculate_original_objective_value ( original_objective_value,solution );
        msolution.set_original_objective ( original_objective_value );
        msolution.set_auxiliary_objective ( mip->getObjValue() );
//     }

// 		msolution.shortlineprint_selected_indices(integer_columns_in_slave,relaxed_out);
// 		msolution.shortlineprint(relaxed_out);

// 		cout <<" --------------" << endl;

		return intretval;
}



int cpfp_executable::get_solution_stage ( My_solution& msolution, double &integer_infeasibility, int stage ){
	
	
	// 	cout << " line: " << __LINE__ << " mip->isProvenOptimal(): " << mip->isProvenOptimal() << endl;
	const double *solution = mip->getColSolution();
	integer_infeasibility= 0;
	
	int intretval = 0;

	
	for ( int i = 0; i < ncols; ++i ) {
		switch(column_types[i]){
			case 0: 
				continue;
				break;
			case 2:
				if (stage ==1) break;
			case 1:
				double roundsol = floor ( solution[i] + 0.5 );
				double frac = fabs ( roundsol - solution[i] );
				
				if ( frac < mip->getIntegerTolerance() ) {
					continue;
				}
				integer_infeasibility+=frac;
				intretval++;
		}
	}

	msolution = My_solution ( ncols,solution );
	double original_objective_value;
	
	calculate_original_objective_value ( original_objective_value,solution );
	msolution.set_original_objective ( original_objective_value );
	msolution.set_auxiliary_objective ( mip->getObjValue() );
	
	return intretval;
}
int cpfp_executable::init_FP_stage2(){
    /** use the same solver interface for FP*/
    int n_general_integer = 0;
    for ( int i = 0; i < ncols; ++i ) {
        if ( column_types[i] == 2 ) {
            n_general_integer++;
        }
    }


//     cout << " line: " << __LINE__ << endl;

    fp_status.n_cols_added = fp_status.general_integer_indices.size();
    fp_status.n_rows_added = 2 * fp_status.general_integer_indices.size();
    int *  	column = new int[1];
    int * 	rows = new int[2];
    double *  	columnelement = new double[1];
    double *  	rowelement = new double[2];
// 	cout << " line: " << __LINE__ << endl;

    rowelement[0] = 1;
    rowelement[1] = -1;
    columnelement[0] = 1;
    double maxub = mip->getInfinity();
	double minlb = -maxub;

    /** replace this one with something else for a tighther bound*/

// 	int xc = fp_status.general_integer_indices.size();
	fp_status.added_row_indices1.resize(fp_status.n_rows_added);
	fp_status.added_row_indices2.resize(fp_status.n_rows_added);
// 	cout << " line: " << __LINE__ <<" fp_status.general_integer_indices.size(): " << xc << endl;

// 	cout << " get_number_of_general_integers:" << get_number_of_general_integers() << endl;
    for ( unsigned i = 0; i < fp_status.general_integer_indices.size(); ++i ) {
// 		cout << " line: " << __LINE__ << " i: " << i << endl;

        /** for all general integer variables add two rows and one column
         * for example if x_1 = a for some FP rounded solution
         * we minimize |x_1 - a|, that is we add one column x'_1
         * and two rows such that |x_1 - a| <= x'_1
         * x_1 - a  <= x'_1  and a - x_1 <= x'_1 ;
         * x_1 - x'_1 <= a   and a <= x_1 + x'_1
         * and minimize x'_1
         */
        /** adding rows for next integer */

        column[0] = fp_status.general_integer_indices[i];
        rows[0] = mip->getNumRows(); /** indices of new rows */

		rows[1] = mip->getNumRows() + 1;

		fp_status.added_row_indices1[i] = rows[0];
		fp_status.added_row_indices2[i] = rows[1];
// 		cout << " line: " << __LINE__ << " i: " << i << endl;

        mip->addRow ( 1,column, columnelement,minlb,maxub ); /** this will be row for a<= x_1 + x'_1 */
        mip->addRow ( 1,column, columnelement,minlb,maxub ); /** this will be row for x_1 - x'_1  <= a*/
// 		cout << " line: " << __LINE__ << " i: " << i << " mip->getNumRows(): " << mip->getNumRows() << endl;


        fp_status.added_col_indices.push_back ( mip->getNumCols() );
// 		cout << " line: " << __LINE__ << " i: " << i << endl;
// 		cout << " line: " << __LINE__ << " rows[0]: "<< rows[0] << endl;
// 		cout << " line: " << __LINE__ << " rows[1]: "<< rows[1] << endl;
// 		cout << " line: " << __LINE__ << " rows[1]: "<< rows[1] << endl;


        mip->addCol ( 2,rows,rowelement,minlb,maxub,0);
// 		cout << " line: " << __LINE__ << " i: " << i << endl;
    }




//     cout << " line: " << __LINE__ << endl;

    delete[] column;
    delete[] rows;
    delete[] columnelement;
    delete[] rowelement;
// 	cout << " line: " << __LINE__ << endl;

// 	const double * objs = mip->getObjCoefficients();
//
// 	for (int i = ncols ; i < mip->getNumCols();++i){
// 		cout << objs[i] << " ";
// 	}
// 	cout << endl;

	init_FP_stage1();
// 	cout << " line: " << __LINE__ << endl;
	return 0;


}
int cpfp_executable::init_FP_stage1()
{
	mip->setObjSense ( 1 );

	FP_initialized = true;
	return 0;
}

int cpfp_executable::undo_FP_stage2(){
// 	cout << " line: " << __LINE__ << " fp_status.n_cols_added: " << fp_status.n_cols_added <<  endl;
// 	cout << " line: " << __LINE__ << " fp_status.n_rows_added: " << fp_status.n_rows_added << " added_row_indices1.size: "<< fp_status.added_row_indices1.size() <<  endl;
// 	cout.flush();
    int * column_indices = new int[fp_status.n_cols_added];
// 	cout << " line: " << __LINE__ << endl;
// 	cout.flush();
    int * row_indices = new int[fp_status.n_rows_added];
// 	cout << " line: " << __LINE__ << endl;
// 	cout.flush();
    for ( int i = 0; i < fp_status.n_cols_added; ++i ) {
// 		cout << " line: " << __LINE__ << endl;

// 		cout << " line: " << __LINE__ << " fp_status.added_col_indices[i]: " << fp_status.added_col_indices[i] << endl;

// 		cout.flush();

		column_indices[i] = fp_status.added_col_indices[i];
// 		cout << " line: " << __LINE__ << " fp_status.added_row_indices1[i]: " << fp_status.added_row_indices1[i] << endl;
// 		cout << " line: " << __LINE__ << " fp_status.added_row_indices2[i]: " << fp_status.added_row_indices2[i] << endl;
// 		cout.flush();



		row_indices[i] = fp_status.added_row_indices1[i];
// 		cout << " line: " << __LINE__ << endl;
// 		cout.flush();

		row_indices[i+fp_status.n_cols_added] = fp_status.added_row_indices2[i];
// 		cout << " line: " << __LINE__ << endl;
// 		cout.flush();

    }
//     cout << " line: " << __LINE__ << endl;
// 	cout.flush();

    mip->deleteCols ( fp_status.n_cols_added, column_indices );
// 	cout << " line: " << __LINE__ << endl;
// 	cout.flush();

	mip->deleteRows ( fp_status.n_rows_added, row_indices );
// 	cout << " line: " << __LINE__ << endl;
// 	cout.flush();

    fp_status.added_col_indices.clear();
    fp_status.n_cols_added = 0;
    fp_status.n_rows_added = 0;
    fp_status.added_row_indices1.clear();
    fp_status.added_row_indices2.clear();
// 	cout << " line: " << __LINE__ << endl;
// 	cout.flush();

	undo_FP_stage1();

    return 0;
}

int cpfp_executable::undo_FP_stage1(){
	mip->setObjSense ( original_objective_direction );

	FP_initialized = false;
	set_objective_function ( original_objective_function );
	return 0;
}

double cpfp_executable::r_gen(){
    static uniform_real_distribution<double> uni ( 0.0,1.0 );
    static mt19937 prng ( seed );
    return uni ( prng );
}



/** basically round it */
int cpfp_executable::round ( vector< double >& starting_solution, vector< double >& rounded_solution, int arg_rounding_type, int stage){

	int n_changed = 0;

// 	static int rest = 0;
// 	int arg_rounding_type = arg_rounding_type2;

	switch (arg_rounding_type){
		case 0:
		case 1:
		case 2:
		case 3:
		case 4:
		case 5:
// 			arg_rounding_type = arg_rounding_type2;
// 			rest =0;
			break;
		case 6:
			arg_rounding_type = floor(r_gen()*4);
// 			rest =0;
			break;
		case 7:
			arg_rounding_type = floor(r_gen()*6);
// 			rest=0;
			break;
		case 8:
// 			rest++;
// 			if (rest<5) arg_rounding_type = 0;
			// 			else arg_rounding_type = floor(r_gen()*4);
			arg_rounding_type = floor(r_gen()*4);
			break;
		case 9:
// 			rest++;
// 			if (rest<5) arg_rounding_type = 0;
			// 			else arg_rounding_type = floor(r_gen()*6);
			arg_rounding_type = floor(r_gen()*6);
			break;
		default:
			arg_rounding_type = 0;
			break;

	}

// 	if (arg_rounding_type2 <0 || arg_rounding_type2>6) arg_rounding_type = floor(r_gen()*6);
// 	else if (arg_rounding_type2 == 6) arg_rounding_type = floor(r_gen()*4);
//

	double threshold =r_gen();
	double temp = threshold;

	switch ( arg_rounding_type ) { 
        case 0:
			threshold = 0.5;
            break;
        case 1:
            break;
        case 2:
            temp= threshold/2 + 0.25;  // threshold = rand (0.25, 0.75) 
			threshold = temp;
            break;
        case 3 :
			if (threshold<=0.5)	temp = 2*(threshold)*(1-threshold);
			else temp = 2*(threshold)*(threshold -1 ) +1;

			threshold = temp;
            break;

    };
// 	static int ok =0;
// 	ok++;
// 	if (ok <=100)
// 	cout  << "ok " << ok << " line: " << __LINE__ << " rounding_type: " << arg_rounding_type << " threshold: "<< threshold << endl;

//     const double *lowerbounds = mip->getColLower();
//     const double *upperbounds = mip->getColUpper();

// 	double temp;
	assert(arg_rounding_type <=5);

// 	if (arg_rounding_type>0) cout << "threshold: " << threshold << " arg_rounding_type: " << arg_rounding_type <<endl;
    for ( int i = 0; i < ncols; ++i ) {
        switch ( column_types[i] ) {
            case 0: /** continious */
                rounded_solution[i] = starting_solution[i];
                //             mip->setObjCoeff ( i,0 );
                break;
            case 1: /** binary*/
                if ( arg_rounding_type <= 3 ) {
					temp = floor ( starting_solution[i]+threshold );
					if (double_inequality(rounded_solution[i],temp, epInt)){
						rounded_solution[i] = temp;
						n_changed++;
					}
                } else {
                    double flip_probability = 0.5;
                    if ( arg_rounding_type == 5 ) {
                        flip_probability =.75;
                    }
                    temp = floor ( starting_solution[i] + 0.5 );
//                     rounded_solution[i] = floor ( starting_solution[i] + 0.5 );
                    if ( r_gen() <=flip_probability ) { /** FLIP*/
                        if ( double_equality ( temp,0,epInt ) ) {

                            rounded_solution[i] = 1;
                        } else {
                            rounded_solution[i] = 0;
                        }
// 						n_changed++;
                    }
                }
                //             if ( double_equality ( rounded_solution[i],0 ) ) { /** minimize x */
                //                 mip->setObjCoeff ( i, 1 );
                //             } else {
                //                 mip->setObjCoeff ( i, -1 );
                //             }
                break;
            case 2: /** general integer */
                if (stage == 1) {

                	/** it was only break; */
					temp = floor ( starting_solution[i]+ 0.5 );
					if (double_inequality(temp, rounded_solution[i],epInt)){
						rounded_solution[i] = temp;
						//n_changed++;
					}
					break;

				}
                if ( arg_rounding_type <= 3 ) {
					temp = floor ( starting_solution[i]+ threshold );
					if (double_inequality(temp, rounded_solution[i],epInt)){
						rounded_solution[i] = temp;
						n_changed++;
					}
                }
                if ( arg_rounding_type == 4 ) { /** JUST FLIP */
					/** check integrality*/

					double f = floor(starting_solution[i]+0.5);
					if (double_equality(f,starting_solution[i], epInt)){
						/** integer */
						if (r_gen() < 0.5){ /** FLIP */

							if (double_equality(starting_solution[i],original_lbs[i],epInt)){
								/** at LB */
								temp = starting_solution[i]+1;
								if (double_inequality(temp, rounded_solution[i])){
									rounded_solution[i] = temp;
									n_changed++;
								}

							}
							else{
								if (double_equality(starting_solution[i],original_ubs[i],epInt)){
									/** at UB */

									temp = starting_solution[i]-1;
									if (double_inequality(temp, rounded_solution[i])){
										rounded_solution[i] = temp;
										n_changed++;
									}
								}
								else{
								/** neither at LB or UB */

									if (r_gen()<0.5)temp= starting_solution[i] + 1;
									else temp = starting_solution[i]-1;


									if (double_inequality(temp, rounded_solution[i])){
										rounded_solution[i] = temp;
										n_changed++;
									}
								}
							}


						}
						else {
							temp=starting_solution[i];
							if (double_inequality(temp, rounded_solution[i])){
								rounded_solution[i] = temp;
								n_changed++;
							}
						}
					}
					else{
						/** not integer */
						double floor_x = floor ( starting_solution[i] );
						double ceil_x = ceil ( starting_solution[i] );
						if (r_gen()<0.5)temp = floor_x;
						else temp = ceil_x;

						if (double_inequality(temp, rounded_solution[i])){
							rounded_solution[i] = temp;
							n_changed++;
						}
					}
                }
                if ( arg_rounding_type == 5 ) { /** aggressive restart*/

                    double lblb = max ( original_lbs[i], starting_solution[i] - VIRTUAL_BOUND );
                    double ubub = min ( original_ubs[i], starting_solution[i] + VIRTUAL_BOUND );
                    temp = floor ( r_gen() * ( ubub - lblb + 1) );

					if (double_inequality(temp, rounded_solution[i])){
						rounded_solution[i] = temp;
						n_changed++;
					}
//                     cout << "\t" << rounded_solution[i] << "\t" << starting_solution[i];
                }
                //             mip->setObjCoeff ( i, 0 );
                //             mip->setObjCoeff ( fp_status.added_col_indices[j],1 );
                //             j++;
                break;
            default:

                temp = floor ( starting_solution[i]+ threshold);
				if (double_inequality(temp, rounded_solution[i],epInt)){
					rounded_solution[i] = temp;
					n_changed++;
				}
				break;
        }
    }



    return n_changed;

}



int cpfp_executable::round_sol ( My_solution& starting_solution, My_solution& rounded_solution, int arg_rounding_type, int stage )
{
	/** quick and dirty*/

	static int ss;
	
	ss++;
	
	
	vector<double> s = starting_solution.get_solution_vector();
	vector<double> r = rounded_solution.get_solution_vector();

	int retval = round(s,r,arg_rounding_type,stage);
// 	cout << "rounding for the "<< ss << "th time rounding type "<< arg_rounding_type << " changed: "<< retval<< endl;

	starting_solution = My_solution(s);
	rounded_solution = My_solution(r);



// 	rounded_solution.shortlineprint_selected_indices(integer_columns_in_slave,rounded_out);
	return retval;

}
/*
int cpfp_executable::round_sol_v2 ( My_solution& starting_solution, My_solution& rounded_solution, int& arg_rounding_type, int stage ){

	vector<double> s_e = starting_solution.get_elements();
	vector<double> r_e = rounded_solution.get_elements();

	vector<int> s_i = starting_solution.get_indices();
	vector<int> r_i = rounded_solution.get_indices();

	for (unsigned i = 0; i < rounded_solution.s)


}*/



int cpfp_executable::iterate_FP ( vector<double> &starting_solution, vector<double> &rounded_solution, int init_rounding_type )
{

    static int zs = 0;
    zs ++;
    int rounding_type = init_rounding_type;
    int retval = 0;
    double d_fractionality= 0;
    double threshold = 0.5;

    /** rounding_type   = 0; standart : threshold = 0.5 ;
     *                  = 1; completely random threshold ;
     *                  = 2; random threshold .25-.75 ;
     *                  = 3; random flip;
     * 					= 4; aggressive restart;
     *
     *                  = 9; random method ;
     */


    static int iz = 0;
    if ( rounding_type == 9 ) {
        double ty = 5*r_gen();
        rounding_type = int ( floor ( ty ) );
        iz = 0;
    }
    /** NOTE: primitive history check */
    int latest_equal = fp_hist.compare_last_with_previous();


    int n_general_integer = 0;

    for ( unsigned i = 0; i < column_types.size(); ++i ) {
        if ( column_types[i] == 2 ) {
            n_general_integer++;
        }
    }


    if ( latest_equal < COMPARE_LAST_N ) {
        iz++;
        rounding_type =3;
        if ( iz >= 100 ) {
            //cout << "cycle detected by " << seed << " " << iz << "th time " << endl;
// 			cout << "n_general_integer: " << n_general_integer << endl;
            iz = 0;
            rounding_type = 4;/** do aggressive restart every 100 flips*/
        }

    } else {

        if ( fp_hist.int_frac_list.size() >= COMPARE_LAST_N*10 ) {
            fp_hist.clear ( COMPARE_LAST_N );
        }
    }

    switch ( rounding_type ) {
    case 0:
        break;
    case 1:
        threshold = r_gen();
        break;
    case 2:
        threshold = 2*threshold -0.25;
        break;
    case 3 :
        threshold = 0;
        break;
    case 4 :
        break;

    };


    const double *lowerbounds = mip->getColLower();
    const double *upperbounds = mip->getColUpper();


    /** ROUND
    * DO NOT SET objective function coefficients yet
    * we may need to change depending on the history check*/
    int j = 0;
    for ( int i = 0; i < ncols; ++i ) {
        switch ( column_types[i] ) {
        case 0: /** continious */
            rounded_solution[i] = starting_solution[i];
//             mip->setObjCoeff ( i,0 );
            break;
        case 1: /** binary*/
            if ( rounding_type < 3 ) {
                rounded_solution[i] = floor ( starting_solution[i]+threshold );
            } else {
                double flip_probability = 0.5;
                if ( rounding_type == 4 ) {
                    flip_probability =.75;
                }

                rounded_solution[i] = floor ( starting_solution[i] + 0.5 );
                if ( r_gen() <=flip_probability ) { /** FLIP*/
                    if ( double_equality ( rounded_solution[i],0 ) ) {
                        rounded_solution[i] = 1;
                    } else {
                        rounded_solution[i] = 0;
                    }

                }
            }
//             if ( double_equality ( rounded_solution[i],0 ) ) { /** minimize x */
//                 mip->setObjCoeff ( i, 1 );
//             } else {
//                 mip->setObjCoeff ( i, -1 );
//             }
            break;
        case 2: /** general integer */
            if ( rounding_type < 3 ) {
                rounded_solution[i] = floor ( starting_solution[i]+ threshold );
            }
            if ( rounding_type == 3 ) { /** JUST FLIP */
                double floor_x = floor ( starting_solution[i] );
                double ceil_x = ceil ( starting_solution[i] );

                if ( double_equality ( floor_x, ceil_x, 0.001 ) ) { /** starting solution is integer +- 1 */
                    if ( r_gen() <=0.5 ) {/** FLIP*/
                        rounded_solution[i] = floor_x + 1;
                    } else {
                        rounded_solution[i] = floor_x - 1;
                    }
                } else { /** current non integer floor or ceiling */
                    if ( r_gen() <=0.5 ) {/** FLIP*/
                        rounded_solution[i] = floor_x;
                    } else {
                        rounded_solution[i] = ceil_x;
                    }
                }
            }
            if ( rounding_type == 4 ) { /** aggressive restart*/
                double lblb = max ( lowerbounds[i], starting_solution[i] - VIRTUAL_BOUND );
                double ubub = min ( upperbounds[i], starting_solution[i] + VIRTUAL_BOUND );
                rounded_solution[i] = floor ( r_gen() * ( ubub - lblb ) );
                cout << "\t" << rounded_solution[i] << "\t" << starting_solution[i];
            }
//             mip->setObjCoeff ( i, 0 );
//             mip->setObjCoeff ( fp_status.added_col_indices[j],1 );
//             j++;
            break;
        default:
            rounded_solution[i] = floor ( starting_solution[i]+ 0.5 );
            break;
        }
    }

    int equality_last_n =fp_hist.compare_last_with_previous ( true );

    vector<double> previous_rounded;


    if ( equality_last_n == 1 ) {
//         cout << "equality_last_n: "<<  equality_last_n << endl;

        previous_rounded = fp_hist.solution_history[fp_hist.solution_history.size()-2].get_solution_vector();


        f_type_restart_1 ( starting_solution, rounded_solution, previous_rounded );
        vector_compare<double> ( rounded_solution,previous_rounded );
    }
    if ( equality_last_n == 2 ) {
//         cout << "equality_last_n: "<<  equality_last_n << endl;

        previous_rounded = fp_hist.solution_history[fp_hist.solution_history.size()-3].get_solution_vector();
        f_type_restart_1 ( starting_solution, rounded_solution, previous_rounded );
    }



    for ( int i = 0; i < ncols; ++i ) {
        switch ( column_types[i] ) {
        case 0: /** continious */
            mip->setObjCoeff ( i,0 );/** */
            break;
        case 1:
            if ( double_equality ( rounded_solution[i],0 ) ) { /** minimize x */
                mip->setObjCoeff ( i, 1 );
            } else {
                mip->setObjCoeff ( i, -1 );
            }
            break;
        case 2:
            mip->setObjCoeff ( i, 0 );
            mip->setObjCoeff ( fp_status.added_col_indices[j],1 );
            j++;
            break;
        }
    }


    if ( j!= fp_status.n_cols_added ) {
        //CPFP_ERROR ( "changing  more or less columns than added ",cout )
        cout << "j = " << j << "  fp_status.n_cols_added: " <<		fp_status.n_cols_added << endl;
    }


    /**change bounds for added rows for general indices */
    for ( unsigned i = 0; i < fp_status.general_integer_indices.size(); ++i ) {
        mip->setRowLower ( fp_status.added_row_indices1[i], rounded_solution[fp_status.general_integer_indices[i]] );
        mip->setRowUpper ( fp_status.added_row_indices2[i], rounded_solution[fp_status.general_integer_indices[i]] );
    }

    fp_inf_itercount_plus();
    int res = optimize();

    if ( res == 1 ) {
		My_solution asd;
        int n_frac = get_solution ( asd,d_fractionality );
        if ( ( n_frac== 0 ) ) {
// 				double_equality(get_integer_infeasibility(mip->getColSolution()),0)){
            /** current solution is integer*/
			best_incumbent = asd;
            fp_hist.clear();
            retval = 1;
        } else {
            starting_solution = getColSolution_asVector();


            /** NOTE: history check later*/
        }
    } else {
        retval = res;
        if ( res == 0 ) {
            /** problem is infeasable */

        } else { /** unknown result */

        }

    }



    return retval;
}

int cpfp_executable::pre_FP ( int iterlim, int resetlim, double timelim, int rounding_type,vector<double> &starting_solution, vector<double> &rounded_solution )
{
    fp_inf.iterlim = iterlim;
    fp_inf.resetlim = resetlim;
    fp_inf.itercount = 0;
    fp_inf.resetcount = 0;
    starting_solution.resize ( ncols );
    rounded_solution.resize ( ncols );
    starting_solution = getColSolution_asVector();
    fp_hist.int_frac_list.clear();
    fp_hist.double_frac_list.clear();
    fp_hist.solution_history.clear();
    return 0;
}



int cpfp_executable::pre_FP_v1 ( int iterlim, int resetlim, int st2_iterlim, int st2_resetlim, int rounding_type,vector<double> &starting_solution, vector<double> &rounded_solution )
{
	fp_inf.iterlim = iterlim;
	fp_inf.resetlim = resetlim;
	fp_inf.itercount = 0;
	fp_inf.resetcount = 0;
	starting_solution.resize ( ncols );
	rounded_solution.resize ( ncols );
	starting_solution = getColSolution_asVector();
	fp_hist.int_frac_list.clear();
	fp_hist.double_frac_list.clear();
	fp_hist.solution_history.clear();

	fp_inf.stage2_iterlim = st2_iterlim;
	fp_inf.stage2_resetlim = st2_resetlim;
	fp_inf.stage2_itercount = 0;
	fp_inf.stage2_resetcount = 0;




	return 0;
}


int cpfp_executable::pre_FP_v1_sol ( int iterlim, int resetlim, int st2_iterlim, int st2_resetlim, int rounding_type, My_solution& starting_solution, My_solution& rounded_solution )
{
	vector<double> s;
	vector<double> r;
	pre_FP_v1(iterlim,resetlim,st2_iterlim,st2_resetlim,rounding_type,s,r);

	starting_solution = My_solution(s);
	rounded_solution = My_solution(r);
	return 0;

}



void cpfp_executable::fp_inf_itercount_plus()
{
    fp_inf.itercount++;
}

My_solution cpfp_executable::get_best_incumbent()
{
    return best_incumbent;
}

double cpfp_executable::get_best_incumbents_objective()
{
//     best_incumbent.calculate_original_objective_value();

    return best_incumbent.get_original_objective();
//     best_incumbent.
}

int cpfp_executable::post_iterate ( vector<double> &starting_solution )
{

    if ( ( fp_inf.iterlim>0 ) && ( fp_inf.itercount > fp_inf.iterlim ) ) {
        return 1;
    }
    if ( fp_inf.itercount > fp_inf.iterlim ) {
        fp_inf.resetcount++;
        fp_inf.itercount = 0;
        starting_solution= getColSolution_asVector();
        return 4;
    }
    if ( fp_inf.resetcount > fp_inf.resetlim ) {
        return 2;
    }
    return 0;
}

void cpfp_executable::reset_fp_inf()
{
    fp_inf.resetcount ++;
    fp_inf.itercount = 0;


}

void cpfp_executable::reset_fp_lim()
{
    fp_inf.resetcount = 0;
    fp_inf.itercount = 0;
}

void cpfp_executable::puke_lp_file ( const char *filename )
{
	mip->writeLp ( filename);
// 	mip->writeLp ( filename,"lp", epInt,2,18);
}

void cpfp_executable::puke_mps_file ( const char *filename )
{
    mip->writeMps ( filename );
}


int cpfp_executable::intoptimize(){
    mip->setDblParam ( OsiPrimalObjectiveLimit,0 );
    mip->setDblParam ( OsiDualTolerance,0 );
    mip->setDblParam ( OsiPrimalTolerance,0 );
    mip->writeLp ( "final" );
    mip->setIntParam ( OsiMaxNumIteration, 100 );

    mip->branchAndBound();

// 	cout << " mip->getIterationCount(); " <<mip->getIterationCount()<<endl;
//     vector<double> sol = getColSolution_asVector();
//
//     My_solution s(sol);
// 	double x;
// 	original_objective_function.calculate_value(x,mip->getColSolution());

// 	mip->setIntParam();




// 	cout << "mip->getObjValue();" << mip->getObjValue() <<endl;

// 	cout << mip->getAuxiliaryInfo();

// 	s.set_original_objective(x);
// 	s.print();
    return 0;

}

int cpfp_executable::f_type_restart_1 ( vector< double >& starting_solution, vector< double >& rounded_solution, vector< double >& previous_rounded ){
    int binary_changed = 0;

    for ( unsigned i = 0; i < column_types.size(); i++ ) {
        if ( column_types[i] == 1 ) {
            double r = r_gen() - 0.47;
// 			std::cout << __LINE__ <<" "  <<r <<std::endl;
            if ( r>0 && ( double_equality (rounded_solution[i], previous_rounded[i], epInt ) )) {
                double sigma = abs ( rounded_solution[i]-starting_solution[i] );
                if ( sigma + r > 0.5 ) { /** FLIP */
                    if ( double_equality ( rounded_solution[i],original_lbs[i] ) ) {
                        rounded_solution[i] = original_ubs[i];
                    } else {
                        rounded_solution[i] = original_lbs[i];
                    }
                    binary_changed++;
                }
            }
        }
    }
//     cout << " n_changed: " << binary_changed << endl;
	
    return binary_changed;
}

int  cpfp_executable::f_type_restart_1 ( My_solution& starting_solution, My_solution& rounded_solution, My_solution& previous_rounded )
{
	/** quick and dirty */

	vector<double> s = starting_solution.get_solution_vector();
	vector<double> r = rounded_solution.get_solution_vector();
	vector<double> p = previous_rounded.get_solution_vector();

	int retval = f_type_restart_1(s,r,p);

	rounded_solution = My_solution(r);


	return retval;
}

int  cpfp_executable::f_type_restart_2 ( My_solution& starting_solution, My_solution& rounded_solution, My_solution& previous_rounded, int iter )
{
	/** quick and dirty */


// 	My_solution copy_r = rounded_solution;
	vector<double> s = starting_solution.get_solution_vector();
	vector<double> r = rounded_solution.get_solution_vector();
	vector<double> p = previous_rounded.get_solution_vector();

	int retval = f_type_restart_2(s,r,p,iter);

// 	cout << " before restart " << endl;
// 	cout << " -------------- " << endl;
// 	rounded_solution.shortlineprint_selected_indices(integer_columns_in_slave,false);
// 	cout << " -------------- " << endl;
// 	cout << " -------------- " << endl;
// 	vector_print<double>(r,cout);
// 	cout << " -------------- " << endl;



	rounded_solution = My_solution(r);

// 	cout << " after restart " << endl;

// 	cout << " -------------- " << endl;
// 	vector_print<double>(r,cout);
// 	cout << " -------------- " << endl;

// 	cout << " -------------- " << endl;

// 	rounded_solution.shortlineprint_selected_indices(integer_columns_in_slave,false);
// 	cout << " -------------- " << endl;


// 	cout << "EQ: " << copy_r.equality_on_selected_indices(rounded_solution, integer_columns_in_slave) << endl;

	return retval;
}


int  cpfp_executable::f_type_restart_2 ( vector< double >& starting_solution, vector< double >& rounded_solution, vector< double >& previous_rounded, int iter ){
    int binary_changed = 0;
    int int_change = 0;
    int min_number_of_change = FP_INT_PARAM[ENM_FP_INT_PARAM_MINFLIP];
//     int n_to_change = min_number_of_change + int ( r_gen() * 2*min_number_of_change );

//     double sigma_min = FP_DBL_PARAM[ENM_FP_DOUBLE_PARAM_FLIP_TRESHOLD];

//     int ind = 0;
    vector<double> d_list;
    priority_queue<pair<double, int>, vector<pair<double, int> >,sort_pair_first<double,int,std::greater<double> > > my_pair_list;

	static int lastiter, rp;

// 	static int called = 0;
// 	called++;

// 	cout << "this is called " << called << " times" <<endl;
// 	usleep(1000000);

	if( lastiter < iter)
	{
		while( ++lastiter < iter )
			rp=(int)(rp*.85);
	}

// 	int n_int = get_number_of_binaries()+get_number_of_general_integers();
	rp = min((get_number_of_binaries()+get_number_of_general_integers())/10, 2*min_number_of_change+rp+1);


	int nInt = integer_columns_in_slave.size();


#ifdef read_from_file
	ofstream reec("s_cre.txt");
	ofstream reer("s_read.txt");

	ifstream ree ("../FP2/ree.txt");
	ofstream ree2 ("ree.txt");

#endif
	double temp;
	for (int j = 0; j < rp; j++){
		int int_index_to_change = (int)floor(r_gen()*(nInt));

		int to_change = integer_columns_in_slave[int_index_to_change];

#ifdef read_from_file

		stringstream s_cre,s_read;
		s_cre << " created int_index_to_change: " << int_index_to_change<< 				" to_change: "<< to_change << 				" lb: " << original_lbs[to_change] <<				" ub: " << original_ubs[to_change] << 				" cv: " << rounded_solution[to_change ];
// 		cout<< " created int_index_to_change: " << int_index_to_change<< 				" to_change: "<< to_change << 				" lb: " << original_lbs[to_change] <<				" ub: " << original_ubs[to_change] << 				" cv: " << rounded_solution[to_change ];

		double lb, ub, old_v, new_v;
		int changed;
		ree >>	int_index_to_change >> lb >> ub >> old_v >> new_v >> changed;
		to_change = integer_columns_in_slave[int_index_to_change];

		s_read << " read    int_index_to_change: " << int_index_to_change<< 		" to_change: "<< to_change << 		" lb: " << lb <<		" ub: " << ub << 		" cv: " << old_v << " uv: " << new_v << " ch: " << changed;

// 		cout << " read    int_index_to_change: " << int_index_to_change<< 		" to_change: "<< to_change << 		" lb: " << lb <<		" ub: " << ub << 		" cv: " << old_v << " uv: " << new_v << " ch: " << changed;
		ree2 << int_index_to_change << "\t" << original_lbs[to_change] <<"\t" << original_ubs[to_change] <<"\t" << rounded_solution[to_change]<< "\t";


#endif

		if ( column_types[to_change] == 0 ) {
			cout << "that is unexpected in file " << __FILE__ << " at line " << __LINE__  << endl;
			j--;
			continue;
		}
		if ( column_types[to_change] == 1 ) {
			if (r_gen()<=0.5)temp = 0;
			else temp = 1;

			#ifdef read_from_file
				s_cre << " uv: " << temp;
				temp = new_v;
			#endif

			if ( double_inequality ( rounded_solution[to_change],temp, epInt ) ) 			{
				rounded_solution[to_change] = temp;
				binary_changed++;
			}
		}
		if (column_types[to_change] == 2){
			double lblb = max ( original_lbs[to_change], starting_solution[to_change] - VIRTUAL_BOUND );
			double ubub = min ( original_ubs[to_change], starting_solution[to_change] + VIRTUAL_BOUND );
			temp = floor ( r_gen() * ( ubub - lblb + 1) );

			#ifdef read_from_file
				s_cre << " uv: " << temp;


// 				cout << int_index_to_change << " " << to_change << " " << temp << " " << new_v <<
				temp = new_v;
			#endif



			if (double_inequality(temp, rounded_solution[to_change], epInt)){
				rounded_solution[to_change] = temp;
				int_change++;
			}
		}

#ifdef read_from_file

	ree2 << rounded_solution[to_change] <<"\t" << binary_changed << " + " <<  int_change <<endl;

		s_cre << " ch: " << binary_changed << "+ " <<  int_change;

		reec << s_cre.str() << endl;
		reer << s_read.str() << endl;
#endif

	}
	
// 	cout << "n_to_change: " << rp << " n_changed: " << binary_changed << " + " << int_change << " = " << binary_changed+int_change << endl;


	return binary_changed+int_change;


}



int cpfp_executable::aggressive_restart ( My_solution& rounded_solution )
{
	vector<double> rs;
	rs.resize(ncols);
	rounded_solution.clear();
	double lblb,ubub,temp,temp2;
	int status;
	int n_changed =0;
	for (unsigned i = 0; i < integer_columns_in_slave.size();++i){
		int t= integer_columns_in_slave[i];
		lblb = max ( original_lbs[t], -VIRTUAL_BOUND );
		ubub = min ( original_ubs[t], + VIRTUAL_BOUND );
		temp = floor ( r_gen() * ( ubub - lblb + 1) );
		temp2 = rounded_solution.get_element(t,status);
		if (double_inequality(temp,temp2,epInt)) {
			rounded_solution.set_element(t,temp);
			n_changed++;
		}
		
	}
	return n_changed;
	
}

int cpfp_executable::aggressive_flip ( My_solution& rounded_solution, double flip_probability )
{
	vector<double> rs;
	rs.resize(ncols);
	rounded_solution.clear();
	double lblb,ubub,temp,temp2;
	int status;
	int n_changed =0;
	for (unsigned i = 0; i < integer_columns_in_slave.size();++i){
		int t= integer_columns_in_slave[i];
		if (column_types[t]==1){
			if (r_gen() < flip_probability){
				if (double_equality(rounded_solution.get_element(t,status), original_lbs[t],epInt)){
					rounded_solution.set_element(t,original_ubs[t]);
				}
				else{
					rounded_solution.set_element(t,original_lbs[t]);
				}
				n_changed++;
			}
		}
		if (column_types[t]==2){
			lblb = max ( original_lbs[t], -VIRTUAL_BOUND );
			ubub = min ( original_ubs[t], + VIRTUAL_BOUND );
			temp = floor ( r_gen() * ( ubub - lblb + 1) );
			temp2 = rounded_solution.get_element(t,status);
			if (double_inequality(temp,temp2,epInt)) {
				rounded_solution.set_element(t,temp);
				n_changed++;
			}
		}
		
	}
	return n_changed;
	
}


int cpfp_executable::f_type_stage_1_set_obj ( vector< double >& rounded_solution ){
    /** set objective for only binary variables */
    for ( int i = 0; i < ncols; ++i ) {
        switch ( column_types[i] ) {
            case 0: /** continious */
            case 2: /** general integer*/

                mip->setObjCoeff ( i,0 );/** */
                break;
            case 1:
                if ( double_equality ( rounded_solution[i],0 ) ) { /** minimize x */
                    mip->setObjCoeff ( i, 1 );
                } else {
                    mip->setObjCoeff ( i, -1 );
                }
                break;
        }
    }
    return 0;


}

My_objective_function cpfp_executable::create_FP_dist_obj_for_solution_stage_1 ( My_solution& my_sol_rounded_solution ){
    vector<double> rounded_s = my_sol_rounded_solution.get_solution_vector();

    //     vector<int> indices = my_sol_rounded_solution.get_indices();
    //     vector<double> elements = my_sol_rounded_solution.get_elements();
    vector<double> objective_s;
    objective_s.resize(ncols,0);
    for ( int i = 0; i < ncols; ++i ) {
        switch ( column_types[i] ) {
            case 0: /** continious */
            case 2: /** general integer*/
// 				mip->setObjCoeff ( i,0 );/** */
				break;
            case 1:

                if ( double_equality ( rounded_s[i],0 ) ) { /** minimize x */
                    objective_s[i] = 1 ;
                } else {
                    objective_s[i] = -1 ;
                }
                break;
        }
    }

    My_objective_function retval(objective_s);

    return retval;

}


My_objective_function cpfp_executable::f_type_stage_1_create_obj_for_FP ( My_solution& my_sol_rounded_solution, double alpha )
{
// 	cout << "STAGE 1: XXXXXXXXXX" <<" line: " << __LINE__ << " alpha: "<< alpha<< endl;
	static int zikkim2 = 0;
	zikkim2 ++;

	My_objective_function distance_obj = create_FP_dist_obj_for_solution_stage_1(my_sol_rounded_solution);

	if (distance_obj.is_empty()) cout << " empty distance_obj at line " << __LINE__ <<endl;

    if (double_equality(alpha,0)){
		return distance_obj;
    }

    if (double_equality(alpha,1)){
		return original_objective_function;
	}

    double distance_obj_multiplier = (1-alpha);

    double other_multiplier = alpha * sqrt(distance_obj.get_size())/original_objective_function.calculate_norm();;



    My_objective_function retval = My_objective_function_multiply_sum(distance_obj_multiplier,distance_obj,other_multiplier,original_objective_function);






// 	ofstream mm("FP_objectives_stage_1.txt",ios_base::app);
// 	mm << "call " << zikkim2 << " alpha: " << alpha << " ";
// 	if (retval.is_empty()) mm << " retval_empty";
// 	else mm << " retval not empty";
// 	retval.shortlineprint(mm);


// 	mm << "---- " << zikkim2 << " -----: " << alpha << " ";
// 	distance_obj.shortlineprint(mm);
// 	mm << "---- " << zikkim2 << " -----: " << alpha << " ";
// 	original_objective_function.shortlineprint(mm);
// 	cout << "STAGE 1: XXXXXXXXXX" <<" line: " << __LINE__ << endl;

//     My_objective_function retval(objective_s);
    return retval;


}


My_objective_function cpfp_executable::create_FP_dist_obj_for_solution_stage_2 ( My_solution& my_sol_rounded_solution )
{
    vector<double> rounded_s = my_sol_rounded_solution.get_solution_vector();

    vector<double> objective_s;
    objective_s.resize(ncols+ fp_status.n_cols_added,0);
    unsigned j =0;
    for ( int i = 0; i < ncols; ++i ) {
        switch ( column_types[i] ) {
            case 0: /** continious */
                //                 mip->setObjCoeff ( i,0 );/** */
                break;
            case 1:
                if ( double_equality ( rounded_s[i],0 ) ) { /** minimize x */
                    objective_s[i] = 1 ;
                } else {
                    objective_s[i] = -1 ;
                }
                break;
            case 2:
//                 mip->setObjCoeff ( i, 0 );

                objective_s[fp_status.added_col_indices[j]]= 1;
                j++;
                break;
        }
    }
    My_objective_function retval(objective_s);
    return retval;

}


My_objective_function cpfp_executable::f_type_stage_2_create_obj_for_FP ( My_solution& my_sol_rounded_solution, double alpha )
{

	static int zikkim = 0;
	zikkim ++;

// 	cout << "SLAVE: XXXXXXXXXX" <<" line: " << __LINE__ << endl;
    My_objective_function distance_obj = create_FP_dist_obj_for_solution_stage_2(my_sol_rounded_solution);

	if (double_equality(alpha,0)){
		return distance_obj;
	}

	if (double_equality(alpha,1)){
		return original_objective_function;
	}

	if (distance_obj.is_empty()) cout << " empty distance_obj at line " << __LINE__ <<endl;

	double distance_obj_multiplier = (1-alpha);

    double other_multiplier = alpha * sqrt(distance_obj.get_size())/original_objective_function.calculate_norm();;



    My_objective_function retval = My_objective_function_multiply_sum(distance_obj_multiplier,distance_obj,other_multiplier,original_objective_function);





// 	ofstream mm("FP_objectives_stage_2.txt",ios_base::app);
// 	mm << "call " << zikkim << " alpha: " << alpha << " ";
//
// 	if (retval.is_empty()) mm << " retval_empty";
// 	else mm << " retval not empty";
//
// 	retval.shortlineprint(mm);
// 	mm << "---- " << zikkim << " -----: " << alpha << " ";
// 	distance_obj.shortlineprint(mm);
// 	mm << "---- " << zikkim << " -----: " << alpha << " ";
// 	original_objective_function.shortlineprint(mm);
	return retval;

}


My_objective_function cpfp_executable::f_type_stage_1_or_2_create_obj_for_FP ( My_solution& my_sol_rounded_solution, double alpha, int stage )
{
	
//     if (my_sol_rounded_solution.is_empty()){
// 		cout << "current_stage: "<< stage << endl;
// 		my_sol_rounded_solution.shortlineprint(cout);
// 	}
	if (stage == 1) return f_type_stage_1_create_obj_for_FP(my_sol_rounded_solution, alpha);
	if (stage == 2) return f_type_stage_2_create_obj_for_FP(my_sol_rounded_solution, alpha);
	
	return f_type_stage_2_create_obj_for_FP(my_sol_rounded_solution, alpha);
	
}


int cpfp_executable::f_type_stage_2_set_obj ( vector< double >& rounded_solution ){
    int j = 0;
    for ( int i = 0; i < ncols; ++i ) {
        switch ( column_types[i] ) {
            case 0: /** continious */
                mip->setObjCoeff ( i,0 );/** */
                break;
            case 1:
                if ( double_equality ( rounded_solution[i],0 ) ) { /** minimize x */
                    mip->setObjCoeff ( i, 1 );
                } else {
                    mip->setObjCoeff ( i, -1 );
                }
                break;
            case 2:
                mip->setObjCoeff ( i, 0 );
                mip->setObjCoeff ( fp_status.added_col_indices[j],1 );
                j++;
                break;
        }
    }
    for ( unsigned i = 0; i < fp_status.general_integer_indices.size(); ++i ) {
        mip->setRowLower ( fp_status.added_row_indices1[i], rounded_solution[fp_status.general_integer_indices[i]] );
        mip->setRowUpper ( fp_status.added_row_indices2[i], rounded_solution[fp_status.general_integer_indices[i]] );
    }
    return 0;
}

int cpfp_executable::f_type_stage_2_set_bounds(My_solution &rd){
	int status; double elem;



	for ( unsigned i = 0; i < fp_status.general_integer_indices.size(); ++i ) {
		elem = rd.get_element(fp_status.general_integer_indices[i],status);
		if(status) cout << "element out of bound. line: "<<__LINE__ <<endl;

		mip->setRowLower ( fp_status.added_row_indices1[i], elem );
		mip->setRowUpper( fp_status.added_row_indices1[i], mip->getInfinity() );

		mip->setRowUpper ( fp_status.added_row_indices2[i], elem );
		mip->setRowLower( fp_status.added_row_indices2[i], -1* mip->getInfinity() );

	}
	return 0;
}


int cpfp_executable::f_type_stage_1_set_obj ( My_solution& my_sol_rounded_solution ){
    vector<double> rounded_s = my_sol_rounded_solution.get_solution_vector();
    return f_type_stage_1_set_obj(rounded_s);
}

int cpfp_executable::f_type_stage_2_set_obj ( My_solution& my_sol_rounded_solution ){
    vector<double> rounded_s = my_sol_rounded_solution.get_solution_vector();
    return f_type_stage_2_set_obj(rounded_s);
}

int cpfp_executable::mip_set_FP_objective_stage_1 ( My_solution& my_sol_rounded_solution, double alpha )
{
    My_objective_function fp_obj =  f_type_stage_1_create_obj_for_FP(my_sol_rounded_solution, alpha);

    set_FP_objective_function(fp_obj);

    return 0;
}

int cpfp_executable::mip_set_FP_objective_stage_2 ( My_solution& my_sol_rounded_solution, double alpha )
{
    My_objective_function fp_obj =  f_type_stage_2_create_obj_for_FP(my_sol_rounded_solution, alpha);

    set_FP_objective_function(fp_obj);


    set_bounds_FP_stage_2(my_sol_rounded_solution);

    return 0;
}

int cpfp_executable::set_bounds_FP_stage_2 ( My_solution& my_sol_rounded_solution )
{
    vector<double> rounded_solution = my_sol_rounded_solution.get_solution_vector();
    for ( unsigned i = 0; i < fp_status.general_integer_indices.size(); ++i ) {
        mip->setRowLower ( fp_status.added_row_indices1[i], rounded_solution[fp_status.general_integer_indices[i]] );
        mip->setRowUpper ( fp_status.added_row_indices2[i], rounded_solution[fp_status.general_integer_indices[i]] );
    }
    return 0;
}


int cpfp_executable::unset_bounds_FP_stage_2 ()
{

    double maxub = mip->getInfinity();
// 	double minlb = -maxub;



    for ( unsigned i = 0; i < fp_status.general_integer_indices.size(); ++i ) {
		double max_change = abs(original_lbs[fp_status.general_integer_indices[i]]) + abs(original_ubs[fp_status.general_integer_indices[i]]);
		max_change = maxub;
        mip->setRowLower ( fp_status.added_row_indices1[i], -max_change );
//         mip->setRowUpper ( fp_status.added_row_indices1[i], max_change );


//         mip->setRowLower ( fp_status.added_row_indices2[i], -max_change );
        mip->setRowUpper ( fp_status.added_row_indices2[i], max_change );
    }
    return 0;
}

int cpfp_executable::get_number_of_general_integers() const{
    return fp_status.general_integer_indices.size();
}

int cpfp_executable::get_number_of_binaries() const{
    return fp_status.binary_indices.size();
}

int cpfp_executable::FP_history_check ( My_solution& rounded_solution, double alpha, int current_stage )
{

// 	get_integer_infeasibility();

// 	fp_hist.int_frac_list.push_back ( n_frac );
// 	fp_hist.double_frac_list.push_back ( d_fractionality );
// 	fp_hist.solution_history.push_back ( starting_solution );


// 	if (fp_hist.compare_last_with_previous()) ;

	return false;
}

int cpfp_executable::f_type_restart_1_or_2 ( My_solution& starting_solution, My_solution& rounded_solution, My_solution& previous_rounded, int stage,int iter ){
// 	bool fcon = true;


	if (stage == 1) return f_type_restart_1(starting_solution,rounded_solution,previous_rounded);
	if (stage ==2) return f_type_restart_2(starting_solution,rounded_solution,previous_rounded,iter);



	return -1;



}

int cpfp_executable::stage_decider(int previous_stage){

	if (previous_stage == 0) {
		cout << " previous_stage : " << previous_stage << endl;
		return -1;
	}

	if (previous_stage == 1) {
		if(fp_inf.itercount++ >= fp_inf.iterlim){
			fp_inf.itercount =0;
			if(fp_inf.resetcount++>= fp_inf.resetlim){
				fp_inf.resetcount = 0;
				return 2;
			}
			return 1;
		}
	}

	if (previous_stage == 2) {
		if(fp_inf.stage2_itercount++ >= fp_inf.stage2_iterlim){
			fp_inf.stage2_itercount =0;
			if(fp_inf.stage2_resetcount++>= fp_inf.stage2_resetlim){
				fp_inf.stage2_resetcount = 0;
				return 1;
			}
			return 2;
		}
	}

	return previous_stage;

}

void cpfp_executable::clear_FP_history(unsigned int keep_latest_n){
	fp_hist.clear(keep_latest_n);
}


int  cpfp_executable::add_to_history ( int n_frac, double d_frac, My_solution sol ){
	fp_hist.int_frac_list.push_back(n_frac);
	fp_hist.double_frac_list.push_back(d_frac);
	fp_hist.solution_history.push_back(sol);
// 	if (fp_hist.solution_history.size() >2000){
// 		fp_hist.solution_history.
// 	}
	return fp_hist.solution_history.size();
}

int cpfp_executable::compare_last_with_previous ( bool use_hist )
{
	return fp_hist.compare_last_with_previous_v1 ( use_hist );
}


int cpfp_executable::get_history (My_solution &z, int n )
{
	if (n > (int)fp_hist.solution_history.size()) return -1;
	if (n < 0) return -2;

	z = fp_hist.solution_history[fp_hist.solution_history.size() - n ];
	return 0;
}

void cpfp_executable::print_solution_wth_column_types (My_solution& sol, ostream& out )
{
	vector<int> indi = sol.get_indices();
	vector<double> elem = sol.get_elements();


	for (unsigned i = 0; i <indi.size(); ++i){

			out << ", X_("<< column_types[indi[i]]<<")" << indi[i] << "= "<< elem[i];
		}
		out << endl;
	return;
}

void cpfp_executable::print_bounds()
{

	cout << "doubleINF: "<<  mip->getInfinity() << endl;
	const double *ubs = mip->getColUpper();
	const double *lbs = mip->getColLower();
	for (int i = 0; i < ncols;++i){
		cout << lbs[i]<< " "		;
	}
	cout<< endl;
	for (int i = 0; i < ncols;++i){
		cout << ubs[i]<< " "		;
	}
	cout<< endl;

}

double cpfp_executable::calculate_pseudo_bound()
{
// 	static double pseudo_bound = -1;

// 	if (pseudo_bound >0) return pseudo_bound;

	double max_ub = 0;
	double max_lb = 0;
// 	const double *ubs = mip->getColUpper();
// 	const double *lbs = mip->getColLower();

	int nrows  = mip->getNumRows();

	const double *row_ubs = mip->getRowUpper();
	const double *row_lbs = mip->getRowLower();

	for (int i = 0; i < nrows;++i){
		double ub = abs(row_ubs[i]);
		double lb = abs(row_lbs[i]);
		if (ub >= mip->getInfinity()/10){
			if (max_ub < ub) max_ub = ub;
		}
		if (lb >= mip->getInfinity()/10){
			if (max_lb < lb) max_lb = lb;
		}
	}

	double rhs = 1;
	if (max_lb > rhs) rhs = max_lb;
	if (max_ub > rhs) rhs = max_ub;

	double sp;
	/** version 2*/
// 	sp = mip->getNumElements()/(ncols*nrows);

// 	if (sp <1)
	sp =1;

	pseudo_bound = 1.0e12;


	double temp = rhs * ncols * sp;
	if (temp < pseudo_bound) pseudo_bound = temp;

	return pseudo_bound;





}



int cpfp_executable::set_pseudo_bounds(double sb){

	const double *ubs = mip->getColUpper();
	const double *lbs = mip->getColLower();

	for (int i = 0; i < ncols;++i){
		if (ubs[i] > sb){
			mip->setColUpper(i,sb);
			extra_mip->setColUpper(i,sb);
			mipCP->setColUpper(i,sb);
			original_ubs[i] = sb;

// 			if (ubs[i]<0) mip->setColUpper(i,-pseudo_bound);
		}
		if (lbs[i] < -sb){
			mip->setColLower(i,-sb);
			extra_mip->setColLower(i,-sb);
			mipCP->setColLower(i,-sb);

			original_lbs[i] = -sb;

		}
	}

	return 0;
}

int cpfp_executable::get_binary_fractionality_of_solution( My_solution &msolution, double& double_infeasibility ){
	vector<int>    indi = msolution.get_indices();
	vector<double> elem = msolution.get_elements();

	double_infeasibility =0;
	int retval =0;
	double roundsol = 0;
	double frac     = 0;
	for (unsigned i = 0 ; i < indi.size(); ++i){
		if ((0<=indi[i]) && (indi[i]<ncols)){
			switch(column_types[indi[i]]){
				case 0:
				case 2:
					continue;
					break;
				case 1:
                    roundsol = floor (elem[i] + 0.5);
					frac     = fabs ( roundsol - elem[i]);
					if (double_equality(frac,0,1e-6) ){
						continue;
					}
					else{
						retval++;
						double_infeasibility += frac;
					}
					break;
				default:
					cout << "should not happen column type :" << column_types[indi[i]] << endl;
					break;
			}
		}
		else{
			cout << "should not happen " << indi[i] << endl; /** ERROR RETURN */
		}

	}

	return retval;

}

int cpfp_executable::get_fractionality_of_solution ( My_solution &msolution, double& double_infeasibility ){
	vector<int>    indi = msolution.get_indices();
	vector<double> elem = msolution.get_elements();

	double_infeasibility =0;
	int retval =0;
	double roundsol = 0;
	double frac     = 0;
	for (unsigned i = 0 ; i < indi.size(); ++i){
		if ((0<=indi[i]) && (indi[i]<ncols)){
			switch(column_types[indi[i]]){
				case 0:
					continue;
					break;
				case 2:

				case 1:

					roundsol = floor (elem[i] + 0.5);
					frac     = fabs ( roundsol - elem[i]);
					if (double_equality(frac,0,1e-6) ){
						continue;
					}
					else{
						retval++;
						double_infeasibility += frac;
					}
					break;
				default:
					cout << "should not happen column type :" << column_types[indi[i]] << endl;
					break;
			}
		}
		else{
			cout << "should not happen " << indi[i] << endl; /** ERROR RETURN */
		}

	}

	return retval;

}


int cpfp_executable::polish_integer_solution ( My_solution &sol, My_solution &new_sol, double& new_obj ){


    /** TODO: rewrite this part in a more efficient way*/
// 	vector<int>    indi = sol.get_indices();
// 	vector<double> elem = sol.get_elements();
//
	vector<double> solvec = sol.get_solution_vector();
// 	int new_ncols = mip->getNumCols();


	double d_frac;
	int i_frac;
	double val;

//	const double *lbs =  mip->getColLower();
//	const double *ubs =  mip->getColUpper();

//	vector<double> vec_lb;
//	vector<double> vec_ub;
//
//	vec_lb.resize(ncols,0);
//	vec_ub.resize(ncols,0);


// 	cout << "new_ncols:" << new_ncols << endl;
// 	for (int i =0; i< new_ncols;++i){
// 		cout<<i<<":"<<lbs[i]<< "_"<< ubs[i] <<"  ";
// 	}
// 	cout << endl;

int retval = -9;
	for (int i = 0 ; i < ncols; ++i){
// 		vec_lb[i] = lbs[i];
// 		vec_ub[i] = ubs[i];
		switch(column_types[i]){
			case 0:
				break;
			case 1:
			case 2:
				/** fix  */
				val = solvec[i];
//				vec_lb[i] = lbs[i];
//				vec_ub[i] = ubs[i];
				mip->setColBounds(i,val,val);
				break;
			default:
				break;
		};
	}



	set_objective_function(original_objective_function);

	int opt_sol = optimize();
		if (opt_sol!= 1) {

		retval = opt_sol;
	}
	else{
		/**  */
		i_frac = get_solution(new_sol,d_frac);
		/** note: it should be i_frac = 0, d_frac = 0 */
		if ((i_frac > 0 ) || (double_inequality(d_frac,0,epInt))){
			/** that should not happen */
			cout << "LINE: " << __LINE__ << " that should not happen " << endl;
			retval =0;
		}
		else{
           // new_sol.print();
            new_obj = new_sol.get_original_objective();
//			calculate_original_objective_value_of_solution(new_obj,new_sol);
//			new_sol.set_original_objective(new_obj);
			best_incumbent = new_sol;
			retval = 1;
		}
	}

    set_objective_function(auxilary_objective_function);



// 	mip->setColSetBounds();
// 		switch(column_types[i]){
//
// 			case 0:
// 				break;
// 			case 1:
// 			case 2:
// 				/** fix  */
// 				// 				val = solvec[i];
// 				mip->setColBounds(i,lbs[i],ubs[i]);
// 				break;
// 			default:
// 				break;
// 		};

	for (int i = 0 ; i < ncols; ++i){
		switch(column_types[i]){
			case 0:
				break;
			case 1:
			case 2:
				/** fix  */
// 				val = solvec[i];
// 				vec_lb[i] = lbs[i];
// 				vec_ub[i] = ubs[i];
				mip->setColBounds(i,original_lbs[i],original_ubs[i]);
				break;
			default:
				break;
		};
	}



// 	cout << "new_ncols:" << new_ncols << endl;
// 	for (int i =0; i< new_ncols;++i){
// 		cout<<i<<":"<<lbs[i]<< "_"<< ubs[i] <<"  ";
// 	}
// 	cout << endl;

// 	cout << " polish_integer_solution completed" <<endl;



	return retval;
}

int cpfp_executable::get_FP_INT_PARAM(FP_INT_PARAM_SET param, int &arg) const{
    if (param < ENM_FP_INT_PARAM_BGN || param > ENM_FP_INT_PARAM_END) return -1;
    arg = FP_INT_PARAM[param];
    return 0;
}

void cpfp_executable::print_FP_parameters ( ostream& out ) const
{
	for (int s = (int)ENM_FP_INT_PARAM_BGN; s < (int)ENM_FP_INT_PARAM_END; ++s)
		out << "FP_INT_PARAM["<<s<<"]: " << FP_INT_PARAM[s] <<endl;
	for (int s = (int)ENM_FP_DOUBLE_PARAM_BGN; s < (int)ENM_FP_DOUBLE_PARAM_END; ++s)
		out << "FP_DBL_PARAM["<<s<<"]: " << FP_DBL_PARAM[s] <<endl;
	
	return;
}

int cpfp_executable::get_FP_DOUBLE_PARAM(FP_DOUBLE_PARAM_SET param, double &arg) const{
	if (param < ENM_FP_DOUBLE_PARAM_BGN || param > ENM_FP_DOUBLE_PARAM_END) return -1;
	arg = FP_DBL_PARAM[param];
	return 0;
}

int cpfp_executable::set_FP_INT_PARAM(FP_INT_PARAM_SET param, int arg) {
    if (param < ENM_FP_INT_PARAM_BGN || param > ENM_FP_INT_PARAM_END) return -1;
    FP_INT_PARAM[param] = arg;
    return 0;
}

int cpfp_executable::set_FP_DOUBLE_PARAM(FP_DOUBLE_PARAM_SET param, double arg) {
    if (param < ENM_FP_DOUBLE_PARAM_BGN || param > ENM_FP_DOUBLE_PARAM_END) return -1;
    FP_DBL_PARAM[param] = arg;
    return 0;
}

int cpfp_executable::getNextIntegerPoint ( My_solution& cont_sol, My_solution& int_sol, int stage )
{
	int minflip = FP_INT_PARAM [ENM_FP_INT_PARAM_MINFLIP];
	int n_changed = 0;
	priority_queue<dipair, vector<dipair>, dipairCmp> q;
	int to_be_changed = (minflip +(int)floor(r_gen() * 2 * minflip))/2;

// 	to_be_changed = 10000;
	double sigma, min_sigma = FP_DBL_PARAM[ENM_FP_DOUBLE_PARAM_FLIP_TRESHOLD];
//     cout << " min_sigma " << min_sigma <<endl;
	
// 	usleep(100000);
// 	int n_bin_changed = 0;
// 	int n_int_changed = 0;
	int status1,status2;
	for ( unsigned j = 0; j < integer_columns_in_slave.size(); ++j ) {
		int i = integer_columns_in_slave[j];
		switch (column_types[i]){
			case 0:
				cout << " THIS SHOULD NOT HAPPEN ANY MORE line: " << __LINE__ << endl;
				break;
			case 2:
				if (stage == 1) break;
			case 1:
				sigma = abs(cont_sol.get_element(i,status1) - int_sol.get_element(i,status2));
				if ((sigma > min_sigma) || double_equality(min_sigma,0)){
					q.push(dipair(sigma, j));
					if( ((int)q.size()) > to_be_changed)
					{
						q.pop();
						min_sigma=q.top().first;
					}
				}
				break;

		}


	}


// 	int q_size = q.size();
// 	cout << "---------------------" << endl;

	while( !q.empty() )
	{
		int j = q.top().second;
		int i = integer_columns_in_slave[j];
// 		cout << j << " " << int_sol.get_element(i,status1)  << " ";

		int status1,status2;
		if( int_sol.get_element(i,status1) > cont_sol.get_element(i,status2) )
			int_sol.element_mm(i); // roundedIntVars[i]--;
		else
			int_sol.element_pp(i); // roundedIntVars[i]++;

// 		cout << int_sol.get_element(i,status1) <<" "<< q.top().first << " "<<cont_sol.get_element(i,status2)  << endl;

		n_changed++;
		q.pop();
	}
// 	static int is = 0;
// 	is ++;
// 	cout << "---------------------------------- here this many times " <<  is<< " toBeChanged: " << to_be_changed << " changed:" << n_changed << " q_size: " << q_size<< endl;


// 	usleep(1000000);

	return n_changed;



}

int cpfp_executable::getDivergentPoint(My_solution &cont_sol, My_solution &int_sol, int stage, int slave_id){
    
//     if (slave_id <= 0) return 0;
    
	int reduced_slave_id = slave_id;
    int n_changed = 0;
    
    
    priority_queue<dipair, vector<dipair>, dipairCmp> q;
   
     
    
    int to_be_changed = ceil(log2(slave_id+1));
    
    //  to_be_changed = 10000;
    double sigma, min_sigma = 0;
    //     cout << " min_sigma " << min_sigma <<endl;
    
    //  usleep(100000);
    //  int n_bin_changed = 0;
    //  int n_int_changed = 0;
    int status1,status2;
    for ( unsigned j = 0; j < integer_columns_in_slave.size(); ++j ) {
        int i = integer_columns_in_slave[j];
        switch (column_types[i]){
            case 0:
                cout << " THIS SHOULD NOT HAPPEN ANY MORE line: " << __LINE__ << endl;
                break;
            case 2:
                if (stage == 1) break;
            case 1:
                sigma = abs(cont_sol.get_element(i,status1) - int_sol.get_element(i,status2));
                if ((sigma > min_sigma) || double_equality(min_sigma,0)){
                    q.push(dipair(sigma, j));
                    if( ((int)q.size()) > to_be_changed)
                    {
                        q.pop();
                        min_sigma=q.top().first;
                    }
                }
                break;
        }
    }
    
    
    
    if (q.size() < to_be_changed){
        
        cout << "less variables to search than available " <<endl;
    }
    
//     bool while_l = true;
//     int flip_index = 1;
    
// 	stringstream cc;
// 	cc << "slave " << slave_id << " "; 
   	
   	for (int k = 0; k < to_be_changed;++k){
		int j = q.top().second;
		int i = integer_columns_in_slave[j];
		
        if (reduced_slave_id%2 == 1){
            int status1,status2;
            if( int_sol.get_element(i,status1) > cont_sol.get_element(i,status2) )
                int_sol.element_mm(i); // roundedIntVars[i]--;
            else
                int_sol.element_pp(i); // roundedIntVars[i]++;
            
			n_changed++;
// 			cc << 1; 
        }
//        	else cc << 0;
		
       	reduced_slave_id/=2;
       
		q.pop();
    }
    
//     cout << cc.str().c_str() <<endl;
    
//     usleep(1000000);
    
    return n_changed;
    
    
}



vector< int > cpfp_executable::get_column_types() const
{
	return column_types;
}

int cpfp_executable::set_bounds_extra_mip ( My_solution& s, int stage )
{
	vector<int> ind = s.get_indices();

	unsigned z = 0;
	int s_ind = ncols +1;
	if (ind.size()>0) s_ind = ind[z];
	int status;

	for (int i = 0;i<ncols;++i){

		switch(column_types[i]){
			case 0:
				continue;
				break;
			case 2:
				if (stage == 1){
					continue;
					break;
				}
			case 1:
				if (i < s_ind){ /** set to original bounds*/
					extra_mip->setColLower(i,original_lbs[i]);
					extra_mip->setColUpper(i,original_ubs[i]);
				}
				else {
					if (i == s_ind){
						extra_mip->setColLower(i,s.get_element(i,status));
						extra_mip->setColUpper(i,s.get_element(i,status));
						z++ ;
						if (z<ind.size()) s_ind = ind[z];
						else s_ind = ncols+1;
					}
					else { /** i > s_ind*/
						z++;
						if (z<ind.size()) s_ind = ind[z];
						else s_ind = ncols+1;

					}
				}
				break;
		}

	}

	return 0;
// 	mip->setColLower();

}

int cpfp_executable::optimize_extra_mip( My_solution& s, int stage )
{
	set_bounds_extra_mip(s,stage);

	if (!initial_solved_extra){
		extra_mip->initialSolve();
		initial_solved_extra = true;
	}
	else{
		extra_mip->resolve();
	}
	OPT_RESULT(extra_mip)
	
}

My_solution cpfp_executable::get_solution_of_extra_mip ( int stage ){

	const double *sol = extra_mip->getColSolution();

	My_solution s(ncols,sol);
	return s;
}

double cpfp_executable::get_lambda_start ( My_solution& s_sol, My_solution& t_sol, int index){
	
	
	int status_s,status_t;
	double x_i_s, x_i_t;
	x_i_s = s_sol.get_element(index,status_s);
	x_i_t = t_sol.get_element(index,status_t);
	assert( (status_t + status_s) == 0);
	
	double plus_minus_point_five = 0.5;
	if (x_i_s < x_i_t) plus_minus_point_five =-0.5;
	
	return ((double_round(x_i_s) + plus_minus_point_five - x_i_s) / (x_i_t - x_i_s));
}

double cpfp_executable::get_lambda_min ( My_solution& s_sol, My_solution& t_sol, int index, double &abs_one_over_dif){

	int status_s,status_t;
	double x_i_s, x_i_t;
	x_i_s = s_sol.get_element(index,status_s);
	x_i_t = t_sol.get_element(index,status_t);
	assert( (status_t + status_s) == 0);
	double plus_minus_point_five = 0.5;
	double dif = x_i_t - x_i_s;
	
	if (dif<0) plus_minus_point_five =-0.5;
	
	abs_one_over_dif = abs(1/dif);
	double lambda_start =  ((double_round(x_i_s) + plus_minus_point_five - x_i_s) / (dif));
	double lambda_min = 0 ; 
	
	if (lambda_start>1) lambda_min = 1.1;
	else{
		if (lambda_start>0)	lambda_min = lambda_start - floor(lambda_start / abs_one_over_dif) * abs_one_over_dif; 
		else lambda_min = lambda_start +  ceil(-lambda_start/abs_one_over_dif)*abs_one_over_dif;
		
	}
	if (lambda_start<=1 && (abs_one_over_dif<1)){
		if (double_inequality(lambda_start,lambda_min))
		cout << "index " << index <<"\t" << "x_i_s: " << x_i_s << "\tx_i_t: " << x_i_t << "\tlam_s: " << lambda_start 
			<< "\tdif: " << dif << "\t1/dif: " << 1/dif << "\tlam_m: " << lambda_min << endl;
// 			else cout << ".";
// 			if (double_inequality(lambda_start,lambda_min)) cout << " XXXXXXXX" <<endl;
// 			else cout << endl;
			
	}
	
	return lambda_min;
	
	
}

int cpfp_executable::FP_line_search_fill_lambda_dif ( My_solution& s_sol, My_solution& t_sol, double alpha_min, double alpha_max )
{
    
	My_solution ss = s_sol;
	My_solution tt = t_sol;
	
	if (double_inequality(alpha_min,0) || double_inequality(alpha_max,1)){
		vector<double> s_v = s_sol.get_solution_vector();
		vector<double> t_v = t_sol.get_solution_vector();
		vector<double>dif = vector_subtract<double>(t_v, s_v);
		if (double_inequality(alpha_min,0)){
			ss = My_solution(vector_sum<double>(s_v,vector_multiply<double>(alpha_min,dif))) ;
		}
		if (double_inequality(alpha_max,1)) {
			tt = My_solution(vector_sum<double>(s_v, vector_multiply<double>(alpha_max,dif)));
		}
	}

	
	for (unsigned i=0;i < integer_columns_in_slave.size(); ++i){
		int status_s,status_t;
		double x_i_s, x_i_t;
		x_i_s = ss.get_element(integer_columns_in_slave[i],status_s);
		x_i_t = tt.get_element(integer_columns_in_slave[i],status_t);
		
// 		mip->get_ub();
		
		assert( (status_t + status_s) == 0);
		double plus_minus_point_five = 0.5;
		double dif = x_i_t - x_i_s;
		if (double_equality(dif,0)) continue;
		
		
		double lambda_start;
		
		double ub_i = original_ubs[integer_columns_in_slave[i]];
		double lb_i = original_lbs[integer_columns_in_slave[i]];
		
// 		assert (ub_i>=x_i_s);
		
		
// 		if (lb_i > x_i_s ) cout << "lb_i "<< lb_i <<" x_i_s "<< x_i_s <<endl;
// 		cout.flush();
// 		assert (lb_i<=x_i_s);
		
		if (dif<0) {
			lambda_start = (double_round(get_max<double,double>(x_i_s, lb_i)) + 0.5 - x_i_s)/(dif);
		}
		else {
			lambda_start = (double_round(get_min<double,double>(x_i_s, ub_i)) - 0.5 - x_i_s)/(dif);
			plus_minus_point_five = -0.5;
		}
		
/*		
		double lambda_start2 =  ((double_round(x_i_s) + plus_minus_point_five - x_i_s) / (dif));
		
		if (double_inequality(lambda_start,lambda_start2)){
			cout << "lb_i "<< lb_i <<" x_i_s "<< x_i_s <<endl;
			cout << "ub_i "<< ub_i <<" dif "<< dif <<endl;
			
			
			usleep(5000000);
		}
		*/
		
		
		
		if (lambda_start >1) continue;
		
// 		if (lambda_start )

		/** we will add this integer with dif d_k and lambda*/
		
		double abs_one_over_dif = abs(1/dif);
		
		double lambda_min = 0; 
		
		if (lambda_start>0)	lambda_min = lambda_start - floor(lambda_start / abs_one_over_dif) * abs_one_over_dif; 
		else lambda_min = lambda_start +  ceil(-lambda_start/abs_one_over_dif)*abs_one_over_dif;

		if (lambda_min > 1) continue;
		
		/** at the bound and to the infeasible direction  */
		int d_k = 1;
		
		if (dif<0) d_k = -1;
		
		if ( (d_k == -1) && (double_equality(double_round(x_i_s), lb_i))) continue;
		if ( (d_k == 1) && (double_equality(double_round(x_i_s), ub_i))) continue;
		
		
		LS_indices.push_back(integer_columns_in_slave[i]);
		LS_one_over_dif.push_back(abs_one_over_dif);
		LS_lambda.push_back(lambda_min);
		LS_d_k.push_back (d_k);
		
		if (lambda_start<=1 /*&& (abs_one_over_dif<1)*/){
			if (double_inequality(lambda_start,lambda_min)){}
				
				cout << "index " << integer_columns_in_slave[i] <<"\t" << "x_i_s: " << x_i_s << "\tx_i_t: " << x_i_t << "\tlam_s: " << lambda_start << "\tdif: " << dif << "\t1/dif: " << 1/dif << "\tlam_m: " << lambda_min << endl;
			// 			else cout << ".";
			// 			if (double_inequality(lambda_start,lambda_min)) cout << " XXXXXXXX" <<endl;
			// 			else cout << endl;
			
		}
		
// 		lambda[i]= get_lambda_min(ss,tt,integer_columns_in_slave[i],dif[i]);
	}
	
   
        

    return 1;
}

int cpfp_executable::FP_line_search_get_next_solution ( My_solution& to_be_updated ){
	static bool z = true;
	
	if (z) {
		cout << "ncols "<< ncols <<endl;
		z =false;
	}
	
	cout << "LS_lambda.size() " << LS_lambda.size() <<endl; 
	
	
	usleep(100000); 
	if (LS_lambda.size() == 0) return -1;
	
	
	bool del_now = true;
	int status1;
	while(del_now && (LS_lambda.size()>0)){
		del_now = false;
		bool no_next_step = false;
		
		int min_of_index = 0;
		double min_of_value = LS_lambda[0];
	
		for(unsigned i = 1; i < LS_lambda.size(); ++i){
			if (min_of_value > LS_lambda[i]){
				min_of_value = LS_lambda[i];
				min_of_index = i;
			}
		}
	
		if (LS_d_k[min_of_index] == 1){
			if ( double_equality(to_be_updated.get_element(LS_indices[min_of_index],status1), original_ubs[LS_indices[min_of_index]])){
				del_now = true;
				cout << "index: " << LS_indices[min_of_index] << " ub: " << original_ubs[LS_indices[min_of_index]] <<" val: " << to_be_updated.get_element(LS_indices[min_of_index],status1) << " LS_d_k: " << LS_d_k[min_of_index]  << " lambda: " << LS_lambda[min_of_index] << " 1/dif: " << LS_one_over_dif[min_of_index] << endl; 
			}
			if (!del_now){
				to_be_updated.element_pp(LS_indices[min_of_index]);
				if ( double_equality(to_be_updated.get_element(LS_indices[min_of_index],status1), original_ubs[LS_indices[min_of_index]])){
					no_next_step = true;	
				}
			}
		}
		else {
			if ( double_equality(to_be_updated.get_element(LS_indices[min_of_index],status1), original_lbs[LS_indices[min_of_index]])){
				del_now = true;
				cout << "index: " << LS_indices[min_of_index] << " lb: " << original_lbs[LS_indices[min_of_index]] <<" val: " << to_be_updated.get_element(LS_indices[min_of_index],status1) << " LS_d_k: " << LS_d_k[min_of_index]  << " lambda: " << LS_lambda[min_of_index] << " 1/dif: " << LS_one_over_dif[min_of_index] << endl; 
				
			}
			if (!del_now){
				to_be_updated.element_mm(LS_indices[min_of_index]);
				if ( double_equality(to_be_updated.get_element(LS_indices[min_of_index],status1), original_lbs[LS_indices[min_of_index]])){
					no_next_step = true;	
				}
			}
			
		}
	
		LS_lambda[min_of_index] += LS_one_over_dif[min_of_index];
	
		if (LS_lambda[min_of_index]> 1 || (del_now) || no_next_step){
			LS_lambda.erase(LS_lambda.begin() + min_of_index);
			LS_indices.erase(LS_indices.begin() + min_of_index);
			LS_d_k.erase(LS_d_k.begin() + min_of_index);
			LS_one_over_dif.erase(LS_one_over_dif.begin() + min_of_index);
		}
		
// 		if (del_now) {
// 			cout << "del_now true" <<endl;
// 		}
	}
		
	return LS_lambda.size();
	
	
	return -10;
	
}




