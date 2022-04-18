#include "my_objective_function.h"

using namespace std;



My_objective_function::My_objective_function() {}
My_objective_function::My_objective_function ( vector<double> arg )
{
    original_size = arg.size();
    for ( int i = 0; i < original_size; ++i ) {
        if ( !double_equality ( arg[i],0 ) ) {
            indices.push_back ( i );
            elements.push_back ( arg[i] );
        }
    }
// 	normalize();
    size = elements.size();
}
My_objective_function::My_objective_function ( const My_objective_function &other )
{
    original_size = other.original_size;
    size= other.size;
    elements = other.elements;
    indices = other.indices;

}

vector<double>  My_objective_function::get_elements() const
{
    return elements;
}
void My_objective_function::set_elements ( vector<double>  arg )
{
    elements = arg;
    return;
}

vector<int> My_objective_function::get_indices() const
{
    return indices;
}
void My_objective_function::set_indices ( vector<int>  arg )
{
    indices = arg;
    return;
}

int My_objective_function::get_index ( int whichindex ) const
{
    return indices[whichindex];
}

double My_objective_function::get_element ( int whichindex ) const
{
    return elements[whichindex];
}



int My_objective_function::get_size() const
{
    return size;
}
void My_objective_function::set_size ( int arg )
{
    size = arg;
    return;
}

int My_objective_function::get_original_size() const
{
    return original_size;
}
void My_objective_function::set_original_size ( int arg )
{
    original_size = arg;
    return;
}

int My_objective_function::get_objective_function ( vector<double> &retval, ostream &asd3 ) const
{

    retval.resize ( original_size,0 );

    for ( int i = 0 ; i<size; ++i ) {
        if ( indices[i]>=original_size ) {
            return -1;
        }
        retval[indices[i]] = elements[i];
    }
    return 0;
}

int My_objective_function::get_objective_function ( double *retval, ostream &asd3 ) const
{

// 	retval.resize(original_size);
    for ( int i = 0 ; i<original_size; ++i ) {
        retval[i] = 0;
    }
    for ( int i = 0 ; i<size; ++i ) {
        if ( indices[i]>=original_size ) {
            return -1;
        }
        retval[indices[i]] = elements[i];
    }
    return 0;
}

My_objective_function My_objective_function::operator= ( const My_objective_function &other )
{
    original_size = other.original_size;
    size= other.size;
    elements.resize(size);
    indices.resize(size);
    for(int i = 0; i < size;++i){
        elements[i] = other.elements[i];
        indices[i] = other.indices[i];
    }
    return *this;

}
void My_objective_function::normalize ( std::ostream &out )
{
    double norm =0;
    if ( norm <  1.0e-20 ) {
        norm =0;
        for ( unsigned i =0; i < elements.size(); ++i ) {
            norm += elements[i]*elements[i];
        }
    }
    divide_coefficients ( sqrt ( norm ) );
}

void My_objective_function::divide_coefficients ( double x )
{
    if ( x<0 ) {
        x*=-1;
    }

    for ( unsigned i = 0; i < elements.size(); ++i ) {
        elements[i]/=x;
    }

    return;
}





     
#ifdef USING_MPI
void My_objective_function::_send_mpi ( int destination, int tag, vector<MPI_Request> &v_request, vector<MPI_Status> v_status, MPI_SEND_TYPE SEND_TYPE )
{

    if ( tag != TAG_SENDING_OBJECTIVE_FUNCTION ) {
        if ( tag != TAG_TESTING ) {
            cout << "ERROR: tag: " << tag << " BUT should be TAG_SENDING_OBJECTIVE_FUNCTION: " <<  TAG_SENDING_OBJECTIVE_FUNCTION << endl;
        }
    }
    v_request.resize ( 4,MPI_REQUEST_NULL );
    v_status.resize ( 4 );

    int* indices2 = ( int * ) malloc ( size * sizeof ( int ) );
    double* elements2 = ( double * ) malloc ( size * sizeof ( double ) );
    for ( int i =0; i < size; ++i ) {
        indices2 [i] = indices[i];
        elements2[i] = elements[i];
    }


    COMMON_DATA_SEND_FUNCTION


    free ( indices2 );
    free ( elements2 );
    return;

}

void My_objective_function::_recv_mpi ( int source, int tag, vector<MPI_Request> &v_request, vector<MPI_Status> v_status, MPI_RECV_TYPE RECV_TYPE )
{
    if ( tag != TAG_SENDING_OBJECTIVE_FUNCTION ) {
        if ( tag != TAG_TESTING ) {
            cout << "ERROR: tag: " << tag << " BUT should be TAG_SENDING_OBJECTIVE_FUNCTION: " <<  TAG_SENDING_OBJECTIVE_FUNCTION << endl;
        }
    }
    v_request.resize ( 4,MPI_REQUEST_NULL );
    v_status.resize ( 4 );

    COMMON_DATA_RECV_FUNCTION;

    return;
}

#else
void My_objective_function::_send_alt1(int destination, int tag, my_signaler comm, ALT1_SEND_TYPE SEND_TYPE ){
	TAG_CONTROL(TAG_SENDING_OBJECTIVE_FUNCTION)
	COMMON_DATA_SERIALIZE_FUNCTION
	//we have s as stringstream
	comm.send_object(destination, tag, s.str(),SEND_TYPE);
	return;

}
void My_objective_function::_recv_alt1(int source, int tag, my_signaler comm, ALT1_RECV_TYPE RECV_TYPE ){
    TAG_CONTROL(TAG_SENDING_OBJECTIVE_FUNCTION)

	int rtag  = tag;
	string sr = comm.recv_object(source, rtag, RECV_TYPE);	
	COMMON_DATA_UNSERIALIZE_FUNCTION(sr)
	return;
}

#endif


    
void My_objective_function::print ( std::ostream & out )
{
    out << "Original size: " << original_size << endl;
    out << "size: " << size << endl;
    for ( int i =0; i < size; ++i ) {
        out << " + "<< elements[i] << " X_" << indices[i];
    }
    out << endl;
    return;
}


int My_objective_function::calculate_value ( double& retval, const double* x )
{
    retval = 0;
    for ( int i =0; i< size ; ++i ) {
        if ( indices[i]>original_size ) {
            return -1;
        }
        retval += ( x[indices[i]] * elements[i] );
    }
    return 0;
}


int My_objective_function::calculate_value_of_vector ( double& retval, vector< double >& sol )
{
	retval = 0;
	for ( int i =0; i< size ; ++i ) {
		if ( indices[i]>original_size ) {
			return -1;
		}
		retval += ( sol[indices[i]] * elements[i] );
	}
	return 0;
}


bool My_objective_function::is_empty()
{
    if ( original_size > 0 ) {
        return false;
    }

    return true;
}



My_objective_function My_objective_function::operator+= ( const My_objective_function& other )
{
    
    
    
    
    /** QUICK AND DIRTY SOLUTION */
    vector<double> v1, v2, v3;
    this->get_objective_function ( v1 );
    other.get_objective_function ( v2 );


    v3 = vector_sum<double> ( v1,v2 );

    My_objective_function asd ( v3 );

    *this = asd;

    return *this;


}




double My_objective_function::calculate_norm(){

    double retval =0;
    
    for (unsigned i = 0;i<elements.size();++i){
        retval+=elements[i]*elements[i];
    }
    return sqrt(retval);
    
    
}


My_objective_function My_objective_function_multiply_sum ( double l1, const My_objective_function& lhs, double r2, const My_objective_function& rhs )
{
    
    
    /** BETTER method DONE*/
    
    vector<int> indices_1,indices_2; 
    vector<double> elements_1,elements_2; 
    
    indices_1 = lhs.get_indices();
    indices_2 = rhs.get_indices();
    
    elements_1 = lhs.get_elements();
    elements_2 = rhs.get_elements();
    
    vector<int> indices_final;
    vector<double> elements_final;
    unsigned i=0, j=0;
    
    while((i<indices_1.size()) && (j<indices_2.size())){
        
        if(indices_1[i] == indices_2[j]){
            indices_final.push_back(indices_1[i]);
            elements_final.push_back( l1*elements_1[i]+ r2*elements_2[j]);
            ++i;
            ++j;
        }
        else{ 
            if (indices_1[i] < indices_2[j]) {
                indices_final.push_back(indices_1[i]);
                elements_final.push_back( l1*elements_1[i]);
                ++i;
            }
            else{
                indices_final.push_back(indices_2[j]);
                elements_final.push_back( r2*elements_2[j]);
                ++j;
            }
        }
    }
    
    while((i<indices_1.size()) ){
        indices_final.push_back(indices_1[i]);
        elements_final.push_back( l1*elements_1[i]);
        ++i;
    }
    while((j<indices_2.size())){
        indices_final.push_back(indices_2[j]);
        elements_final.push_back( r2*elements_2[j]);
        ++j;           
    }
    
    My_objective_function retval2;
    retval2.set_elements(elements_final);
    retval2.set_indices(indices_final);
    retval2.set_original_size(lhs.get_original_size());
    retval2.set_size(indices_final.size());
    
    return retval2;
    
    
    /**  */ 
    
    /** QUICK AND DIRTY SOLUTION */
    
    vector<double> v1, v2, v3;
    lhs.get_objective_function ( v1 );
    rhs.get_objective_function ( v2 );
    
    v3 = vector_multiply_sum<double> ( l1,v1,r2,v2 );
    
    
    My_objective_function retval ( v3 );
    
    return retval;
    
    
    
}



void My_objective_function::shortlineprint ( ostream& out )
{
	if (size == 0) out << " size: " << size <<" elements.size(): " << elements.size() << " indices.size(): " << indices.size();
	for( int i =0; i < size; ++i ) {
		out << setw(16)<< "["<< indices[i] << "]="<< elements[i];
	}
	out << endl;
	return;
}






// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
