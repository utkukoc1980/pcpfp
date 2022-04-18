#include "my_solution.h"

using namespace std;

vector<double>  My_solution::get_elements() const {
	return elements;
}


double My_solution::get_element(int index,int &status ) const{
	if (index < 0 || index > original_size ) status = 1;
	status = 0;
	for (unsigned i = 0 ;i < indices.size();++i){
		if (indices[i]==index) return elements[i];
	}
	return 0;

}


void My_solution::set_elements( vector<double>  arg ) {
	elements = arg;
	return;
}

int My_solution::set_element ( int index, double arg )
{
	if (index < 0 || index > original_size ) return 1;
	
	int erase_position = -1;
	int insert_position = -1;
	

	for (unsigned i = 0 ;i < indices.size();++i){
		if (indices[i]==index) {
			elements[i] = arg;
			if (double_equality(arg,0)) erase_position = i;
			break;
		}
		else if (indices[i]>index) {
			insert_position = i;
			break;
		}
	}
	
	if (insert_position>=0 && double_inequality(arg,0) ){
		std::vector<int>::iterator it = indices.begin()+insert_position;
		std::vector<double>::iterator dit = elements.begin()+insert_position;
		indices.insert(it,index);
		elements.insert(dit,arg);
	}
	if(erase_position>=0){
		std::vector<int>::iterator it = indices.begin()+erase_position;
		std::vector<double>::iterator dit = elements.begin()+erase_position;
		indices.erase(it);
		elements.erase(dit);
	}
	size = indices.size();

	return 0;
// 	else

}

int My_solution::element_pe ( int index, double arg )
{
	if (index < 0 || index > original_size ) return 1;
	if (double_equality(arg,0)) return 2;
	
	int insert_position = -1;
	int erase_position = -1;
	bool set_done =false;
	for (unsigned i = 0 ;i < indices.size();++i){
		if (indices[i]==index) {
// 			cout << "element[i]: " << elements[i]<<endl;
			elements[i] += arg;
// 			cout << "element[i]: " << elements[i]<<endl;
			set_done = true;
			if (double_equality(elements[i],0)) erase_position = i;
			break;
		}
		else if (indices[i]>index) {
			set_done = true;
			insert_position = i;
			break;
		}
	}
	if(!set_done){
		insert_position = indices.size();
	}
	
	
// 	cout << "insert_position: " << insert_position << " erase_position: " << erase_position << endl;
	if (insert_position>=0){
// 		cout << "insert_position: " << insert_position << " erase_position: " << erase_position << endl;
		
		std::vector<int>::iterator it = indices.begin()+insert_position;
		std::vector<double>::iterator dit = elements.begin()+insert_position;
		indices.insert(it,index);
		elements.insert(dit,arg);
	}
	if(erase_position>=0){
// 		cout << "insert_position: " << insert_position << " erase_position: " << erase_position << endl;
		
		std::vector<int>::iterator it = indices.begin()+erase_position;
		std::vector<double>::iterator dit = elements.begin()+erase_position;
		indices.erase(it);
		elements.erase(dit);
	}
	size = indices.size();
	
	
	return 0;
}

int My_solution::element_me ( int index, double arg )
{
	return element_pe(index,-arg);
}
int My_solution::element_mm ( int index )
{
	return element_pe(index,-1);
}

int My_solution::element_pp ( int index )
{
	return element_pe(index,1);
}





vector<int> My_solution::get_indices() const {
	return indices;
}
void My_solution::set_indices( vector<int>  arg ) {
	indices = arg;
	return;
}

int My_solution::get_size() const {
	return size;
}
void My_solution::set_size( int arg ) {
	size = arg;
	return;
}

int My_solution::get_original_size() const {
	return original_size ;
}
void My_solution::set_original_size( int arg ) {
	original_size = arg;
	return;
}

double My_solution::get_original_objective() const {
	return original_objective;
}
void My_solution::set_original_objective( double  arg ) {
	original_objective = arg;
	return;
}
double My_solution::get_auxiliary_objective() const {
	return auxiliary_objective;
}
void My_solution::set_auxiliary_objective( double arg ) {
	auxiliary_objective = arg;
	return;
}


My_solution::My_solution():size(0) {}
My_solution::My_solution( vector<double> arg ) {
	original_size = arg.size();
	for( int i = 0; i < original_size; ++i ) {
		if( !double_equality( arg[i],0 ) ) {
			indices.push_back( i );
			elements.push_back( arg[i] );
		}
	}
	// 	normalize();
	original_objective = 0;
	auxiliary_objective = 0;
	size = elements.size();

}

My_solution::My_solution( int ncols, const double *arg ) {
	original_size = ncols;
	for( int i = 0; i < original_size; ++i ) {
		if( !double_equality( arg[i],0 ) ) {
			indices.push_back( i );
			elements.push_back( arg[i] );
		}
	}
	// 	normalize();
	size = elements.size();
}


My_solution::My_solution( const My_solution &other ) {
	original_size = other.original_size;
	size= other.size;
	elements = other.elements;
	indices = other.indices;
	original_objective = other.original_objective;
	auxiliary_objective = other.auxiliary_objective;
}

My_solution My_solution::operator=( const My_solution &other ) {
	original_size = other.original_size;
	size= other.size;
	elements = other.elements;
	indices = other.indices;
	original_objective = other.original_objective;
	auxiliary_objective = other.auxiliary_objective;
	return *this;
}


void My_solution::print( std::ostream & out ) {
	out << "Org_obj: " << original_objective << " | ";
	// 	out << "size: " << size << endl;
	if (size <=0){
		out << endl;
	}
	for( int i =0; i < size; ++i ) {
		out << ", X_" << indices[i] << "= "<< elements[i];
	}
	out << endl;
	return;
}



void My_solution::shortlineprint( std::ostream & out ) {
// 	out << "Org_obj: " << original_objective << " | ";
		
	
	out << "original_size: " << original_size << endl;
	out << "size: " << size << endl;
		
	for( int i =0; i < size; ++i ) {
		out << "["<< indices[i] << "]="<< elements[i];
	}
	out << endl;
	return;
}


void My_solution::shortlineprint_selected_indices( vector<int> selected_indices,bool pos, std::ostream & out ) {
	// 	out << "Org_obj: " << original_objective << " | ";
	// 	out << "size: " << size << endl;
// 	unsigned s_ind = 0;
	unsigned m_ind = 0;
	
// 	static int it = 0;
// 	it++;
	
// 	vector<double> el = get_solution_vector();
// 	for (int i =0; i < selected_indices.size();++i){
// 		out << el[selected_indices[i]] << " ";
// 		
// 	}
// 	out<<endl;
// // 	vector_print<double>(el,out);
	
// 	out <<"call " << it <<"\t";
	if (!pos){
		for (unsigned i =0; i < selected_indices.size();++i){
			int s_ind = selected_indices[i];
			if (m_ind < indices.size()){
				if (s_ind < indices[m_ind]){
					out << ".";
	// 				continue;
				}
				else{ 
					if (s_ind == indices[m_ind]){
						out << elements[m_ind] /*<< " "*/;
						m_ind++;
					}
					else{ /** s_ind > m_ind*/
						m_ind++;
						i--;
					}
				}
			}
			else{
				out << ".";
			}
		}
	}
	else{
		for (unsigned i =0; i < selected_indices.size();++i){
			int s_ind = selected_indices[i];
			if (m_ind < indices.size()){
				if (s_ind < indices[m_ind]){
// 					out << 0 << " ";
					// 				continue;
				}
				else{ 
					if (s_ind == indices[m_ind]){
						out << "["<< indices[m_ind] << "]" << elements[m_ind] << " ";
						m_ind++;
					}
					else{ /** s_ind > m_ind*/
						m_ind++;
						i--;
					}
				}
			}
			else{
// 				out << "x" << " ";
			}
		}
	}
	out << endl;
	return;
}


void My_solution::normalize( std::ostream &out ) {
	double norm =0;
	if( norm <  1.0e-20 ) {
		norm =0;
		for( unsigned i =0; i < elements.size(); ++i ) {
			norm += elements[i]*elements[i];
		}
	}
	divide_coefficients( sqrt( norm ) );
}

void My_solution::divide_coefficients( double x ) {
	if( x<0 ) {
		x*=-1;
	}

	for( unsigned i = 0; i < elements.size(); ++i ) {
		elements[i]/=x;
	}

	return;
}
vector< double > My_solution::get_solution_vector() const {
	vector<double> retval;
// 	std::cout << "original_size "<< original_size <<std::endl;
// 	std::cout << "size "<< size <<std::endl;
	retval.resize( original_size,0 );

	for( int i =0; i < size; ++i ) {
		retval[indices[i]] = elements[i];
	}
	return retval;
}


#ifdef USING_MPI

void My_solution::_send_mpi( int destination, int tag, vector<MPI_Request> &v_request, vector<MPI_Status> v_status, MPI_SEND_TYPE SEND_TYPE ) {

	if( tag != TAG_SENDING_SOLUTION ) {

		if( tag != TAG_TESTING ) {
			cout << "ERROR: tag: " << tag << " BUT should be TAG_SENDING_SOLUTION: " <<  TAG_SENDING_SOLUTION << endl;
		}
	}
	v_request.resize( 6,MPI_REQUEST_NULL );
	v_status.resize( 6 );

	int* indices2 = ( int * ) malloc( size * sizeof( int ) );
	double* elements2 = ( double * ) malloc( size * sizeof( double ) );
	for( int i =0; i < size; ++i ) {
		indices2 [i] = indices[i];
		elements2[i] = elements[i];
	}

	COMMON_DATA_SEND_SOLUTION


	free( indices2 );
	free( elements2 );
	return;
}

void My_solution::_recv_mpi( int source, int tag, vector<MPI_Request> &v_request, vector<MPI_Status> v_status, MPI_RECV_TYPE RECV_TYPE ) {
	v_request.resize( 6,MPI_REQUEST_NULL );
	v_status.resize( 6 );
	if( tag != TAG_SENDING_SOLUTION ) {
		if( tag != TAG_TESTING ) {
			cout << "ERROR: tag: " << tag << " BUT should be TAG_SENDING_SOLUTION: " <<  TAG_SENDING_SOLUTION << endl;
		}
	}
	COMMON_DATA_RECV_SOLUTION
	return;
}
#else
void My_solution::_send_alt1(int destination, int tag, my_signaler comm, ALT1_SEND_TYPE SEND_TYPE ){
	TAG_CONTROL(TAG_SENDING_SOLUTION)
	COMMON_DATA_SERIALIZE_SOLUTION 
	//we have s as stringstream
	comm.send_object(destination, tag, s.str(),SEND_TYPE);
	return;

}
void My_solution::_recv_alt1(int source, int tag, my_signaler comm, ALT1_RECV_TYPE RECV_TYPE ){
	TAG_CONTROL(TAG_SENDING_SOLUTION)
	
	int rtag  = tag;
	//cout << " printing at " << __FILE__ << " and line "<< __LINE__ <<endl; 
	string sr = comm.recv_object(source, rtag, RECV_TYPE);	
	//cout << " printing at " << __FILE__ << " and line "<< __LINE__ <<endl; 
	
	COMMON_DATA_UNSERIALIZE_SOLUTION(sr);
	
	return;
}
#endif
bool My_solution::operator== ( const My_solution& rhs ) const{
	if (size != rhs.size) return false;
	
	for (int i = 0 ; i < size; ++i){
		if (indices[i] != rhs.indices[i]) return false;
		if (double_inequality(elements[i],rhs.elements[i])) return false;
	}
	
	return true;
	
}

bool My_solution::operator< ( const My_solution& rhs ) const{
	if (size < rhs.size) return true;
	if (size > rhs.size) return false;
	
	for (int i = 0 ; i < size; ++i){
		if (indices[i] < rhs.indices[i]) return true;
		if (indices[i] > rhs.indices[i]) return false;
	}

	for (int i = 0 ; i < size; ++i){
		if (double_equality(elements[i],rhs.elements[i])) continue;
		else return (elements[i] < rhs.elements[i]);
	}
	return false;
}



bool My_solution::equality_on_selected_indices ( const My_solution& rhs, vector< int > selected_indices ) const
{
	if (selected_indices.size() == 0 || (int)selected_indices.size() == original_size) {
// 		cout << " ZZZZZZZZZZZZZZZZZZZ " << endl;	
		return (*this==rhs);
	}
	int index_own = 0;
	int index_rhs = 0;
// 	int index_selected = 0;
	
	bool while_own= true;
	bool while_rhs= true;
	
	for (unsigned i = 0; i < selected_indices.size(); ++i){
		int ic = selected_indices[i];
		while_own = true;
		
		while(while_own){
			while_own = false;
			if (index_own< size){
				if (indices[index_own]<ic){
					index_own++;
					while_own = true;
				}
			}
		}
// 		while (indices[index_own]<ic && index_own<size) index_own++;
		/** indices[index_own] >= ic index to be compared */ 

		while_rhs= true;
		
		while(while_rhs){
			while_rhs = false;
			if (index_rhs<rhs.size){
				if (rhs.indices[index_rhs]<ic ){
					index_rhs++;
					while_rhs= true;
				}
			}
		}
// 		while (rhs.indices[index_rhs]<ic && index_rhs<rhs.size) index_rhs++;
		/** rhs.indices[index_rhs] >= ic index to be compared */ 
		
		
		
		if (index_own>=size && index_rhs>=rhs.size) 
			return true; /** both solutions have 0s after index_own and index_rhs*/
		
		
		if (index_own>=size && index_rhs<rhs.size) {/** own done check for rhs with selected_indices*/
			if (rhs.indices[index_rhs] == ic) 
				return false; /** own done rhs has a nonzero in selected_indices*/
			else 
				continue; /** else continue with the next ic */
		}
		
		if (index_own<size && index_rhs>=rhs.size) {/** rhs done check for own with selected_indices*/
			if (indices[index_own] == ic) 
				return false; /** rhs done own has a nonzero in selected_indices*/ 
			else
				continue;/** else continue with the next ic */
		}
		
		
		/** NONE DONE YET **/

		/** one is nonzero other is zero */ 
		if (indices[index_own] == ic && rhs.indices[index_rhs] > ic) 
			return false; /** own has a nonzero value on index ic rhs is zero*/
		if (indices[index_own] > ic && rhs.indices[index_rhs] == ic) 
			return false; /** rhs has a nonzero value on index ic own is zero*/
		
		if (indices[index_own] > ic && rhs.indices[index_rhs] > ic) 
			continue; /** both are 0*/ 
		
			
		/** both indices are equal to ic => both are nonzero, compare elements*/
		if( double_inequality(elements[index_own],rhs.elements[index_rhs])) 
			return false; /** unequal elements return false */
		
		
	}
	
	return true;
	
}


bool My_solution::is_empty()
{
	
	if (original_size <= 0) return true;
	
	if ( original_size > 0 ) {
		return false;
	}
	
	return true;
}

void My_solution::clear()
{
	indices.clear();
	elements.clear();
	size = 0;
	original_objective = -1e99;
	auxiliary_objective = -1e99;
	
}




// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4;
