#include "myserializablerowcut.h"
using namespace std;




void mySerializableRowCut::internal_sort() {
	vector<pair<int,double> > v;
	for( int i = 0; i < size; ++i ) {
		pair<int,double> asd( indices[i],elements[i] );
		v.push_back( asd );
	}
	sort( v.begin(),v.end() );
	for( int i = 0; i < size; ++i ) {
		indices[i] = v[i].first;
		elements[i] = v[i].second;
	}
	generate_OsiRowCut();
	return;
}

mySerializableRowCut::mySerializableRowCut() {};
mySerializableRowCut::mySerializableRowCut( const mySerializableRowCut &other ) {
	lb = other.lb;
	ub = other.ub;
	sense = other.sense;
	size = other.size;
// 	original_size = other.original_size;
	elements = other.elements;
	indices = other.indices;
	norm = other.norm;
	for( int i = CUT_EFFICIENCY_BGN; i < CUT_EFFICIENCY_END; ++i ){
		efficiency[i] = other.efficiency[i];
	}
}


mySerializableRowCut& mySerializableRowCut::operator= ( const mySerializableRowCut& other )
{
	lb = other.lb;
	ub = other.ub;
	sense = other.sense;
	size = other.size;
	// 	original_size = other.original_size;
	elements = other.elements;
	indices = other.indices;
	norm = other.norm;
	for( int i = CUT_EFFICIENCY_BGN; i < CUT_EFFICIENCY_END; ++i ){
		efficiency[i] = other.efficiency[i];
	}
	return *this;
}



mySerializableRowCut::mySerializableRowCut( vector<int> ind, vector<double> ele, double llb, double uub ) {
	lb = llb;
	ub = uub;
	size = ind.size();
	elements = ele ;

	indices = ind;
	sense = 'U';
	generate_OsiRowCut();
	reset_efficiency();
}

mySerializableRowCut::mySerializableRowCut( OsiRowCut* cut2 ) {

	lb = cut2->lb();
	ub = cut2->ub();
	sense = cut2->sense();

	packedvector = CoinPackedVector( cut2->mutableRow() );
	size = packedvector.getNumElements();

	int *indices2;
	double *elements2;

	indices2  = packedvector.getIndices();
	elements2 = packedvector.getElements();

	elements.resize( size );
	indices.resize( size );

	for( int i = 0; i<size; ++i ) {
		elements[i] = elements2[i];
		indices[i] = indices2[i];
	}
	norm = 0;
	internal_sort();
	normalize();
	reset_efficiency();
}

mySerializableRowCut::mySerializableRowCut( OsiRowCut cut ) {

	rowcut = OsiRowCut( cut );
	lb = rowcut.lb();
	ub = rowcut.ub();
	sense = rowcut.sense();


	packedvector = rowcut.mutableRow();

	size = packedvector.getNumElements();

	int *indices2;
	double *elements2;

	indices2  = packedvector.getIndices();
	elements2 = packedvector.getElements();

	elements.resize( size );
	indices.resize( size );

	for( int i = 0; i<size; ++i ) {
		elements[i] = elements2[i];
		indices[i] = indices2[i];
	}
	norm = 0;
	internal_sort();
	normalize();
	reset_efficiency();
}



void mySerializableRowCut::print( std::ostream &out,bool eff ) {
	if( size < 1 ) {
		out << "size < 1 nothing to print " << endl;

		out << "elements.size() " << elements.size() << endl;
		out << "indices.size() " << indices.size() << endl;
		out << "lb "<< lb << endl;
		out << "ub "<< ub << endl;
		return;
	}

	if( lb > -DBL_MAX )
		out << lb << " <= ";

	for( int i = 0; i<size; ++i ) {
		if( elements[i] > 0 )
			out << "+" << elements[i] << " * x_" << indices[i] << " " ;
		else
			out << elements[i] << " * x_" << indices[i] << " " ;
	}

	if (eff){
		if( ub < DBL_MAX ) {
			out << " <= " << ub << " Eff: ";
			for( int i = 0 ; i < CUT_EFFICIENCY_END; ++i ) {
				out << efficiency[i] << " ";
			}
			out << "\t NORM : "<< norm << endl;
		}
		else {
			out << " Eff: ";
			for( int i = 0 ; i < CUT_EFFICIENCY_END; ++i ) {
				out << efficiency[i] << " ";
			}
			out << "\t NORM : "<< norm << endl;
		}
	}
	else{
		if( ub < DBL_MAX ) {
			out << " <= " << ub  << "\t NORM : "<< norm << endl;
		}
		else {
			out  << "\t NORM : "<< norm << endl;
		}
	}
	return;
}

CoinPackedVector mySerializableRowCut::generate_CoinPackedVector( std::ostream &out ) {
	int *indices2 = new int[size];
	double *elements2 = new double[size];
	for( int i = 0; i<size; ++i ) {
		elements2[i] = elements[i];
		indices2[i] = indices[i];
	}
// 			out << "before calling CoinPackedVector(size,indices2,elements2);" << endl;
	packedvector = CoinPackedVector( size,indices2,elements2 );
// 			out << "after calling CoinPackedVector(size,indices2,elements2);" << endl;

	delete[] indices2;
	delete[] elements2;

// 			out << "deleted indices2 and elements2  returning " << endl;
	return packedvector;

}

OsiRowCut* mySerializableRowCut::generate_OsiRowCutPointer() {

	rowcut.setRow( generate_CoinPackedVector() );

	/* make this more complicated wrt sense */
	rowcut.setLb( lb );
	rowcut.setUb( ub );
	return &rowcut;
}

OsiRowCut mySerializableRowCut::generate_OsiRowCut( std::ostream &out ) {
// 	out << "setRow "<< endl;
	rowcut.setRow( generate_CoinPackedVector( out ) );
// 	out << "setLB "<< endl;

	/* make this more complicated wrt sense */
	rowcut.setLb( lb );
// 	out << "setUB "<< endl;

	rowcut.setUb( ub );

// 	out << "return "<< endl;
// 	out << "rowcut.print() "<< endl;

// 	rowcut.print();
// 	out << "rowcut.print() end"<< endl;

// 	out << "this->print(out)" << endl;
// 	this->print();
	return rowcut;
}

void mySerializableRowCut::normalize( std::ostream &out ) {
// 	return;
	if( norm <  1.0e-20 ) {
		norm =0;
		for( unsigned i =0; i < elements.size(); ++i ) {
			norm += elements[i]*elements[i];
		}
	}
// 	double z = norm;
// 	divide_coefficients( sqrt( norm ) );
// 	norm = z;
// 	divide_coefficients( sqrt( norm ) );
    return;
}


void mySerializableRowCut::divide_coefficients( double x ) {
	if( x<0 ) {
		x*=-1;
	}

	for( unsigned i = 0; i < elements.size(); ++i ) {
		elements[i]/=x;
	}

	if( lb > -1e30 ) {
		lb/=x;
	}

	if( ub < 1e30 ) {
		ub/=x;
	}
	
	norm /= ( x*x );

	return;
}

void mySerializableRowCut::multiply_coefficients ( double x ){
	for( unsigned i = 0; i < elements.size(); ++i ) {
		elements[i]*=x;
	}
	if( lb > -1e30) {
		lb*=x;
	}
	
	if( ub < 1e30 ) {
		ub*=x;
	}
	norm *= ( x*x );
	return;
}
bool mySerializableRowCut::operator<( const mySerializableRowCut &rhs ) const {


	if( size != rhs.size )
		return ( size < rhs.size );

	for( int i = 0; i < size; ++i ) {
		if( indices[i] != rhs.indices[i] )
			return ( indices[i] < rhs.indices[i] );
		if( fabs( elements[i] - rhs.elements[i] ) > __CUT_EQUALITY_TOLERANCE)
			return ( elements[i] < rhs.elements[i] );
	}
	
	
// 	if( fabs( lb - rhs.lb ) > __CUT_EQUALITY_TOLERANCE )
// 		return ( lb < rhs.lb );
// 	if( fabs( ub - rhs.ub ) > __CUT_EQUALITY_TOLERANCE )
// 		return ( ub < rhs.ub );
// 	
// 	if( fabs( norm - rhs.norm ) < __CUT_EQUALITY_TOLERANCE )
// 		return ( norm<rhs.norm );

	return false;

}

// NEED TO CHECK IF LBs or UBs = inf and sense is not an issue;
bool mySerializableRowCut::operator==( const mySerializableRowCut &rhs ) const {

	
	if( size != rhs.size )
		return ( false );
	
	double ratio = elements[0] / rhs.elements[0];
	for( int i = 0; i < size; ++i ) {
		if( indices[i] != rhs.indices[i] )
			return ( false );
		if( double_inequality(ratio, elements[i]/rhs.elements[i], __CUT_EQUALITY_TOLERANCE ) )
			return ( false );
	}

// 			if (sense != rhs.sense) return (false);
// 	if( fabs( lb - rhs.lb ) > __CUT_EQUALITY_TOLERANCE )
// 		return ( false );
// 	if( fabs( ub - rhs.ub ) > __CUT_EQUALITY_TOLERANCE )
// 		return ( false );
// 	if( fabs( norm - rhs.norm ) > __CUT_EQUALITY_TOLERANCE )
// 		return ( false );

	return true;
}

mySerializableRowCut::mySerializableRowCut( double arg_lb,
        double arg_ub,
        char arg_sense,
        int arg_size,
        vector<double> arg_elements,
        vector<int> arg_indices	):
	lb( arg_lb ),
	ub( arg_ub ),
	sense( arg_sense ),
	size( arg_size ),
	indices( arg_indices ),
	elements( arg_elements ) {}

double mySerializableRowCut::get_lb() const {
	return lb;
}
void mySerializableRowCut::set_lb( double arg ) {
	lb = arg;
	return;
}

double mySerializableRowCut::get_ub() const {
	return ub;
}
void mySerializableRowCut::set_ub( double arg ) {
	ub = arg;
	return;
}

double mySerializableRowCut::get_norm() const {
	return norm;
}
void mySerializableRowCut::set_norm( double arg ) {
	norm = arg;
	return;
}

char mySerializableRowCut::get_sense() const {
	return sense;
}
void mySerializableRowCut::set_sense( char arg ) {
	sense = arg;
	return;
}

int mySerializableRowCut::get_size() const {
	return size;
}
void mySerializableRowCut::set_size( int arg ) {
	size = arg;
	return;
}

// int mySerializableRowCut::get_original_size() const { return original_size;}
// void mySerializableRowCut::set_original_size(int arg) { original_size= arg; return; }

vector<double>  mySerializableRowCut::get_elements() const {
	return elements;
}
void mySerializableRowCut::set_elements( vector<double>  arg ) {
	elements = arg;
	return;
}

vector<int> mySerializableRowCut::get_indices() const {
	return indices;
}
void mySerializableRowCut::set_indices( vector<int>  arg ) {
	indices = arg;
	return;
}

CoinPackedVector mySerializableRowCut::get_packedvector() const {
	return packedvector;
}
void mySerializableRowCut::set_packedvector( CoinPackedVector arg ) {
	packedvector = arg;
	return;
}

OsiRowCut mySerializableRowCut::get_rowcut() const {
	return rowcut;
}
void mySerializableRowCut::set_rowcut( OsiRowCut arg ) {
	rowcut = arg;
	return;
}

void mySerializableRowCut::reset_efficiency() {
	for( int i =CUT_EFFICIENCY_BGN; i < CUT_EFFICIENCY_END; ++i ) {
		efficiency[i] = 0;
	}
	return;
}

void mySerializableRowCut::set_efficiency( double arg[CUT_EFFICIENCY_END] ) {
	for( int i =CUT_EFFICIENCY_BGN; i < CUT_EFFICIENCY_END; ++i ) {
		efficiency[i] = arg[i];
	}
	return;
}

void mySerializableRowCut::set_efficiency(CUT_EFFICIENCY_ENM enm , double arg){
	efficiency[enm] = arg;
	return;
}

double mySerializableRowCut::get_efficiency( CUT_EFFICIENCY_ENM enm ) const {
	return efficiency[enm];
}
double* mySerializableRowCut::get_efficiency() {
	return efficiency;
}
#ifdef USING_MPI

void mySerializableRowCut::_send_mpi( int destination, int tag, vector<MPI_Request> &v_request, vector<MPI_Status> v_status, MPI_SEND_TYPE SEND_TYPE ) {
	if( tag != TAG_SENDING_ROWCUT ) {

		if( tag != TAG_TESTING ) {
			cout << "ERROR rank: " <<MPI::COMM_WORLD.Get_rank() << " : tag: " << tag << " BUT should be TAG_SENDING_ROWCUT: " <<  TAG_SENDING_ROWCUT << " FILE: " << __FILE__<< " LINE: " << __LINE__ << endl ;
		}
	}
	v_request.resize( 8,MPI_REQUEST_NULL );
	v_status.resize( 8 );

	int* indices2 = ( int * ) malloc( size * sizeof( int ) );
	double* elements2 = ( double * ) malloc( size * sizeof( double ) );
	for( int i =0; i < size; ++i ) {
		indices2 [i] = indices[i];
		elements2[i] = elements[i];
	}

	COMMON_DATA_SEND_ROWCUT


	free( indices2 );
	free( elements2 );
	return;
}


void mySerializableRowCut::_recv_mpi( int source, int tag, vector<MPI_Request> &v_request, vector<MPI_Status> v_status, MPI_RECV_TYPE RECV_TYPE ) {
	v_request.resize( 8,MPI_REQUEST_NULL );
	v_status.resize( 8 );
	if( tag != TAG_SENDING_ROWCUT ) {
		if( tag != TAG_TESTING ) {
			cout << "ERROR: tag: " << tag << " BUT should be TAG_SENDING_ROWCUT: " <<  TAG_SENDING_ROWCUT << " FILE: " << __FILE__<< " LINE: " << __LINE__ << endl ;
		}
	}
	COMMON_DATA_RECV_ROWCUT

	return;
}



#else
void mySerializableRowCut::_send_alt1(int destination, int tag, my_signaler comm, ALT1_SEND_TYPE SEND_TYPE ){
	TAG_CONTROL(TAG_SENDING_ROWCUT)
	COMMON_DATA_SERIALIZE_ROWCUT
	//we have s as stringstream
	comm.send_object(destination, tag, s.str(),SEND_TYPE);
	return;

}
void mySerializableRowCut::_recv_alt1(int source, int tag, my_signaler comm, ALT1_RECV_TYPE RECV_TYPE ){
    TAG_CONTROL(TAG_SENDING_ROWCUT)

	
	int rtag  = tag;
	string sr = comm.recv_object(source, rtag, RECV_TYPE);	
	COMMON_DATA_UNSERIALIZE_ROWCUT(sr)
	return;
}

#endif
void mySerializableRowCut::calculate_efficiency( const double* x ) {

	reset_efficiency();
	efficiency[CUT_EFFICIENCY_VIOLATION] = -ub;

	double bar_a = 0;
	double multiplication_of_aj= 1;
	for( int i = 0; i < size ; ++i ) {

		efficiency[CUT_EFFICIENCY_NORM_ALPHA] += ( elements[i]*elements[i] );
		efficiency[CUT_EFFICIENCY_VIOLATION] += ( elements[i]*x[indices[i]] );
		if( double_inequality(x[indices[i]],0)){
			bar_a += ( elements[i]*elements[i]);
		}
		multiplication_of_aj*= elements[i];
	}
	efficiency[CUT_EFFICIENCY_NORM_ALPHA] = sqrt( efficiency[CUT_EFFICIENCY_NORM_ALPHA] );
	efficiency[CUT_EFFICIENCY_RELATIVE_VIOLATION] = efficiency[CUT_EFFICIENCY_VIOLATION]/abs( ub );
	efficiency[CUT_EFFICIENCY_DISTANCE] = efficiency[CUT_EFFICIENCY_VIOLATION] / efficiency[CUT_EFFICIENCY_NORM_ALPHA];
	efficiency[CUT_EFFICIENCY_ADJUSTED_DISTANCE] = efficiency[CUT_EFFICIENCY_VIOLATION] /
													(1+sqrt(bar_a));

	efficiency[CUT_EFFICIENCY_DISTANCE_VARIANT] = pow(
										efficiency[CUT_EFFICIENCY_VIOLATION], CUT_EFFICIENCY_K) / pow(abs(multiplication_of_aj),(CUT_EFFICIENCY_K/size));


	return;
}

// void mySerializableRowCut::calculate_efficiency_v(const vector<double> x ) {
// 	reset_efficiency();
// 	efficiency[CUT_EFFICIENCY_VIOLATION] = -ub;
// 	for( int i = 0; i < size ; ++i ) {
// 		efficiency[CUT_EFFICIENCY_NORM_ALPHA] += ( elements[i]*elements[i] );
// 		efficiency[CUT_EFFICIENCY_VIOLATION] += ( elements[i]*x[indices[i]] );
// 	}
// 	efficiency[CUT_EFFICIENCY_NORM_ALPHA] = sqrt( efficiency[CUT_EFFICIENCY_NORM_ALPHA] );
// 	efficiency[CUT_EFFICIENCY_RELATIVE_VIOLATION] = efficiency[CUT_EFFICIENCY_VIOLATION]/abs( ub );
// 	efficiency[CUT_EFFICIENCY_DISTANCE] = efficiency[CUT_EFFICIENCY_VIOLATION] / efficiency[CUT_EFFICIENCY_NORM_ALPHA];
// 	return;
// }

void mySerializableRowCut::calculate_efficiency_wrt_objective(My_objective_function &obj, bool is_org) {
	unsigned  add = 0;
	if (!is_org) add++;

	for (int i = 0; i < obj.get_size();++i){
		efficiency[CUT_EFFICIENCY_NORM_C_ORG + add] +=  (obj.get_element(i)*obj.get_element(i));
	}
	efficiency[CUT_EFFICIENCY_NORM_C_ORG + add] = sqrt(efficiency[CUT_EFFICIENCY_NORM_C_ORG + add] );
	int i = 0;
	int j = 0;

	while( (i<size )&& (j < obj.get_size())){
		if(indices[i] == obj.get_index(j)){
			efficiency[CUT_EFFICIENCY_OBJECTIVE_PARALLELISM_WRT_ORG+add] +=
			(elements[i]*obj.get_element(j));
			++i;++j;
		}
		else{
			if (indices[i] > obj.get_index(j)){
				++j;
			}
			else {
				++i;
			}
		}
	}
	efficiency[CUT_EFFICIENCY_OBJECTIVE_PARALLELISM_WRT_ORG+add] /= (efficiency[CUT_EFFICIENCY_NORM_C_ORG+add]*efficiency[CUT_EFFICIENCY_NORM_ALPHA]);

	efficiency[CUT_EFFICIENCY_EXPECTED_IMPROVEMENT_WRT_ORG+add] =
		efficiency[CUT_EFFICIENCY_NORM_C_ORG+add]*
		efficiency[CUT_EFFICIENCY_OBJECTIVE_PARALLELISM_WRT_ORG+add] *
		efficiency[CUT_EFFICIENCY_DISTANCE];


	efficiency[CUT_EFFICIENCY_SUPPORT] = 1- (size/ ((double) obj.get_original_size()));

	return;
}

void mySerializableRowCut::update_efficiency_integral_support(const vector<int> &column_types, int n_binary_plus_n_int){
	int integral_support = 0;
	for (int i = 0; i < size; ++i){
		if(column_types[indices[i]] > 0) integral_support++;
	}
	efficiency[CUT_EFFICIENCY_INTEGRAL_SUPPORT] = (((double)integral_support)/((double)n_binary_plus_n_int) );
	return;
}


