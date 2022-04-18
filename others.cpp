#include "others.h"

Seeder::Seeder( unsigned arg_seed ):seed( arg_seed ) {
	my_seed_seeder.seed(arg_seed);
};

Seeder &Seeder::operator=( const Seeder &other){
	seed = other.seed;
	my_seed_seeder.seed(seed);
return *this;	
};

unsigned Seeder::nextSeed() {
	return ++seed;
};
unsigned Seeder::nextSeed_lcm() {
	unsigned retval = seed;
	seed = my_seed_seeder();
	return retval;
};
unsigned Seeder::lth_seed_lcm( int l ) {
	for( int i =0; i<l-1; ++i ) {
		nextSeed_lcm();
	}
	return nextSeed_lcm();
}

Objective_type_distributor::Objective_type_distributor() {
	for( int i =OBJECTIVE_TYPE_BGN; i<OBJECTIVE_TYPE_END; ++i ) {
		number_of_instances_of_type[i]=0;
	}
	total_number_of_slaves = 0;
	perturbation =0;
}
Objective_type_distributor::Objective_type_distributor( const Objective_type_distributor &other ) {
	for( int i =OBJECTIVE_TYPE_BGN; i<OBJECTIVE_TYPE_END; ++i ) {
		number_of_instances_of_type[i]=other.number_of_instances_of_type[i];
	}
	total_number_of_slaves = other.total_number_of_slaves;
	perturbation = other.perturbation;
}

Objective_type_distributor &Objective_type_distributor::operator=( const Objective_type_distributor &other ) {
	for( int i =OBJECTIVE_TYPE_BGN; i<OBJECTIVE_TYPE_END; ++i ) {
		number_of_instances_of_type[i]=other.number_of_instances_of_type[i];
	}
	total_number_of_slaves = other.total_number_of_slaves;
	perturbation = other.perturbation;
	return *this;
}

int Objective_type_distributor::add_instance( OBJECTIVE_TYPE_ENUM ENM, int number ) {
	number_of_instances_of_type[ENM] += number;
	total_number_of_slaves +=number;
	return total_number_of_slaves;
}

int  Objective_type_distributor::I_am_running_what_objective( int slave_id ) const{
	int subtotal = 0;

	if (slave_id < 0) return -2;
	for( int i = OBJECTIVE_TYPE_BGN; i <OBJECTIVE_TYPE_END; ++i ) {
		subtotal += number_of_instances_of_type[i];
		if( subtotal > (slave_id) ) {
			return i	;
		}
	}
	return -1;
}

int Objective_type_distributor::I_am_running_what_objective_order ( int slave_id ) const{
	int subtotal = 0;
	int pre_subtotal = -1;
	if (slave_id < 0) return -2;
	for( int i = OBJECTIVE_TYPE_BGN; i <OBJECTIVE_TYPE_END; ++i ) {
 		
		if (i>OBJECTIVE_TYPE_BGN) pre_subtotal += number_of_instances_of_type[i-1];
		
		subtotal += number_of_instances_of_type[i];
		
		
		if( subtotal > (slave_id) ) {
			return slave_id - pre_subtotal;
			return i	;
		}
	}
	return -1;
}

int Objective_type_distributor::get_number_of_instances_of_type( OBJECTIVE_TYPE_ENUM ENM ) const{
	return number_of_instances_of_type[ENM];
}
int Objective_type_distributor::get_number_of_instances_that_need_AC() const{
	int retval =0;
	for( int i =OBJECTIVE_TYPE_BGN+1; i<OBJECTIVE_TYPE_END; ++i ) {
		retval+= number_of_instances_of_type[i];

	}
	return retval;
}
int Objective_type_distributor::get_number_of_instances_that_does_not_need_AC() const{
	return number_of_instances_of_type[OBJECTIVE_TYPE_ORIGINAL];
}


vector<int> Objective_type_distributor::get_type_per_slave() const{
	vector<int> retval;
	for( int i = OBJECTIVE_TYPE_BGN; i<OBJECTIVE_TYPE_END; ++i ) {
		for( int j =0; j< number_of_instances_of_type[i]; ++j ) {
			retval.push_back( i );
		}
	}
	return retval;
}
double Objective_type_distributor::get_perturbation() const {
	return perturbation;
}
void Objective_type_distributor::set_perturbation( double arg ) {
	perturbation = arg;
	return;
}
void Objective_type_distributor::print( ostream &out ) {
	//out << "total_number_of_slaves: " << total_number_of_slaves;
	for( int i = OBJECTIVE_TYPE_BGN; i<OBJECTIVE_TYPE_END; ++i ) {
		out << number_of_instances_of_type[i] <<" ";
	}
	out << " sum: " << total_number_of_slaves << " pert: " << perturbation << endl;

}

void Objective_type_distributor::shortlineprint( ostream &out, int line  ) {
	//out << "total_number_of_slaves: " << total_number_of_slaves;
	stringstream s;
	s << " ";
	if(number_of_instances_of_type[OBJECTIVE_TYPE_ORIGINAL]>0)
		s << number_of_instances_of_type[OBJECTIVE_TYPE_ORIGINAL] << "O";
	
	int n_p= 0;
	n_p+= number_of_instances_of_type[OBJECTIVE_TYPE_PERTURBED_FROM_HIT_RUN];
	n_p+= number_of_instances_of_type[OBJECTIVE_TYPE_PERTURBED_FROM_LONG_DIKIN];
	n_p+= number_of_instances_of_type[OBJECTIVE_TYPE_PERTURBED_FROM_SHORT_DIKIN];

	int n_r= 0;
	n_r+= number_of_instances_of_type[OBJECTIVE_TYPE_RANDOM_FROM_HIT_RUN];
	n_r+= number_of_instances_of_type[OBJECTIVE_TYPE_RANDOM_FROM_LONG_DIKIN];
	n_r+= number_of_instances_of_type[OBJECTIVE_TYPE_RANDOM_FROM_SHORT_DIKIN];
	
	
	int n_w= 0;
	n_w+= number_of_instances_of_type[OBJECTIVE_TYPE_RANDOM_WALK_HIT_RUN];
	n_w+= number_of_instances_of_type[OBJECTIVE_TYPE_RANDOM_WALK_LONG_DIKIN];
	n_w+= number_of_instances_of_type[OBJECTIVE_TYPE_RANDOM_WALK_SHORT_DIKIN];
	
// 	if(number_of_instances_of_type[OBJECTIVE_TYPE_PERTURBED_FROM_HIT_RUN]>0)
// 		s << number_of_instances_of_type[OBJECTIVE_TYPE_ORIGINAL] << "O";
	if (n_p > 0)
		s << n_p << "P";
	if (n_r > 0)
		s << n_r << "R";
	if (n_w > 0)
		s << n_w << "W";
	
	if (line >0) s<< " ln: " << line;
	s << "  ";
	out << s.str().c_str();
	
// 	for( int i = OBJECTIVE_TYPE_BGN; i<OBJECTIVE_TYPE_END; ++i ) {
// 		out << number_of_instances_of_type[i] <<" ";
// 	}
// 	out << " sum: " << total_number_of_slaves << " pert: " << perturbation << endl;
	
}

string Objective_type_distributor::string_print(int line) {
	//out << "total_number_of_slaves: " << total_number_of_slaves;
	stringstream s;
	if(number_of_instances_of_type[OBJECTIVE_TYPE_ORIGINAL]>0)
		s << number_of_instances_of_type[OBJECTIVE_TYPE_ORIGINAL] << "O";
	
	int n_p= 0;
	n_p+= number_of_instances_of_type[OBJECTIVE_TYPE_PERTURBED_FROM_HIT_RUN];
	n_p+= number_of_instances_of_type[OBJECTIVE_TYPE_PERTURBED_FROM_LONG_DIKIN];
	n_p+= number_of_instances_of_type[OBJECTIVE_TYPE_PERTURBED_FROM_SHORT_DIKIN];
	
	int n_r= 0;
	n_r+= number_of_instances_of_type[OBJECTIVE_TYPE_RANDOM_FROM_HIT_RUN];
	n_r+= number_of_instances_of_type[OBJECTIVE_TYPE_RANDOM_FROM_LONG_DIKIN];
	n_r+= number_of_instances_of_type[OBJECTIVE_TYPE_RANDOM_FROM_SHORT_DIKIN];
	int n_w= 0;
	n_w+= number_of_instances_of_type[OBJECTIVE_TYPE_RANDOM_WALK_HIT_RUN];
	n_w+= number_of_instances_of_type[OBJECTIVE_TYPE_RANDOM_WALK_LONG_DIKIN];
	n_w+= number_of_instances_of_type[OBJECTIVE_TYPE_RANDOM_WALK_SHORT_DIKIN];
	
	
// 	if(number_of_instances_of_type[OBJECTIVE_TYPE_PERTURBED_FROM_HIT_RUN]>0)
// 		s << number_of_instances_of_type[OBJECTIVE_TYPE_ORIGINAL] << "O";
	if (n_p > 0)
		s << n_p << "P";
	if (n_r > 0)
		s << n_r << "R";	
	if (n_w > 0)
		s << n_w << "W";
	if (line >0) s<< " ln: " << line;
	
	return s.str();
	
	// 	for( int i = OBJECTIVE_TYPE_BGN; i<OBJECTIVE_TYPE_END; ++i ) {
	// 		out << number_of_instances_of_type[i] <<" ";
	// 	}
	// 	out << " sum: " << total_number_of_slaves << " pert: " << perturbation << endl;
	
}
// kate: indent-mode cstyle; indent-width 4; replace-tabs off; tab-width 4;
