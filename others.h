#ifndef OTHERS_H
#define OTHERS_H


#include "parallel_cpfp_stuff.h"

#include <random>



class Seeder{
	unsigned seed;
	minstd_rand0 my_seed_seeder;
public:
	Seeder(unsigned arg_seed = my_default_seed);
    Seeder &operator=(const Seeder &other);
    
	unsigned nextSeed();
	unsigned nextSeed_lcm();
	unsigned lth_seed_lcm(int l);

};

class Objective_type_distributor{
        int total_number_of_slaves;
        int number_of_instances_of_type[OBJECTIVE_TYPE_END];
		double perturbation;
	public:
		Objective_type_distributor();
        Objective_type_distributor(const Objective_type_distributor &other);

		Objective_type_distributor &operator=(const Objective_type_distributor &other);

        int add_instance(OBJECTIVE_TYPE_ENUM ENM, int number);

        int I_am_running_what_objective(int slave_id) const;
        int get_number_of_instances_of_type(OBJECTIVE_TYPE_ENUM ENM) const;
        int get_number_of_instances_that_need_AC() const;
        int get_number_of_instances_that_does_not_need_AC() const;
		vector<int> get_type_per_slave() const;
		double get_perturbation() const;
		void set_perturbation(double arg);
		void print(ostream &out =std::cout);
		void shortlineprint(ostream &out =std::cout, int line =0);
        string string_print(int line=0);
		int  I_am_running_what_objective_order(int slave_id) const;
		
};


#endif
// kate: indent-mode cstyle; indent-width 1; replace-tabs on;
