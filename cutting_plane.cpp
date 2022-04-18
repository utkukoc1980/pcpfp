
#include <iostream>
#include <string>


#include "parallel_cpfp_stuff.h"
#include "myserializablerowcut.h"
// #include "/home/utkukoc/coin-Cgl/build/include/coin/OsiCuts.hpp"
// #include "OsiCuts.hpp"
//
// #include "OsiSolverInterface.hpp"
// #include "OsiSolverParameters.hpp"
// #include "OsiClpSolverInterface.hpp"
// #include "CglRedSplit.hpp"
// //
// #include "CoinPackedVectorBase.hpp"
// #include "CglKnapsackCover.hpp"
// #include "CglSimpleRounding.hpp"
//
//
// #include "CglGomory.hpp"

using namespace std;



int mainCP (int argc, char* argv[]){

	char *f_name_lp = argv[1];
	OsiClpSolverInterface *clp = new OsiClpSolverInterface;
	int ncols;

	CglGomory gomory_cut_generator;
	CglKnapsackCover knapsack_cut_generator;
	CglRedSplit redsplit_cut_generator;
	CglSimpleRounding simplerounding_cut_generator;



	if(strcmp(&(f_name_lp[strlen(f_name_lp)-3]), ".lp") == 0) {
		clp->readLp(f_name_lp);
	}
	else {
		if(strcmp(&(f_name_lp[strlen(f_name_lp)-4]), ".mps") == 0) {
			clp->readMps(f_name_lp);
		}
		else {
			printf("### ERROR: unrecognized file type\n");
			exit(1);
		}
	}


	ncols = clp->getNumCols();
	clp->initialSolve();
	OsiCuts cuts;
	cout << "generating cuts" << endl;

	gomory_cut_generator.generateCuts(*clp, cuts);
	knapsack_cut_generator.generateCuts(*clp, cuts);

	redsplit_cut_generator.generateCuts(*clp, cuts);
	simplerounding_cut_generator.generateCuts(*clp, cuts);

	vector<mySerializableRowCut> my_cut_list;
	int number_of_cuts = cuts.sizeCuts() ;

	cout << cuts.sizeCuts() << endl;

	for (int i = 0; i < number_of_cuts; ++i){
		mySerializableRowCut cut3(cuts.rowCut(i));
		my_cut_list.push_back(cut3);
	}


	cout << cuts.sizeCuts() << endl;

	for (int i = 0; i < number_of_cuts; ++i){
		my_cut_list[i].print();
	}

	CoinPackedVector v;
	OsiRowCut cut;



	cout << number_of_cuts<< " cuts generated" << endl;

	return 0;

}


// kate: indent-mode cstyle; indent-width 1; replace-tabs on; 
