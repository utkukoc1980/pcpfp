// Feasibility Pump
// Proof-of-concept implementation with ILOG Cplex 9/Concert Technology
//
// Reference papers:
//
// A Feasibility Pump heuristic for general Mixed-Integer Problems
// L. Bertacco, M. Fischetti, A. Lodi.
// Technical Report OR/05/5 DEIS Universita' di Bologna, 2005
//
// The Feasibility Pump.
// M. Fischetti, F. Glover, A. Lodi.
// Mathematical Programming 2005
// Digital Object Identi¯er (DOI) 10.1007/s10107-004-0570-3.
//
// Livio Bertacco (livio@bertacco.it), July 8th, 2005
// PhD Candidate, University of Padua
// All rights reserved.




#include <cmath>
#include "myfp.h"



using namespace std;

static char presolvevarname[255];

const char* my_FP::getVarName(IloNumVar var)
{
   const char* name = var.getName();
   if( name == 0 )
   {
      sprintf(presolvevarname, "<presolved_%d>", var.getId());
      return presolvevarname;
   }
   else
      return name;
}

void my_FP::findIntVars()
{
   for( int i = 0; i < nVars; i++ )
      switch(vars[i].getType())
      {
      case IloNumVar::Bool:
      case IloNumVar::Int:
         double l = ceil(vars[i].getLB()-epInt);
         double u = floor(vars[i].getUB()+epInt);
         if( l != vars[i].getLB() )
         {

					#ifdef DEBUG_XX
								  cout << "Tightening LB of " << getVarName(vars[i]) << " from " << vars[i].getLB() << " to " << l << endl;
						#endif
            vars[i].setLB(l);
         }
         if( u != vars[i].getUB() )
         {

					#ifdef DEBUG_XX
								             cout << "Tightening UB of " << getVarName(vars[i]) << " from " << vars[i].getUB() << " to " << u << endl;
						#endif
            vars[i].setUB(u);
         }
         if( l < u )
         {
            intVars.add(vars[i]);
            auxintvars.add(auxvars[i]);
         }
      }
   nIntVars = intVars.getSize();
   nBinVars = 0;
   isBin.resize(nIntVars+1);
   for( int j = 0; j < nIntVars; j++ )
   {
      isBin[j] = (intVars[j].getUB() - intVars[j].getLB()) == 1;
      nBinVars += isBin[j];
   }
}



void my_FP::cplexSolve()
{
   cplex.solve();
   cplex.getValues(relaxedIntVars, intVars);
}

double my_FP::round()
{
   double totinf = 0;
   changed = 0;
   double thrs;
   int etrsrnd = trsrnd; //effective randomization used
   if( etrsrnd == 4 )
      etrsrnd = (restarts-s1restarts<5) ? 0 : rng.getInt(4); //if auto, use standard initially, then choose randomly
   if( etrsrnd )
   {
      thrs = rand() / (((double)RAND_MAX)+1);
      if( etrsrnd > 2 )
      {
         /* random quadratic threshold */
         if( thrs <= 0.5)
            thrs = 2*thrs*(1-thrs);
         else
            thrs = 2*thrs*(thrs-1) + 1;
      }
      else {/* random threshold */
         if( etrsrnd == 2 )
            thrs = thrs/2+.25;
      }
   }
   else
      thrs = 0.5; /* standard */



   for( int i = 0; i < nIntVars; i++ )
   {
      if( stage > 1 || isBin[i] )
      {
         double r = floor(relaxedIntVars[i]+thrs);
         if( r < intVars[i].getLB() || r > intVars[i].getUB())
            throw string("relax out of bound");
         double inf = fabs(relaxedIntVars[i]-r);
         totinf += inf;
         if( roundedIntVars[i] != r)
         {
            roundedIntVars[i] = r;
            changed++;
         }
      }
      else
         roundedIntVars[i] = floor(relaxedIntVars[i]+0.5);
   }
   return totinf;
}
