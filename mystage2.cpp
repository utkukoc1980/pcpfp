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

// #include <iostream>
// #include <sstream>
// #include <iomanip>
// #include <vector>
// #include <queue>
// #include <cmath>

using namespace std;


#include "myfp.h"


void my_FP::introduceDelta(int i)
{
   nDeltaVars++;
   deltaVars[i] = new IloNumVar(env, 0.0, IloInfinity, IloNumVar::Float, (string("delta-")+getVarName(intVars[i])).c_str());
   model.add(*deltaVars[i]);
   vars.add(*deltaVars[i]);
   deltaCnstrUB[i] = new IloRange(env, -IloInfinity, intVars[i]-*deltaVars[i], IloInfinity);
   model.add(*deltaCnstrUB[i]);
   rngs.add(*deltaCnstrUB[i]);
   deltaCnstrLB[i] = new IloRange(env, -IloInfinity, *deltaVars[i]+intVars[i], IloInfinity);
   model.add(*deltaCnstrLB[i]);
   rngs.add(*deltaCnstrLB[i]);
}

void my_FP::setNewObjS2()
{
   foconst = 0;
   IloExpr expr(env);
   double expr_norm = 0.0;
   for( int i = 0; i < nIntVars; i++ )
   {
      if( !deltaVars[i] )
      {
         if( roundedIntVars[i] == intVars[i].getLB() )
         {
            expr += intVars[i];
            foconst -= intVars[i].getLB();
         }
         else if( roundedIntVars[i] == intVars[i].getUB() )
         {
            expr -= intVars[i];
            foconst += intVars[i].getUB();
         }
         else
            introduceDelta(i);
      }
      if( deltaVars[i] )
      {
         expr += (*deltaVars[i]);
         deltaCnstrUB[i]->setUB(roundedIntVars[i]);
         deltaCnstrLB[i]->setLB(roundedIntVars[i]);
      }
      expr_norm += 1.0;

   }
   dist = expr;
   if( stage != 3)
   {
      expr = (1 - alpha) * expr + alpha * sqrt(expr_norm) * origObjExpr / origobj_norm;
      expr.normalize();
   }

   activeObj.setExpr(expr);
   activeObj.setSense(IloObjective::Minimize);
   expr.end();
}


void my_FP::getNextIntegerPtS2()
{
   if( !changed )
   {
      priority_queue<dipair, vector<dipair>, dipairCmp> q;
      double sigma, min_sigma=trsmin;
      // populate queue with top <toBeChanged> biggest sigma (and sigma>trsld)
      unsigned int toBeChanged=(minChange+rng.getInt(minChange+minChange))/2;
      for( int i = 0; i < nIntVars; i++ )
      {
         sigma=fabs(roundedIntVars[i]-relaxedIntVars[i]);
         if( sigma > min_sigma )
         {
            q.push(dipair(sigma, i));
            if( q.size() > toBeChanged)
            {
               q.pop();
               min_sigma=q.top().first;
            }
         }
      }
      // change the vars in the queue
      while( !q.empty() )
      {
         int i = q.top().second;
         if( roundedIntVars[i] > relaxedIntVars[i] )
            roundedIntVars[i]--;
         else roundedIntVars[i]++;
         changed++;
         q.pop();
      }
   }
}

bool my_FP::isCoincidentS2()
{
   for( int i = 0; i < nIntVars; i++ )
      if( fabs(relaxedIntVars[i]-roundedIntVars[i]) > epInt)
         return false;
   return true;
}

void my_FP::restartS2()
{
   static int lastiter, rp;
   restarts++;

   if( lastiter < iter)
   {
      while( ++lastiter < iter )
         rp=(int)(rp*.85);
   }
   rp = min(nIntVars/10, 2*minChange+rp+1);
   for( int j = 0; j < rp; j++)
   {
      long mlh = 1<<30;
      int i = rng.getInt(nIntVars);
      double r;
      if( intVars[i].getUB()-intVars[i].getLB() < 1e15 )
         r = floor(intVars[i].getLB() + ( 1 + intVars[i].getUB() - intVars[i].getLB() ) * rng.getFloat());
      else if( roundedIntVars[i]-intVars[i].getLB() < mlh )
         r = intVars[i].getLB() + rng.getInt( -1+mlh+mlh );
      else if( -roundedIntVars[i]+intVars[i].getUB() < mlh )
         r = intVars[i].getUB() - rng.getInt(-1+mlh+mlh);
      else
         r = roundedIntVars[i]+rng.getInt(-1+mlh+mlh)-mlh;
      if( r != roundedIntVars[i] )
      {
         changed++;
         roundedIntVars[i] = r;
      }
   }
}
