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


using namespace std;


#include "myfp.h"


void my_FP::setNewObj()
{
    foconst=0;
    IloExpr expr ( env );
    double expr_norm = 0.0;
    for ( int i = 0; i < nIntVars; i++ ) {
        if ( isBin[i] ) {

            if ( roundedIntVars[i] == intVars[i].getLB() ) {
                expr += intVars[i];
                foconst -= intVars[i].getLB();
            } else if ( roundedIntVars[i] == intVars[i].getUB() ) {
                expr -= intVars[i];
                foconst += intVars[i].getUB();
            } else {
                throw string ( "Hei, we have a binary var that, once rounded, is not at a bound. This shouldn't happen!" );
            }
            expr_norm += 1.0;
        }
    }
    dist = expr;
    expr = ( 1 - alpha ) * expr + alpha * sqrt ( expr_norm ) * origObjExpr / origobj_norm;
    expr.normalize();
    activeObj.setExpr ( expr );
    activeObj.setSense ( IloObjective::Minimize );
    expr.end();
}

void my_FP::getNextIntegerPt()
{
    if ( !changed ) {
        priority_queue<dipair, vector<dipair>, dipairCmp> q;
        unsigned int toBeChanged= ( minChange+rng.getInt ( minChange+minChange ) ) /2;
        // populate queue with top <toBeChanged> biggest sigma (and sigma>trsld)
        double sigma, min_sigma=trsmin;
        for ( int i = 0; i < nIntVars; i++ ) {
            if ( isBin[i] ) {
                sigma=fabs ( roundedIntVars[i]-relaxedIntVars[i] );
                if ( sigma > min_sigma || min_sigma == 0 ) {
                    q.push ( dipair ( sigma, i ) );
                    if ( q.size() > toBeChanged ) {
                        q.pop();
                        min_sigma=q.top().first;
                    }
                }
            }
        }
        // change the vars in the queue
        while ( !q.empty() ) {
            int i = q.top().second;
            if ( roundedIntVars[i] > intVars[i].getLB() + 0.5 ) {
                roundedIntVars[i]--;
            } else {
                roundedIntVars[i]++;
            }
            changed++;
            q.pop();
        }
    }
}

bool my_FP::isCoincident()
{
    for ( int i = 0; i < nIntVars; i++ )
        if ( isBin[i] && fabs ( relaxedIntVars[i]-roundedIntVars[i] ) > epInt ) {
            return false;
        }
    return true;
}

void my_FP::restart(){
    restarts++;
    for ( int i = 0; i < nIntVars; i++ ) {
        if ( isBin[i] ) {
            double r=rng.getFloat()-.47;
// 		 r = -0.03;
            if ( r > 0 && prevRoundedIntVars[i] == roundedIntVars[i] ) {
                double sigma=fabs ( roundedIntVars[i]-relaxedIntVars[i] );
                if ( sigma + r > 0.5 ) {
                    if ( roundedIntVars[i] == intVars[i].getUB() ) {
                        roundedIntVars[i]--;
                    } else {
                        roundedIntVars[i]++;
                    }
                    changed++;
                }
            }
        }
    }
}


// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
