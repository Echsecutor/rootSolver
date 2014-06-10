/**
 * @file testMRS.cpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0.2014-05-21
 *
 * @section DESCRIPTION
 *
 * Testcase generator for the multi root solver.
 *
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License and a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#define DEBUG 11

#include "preProDebugFlags.h"


#include <iostream>
#include <string>
#include <time.h>

#include <cstdlib>//rand

using namespace std;

//actual root finding:
#include "multiRootSolver.hpp"

//trial function
#include "polynomials.hpp"

using namespace root_solver;

void help(){
  cout << "\nusage: test [-h] [-r] [-s SEED] [-m MAX_TOTAL_DEGREE] [-p PRECISION]"<<endl;
  cout<<"-h\tthis help\n";
  cout<<"-r\trandomise using time\n";
  cout<<"-s\tseed the random number generator with SEED to reproduce a known test case (default 42)\n";
  cout<<"-m\tset the maximal degree of polynomials used to MAX_TOTAL_DEGREE (default 9)\n";
  cout<<"-p\tset the precision goal to PRECISION, i.e. ||f||<PRECISION is counted as a root (default 1e-8)\n";
  cout <<endl;

  exit(0);
}


////////////////////////////////////////////////////////////////////////////////
/// main ;)
///
int main(int args, char *arg[]){

  typedef double realT;

  const int dim=3;

  realT epsilon=1e-8;

  int seed=42;
  int maxTotalDegree=9;

  for (int i=1;i<args;i++){
    if(strcmp(arg[i],"-h")==0||strcmp(arg[1],"--help")==0){
      help();
    }else if(strcmp(arg[i],"-r")==0){
      cout << __FILE__<<" : Seeding RNG with time, hence producing 'unique' testcase"<<endl;
      seed=time(NULL) ^ clock();
      //      cout << __FILE__ << " : seed="<<seed<<endl;
    }else if(strcmp(arg[i],"-s")==0){
      seed=atoi(arg[i+1]);
      i++;
    }else if(strcmp(arg[i],"-m")==0){
      maxTotalDegree=atoi(arg[i+1]);
      cout << __FILE__ << " : Set maximum degree of polynomials to "<<maxTotalDegree<<endl;
      i++;
    }else if(strcmp(arg[i],"-p")==0){
      epsilon=atof(arg[i+1]);
      cout << __FILE__ << " : Set precision to "<<epsilon<<endl;
      i++;
    }else{
      cout << __FILE__ << " : Unrecognized option " << arg[i] <<endl;
    }
  }

  cout << __FILE__ << " : Seeding RNG with seed="<<seed<<endl;
  srand(seed);

  polyParams<realT>* para = new polyParams<realT>();
  para->maxDegree = maxTotalDegree;
  polynomials<dim,realT> *F;
  F = new polynomials<dim,realT>(para);
  F->randomiseTestcase();
  multiRootSolver<root_solver::complex<realT>,dim> MRS(F);
  MRS.setStartPoint(F->guessStartPoint());

  int step=1;
  MRS.setPrecisionGoal(epsilon);
  solver_state state=MRS.step();
  while (state >= CONTINUE){
    cout << __FILE__ << " : After " << step << " iterations the solver achieved |f| = " << MRS.getAbsF() << ", state = " << state<<endl<<endl;
    step++;
    state=MRS.step();
  }
  cout << __FILE__ << " : After " << step << " iterations the solver achieved |f| = " << MRS.getAbsF() <<endl<<endl;

  if (MRS.getState()==SUCCESS){
    cout << __FILE__ << " : The solver found an approximate root f(z) =\n" <<MRS.getLastValue()<< " \nat z =\n" << MRS.getLastPoint()<<endl;
  }else{
    cout << __FILE__ << " : The solver got stuck after " << step << " iterations at z = \n" << MRS.getLastPoint() << "\nwhere f(z) =\n" << MRS.getLastValue()<<endl;
    cout << __FILE__ << " : Last step size = " << MRS.getAbsLastStep() << " last abs value change = " << MRS.getLastAbsValueChange()<<endl;
  }

  cout << __FILE__ << " : The solver used " << F->getFunctionCallsCounter() << " function calls and computed the jacobian " << F->getJacobianCallsCounter() << " times."<<endl;

  return 0;
}
