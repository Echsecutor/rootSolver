/**
 * @file testSRS.cpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0.2014-06-11
 *
 * @section DESCRIPTION
 *
 * Testcase generator for the single root solver.
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
#include <string.h>
#include <time.h>

#include <cstdlib>//rand


using namespace std;

//actual root finding:
#include "singleRootSolver.hpp"

//trial function
#include "polynomials.hpp"

#include "complex.hpp"
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

  realT epsilon=1e-10;

  int seed=42;
  int maxTotalDegree=42;

  for (int i=1;i<args;i++){
    if(strcmp(arg[i],"-h")==0||strcmp(arg[1],"--help")==0){
      help();
    }else if(strcmp(arg[i],"-r")==0){
      cout << __FILE__ << " : Seeding RNG with time, hence producing 'unique' testcase"<<endl;
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
      cout << "use -h for help"<<endl;
    }
  }

  cout << __FILE__ << " : Seeding RNG with seed="<<seed<<endl;
  srand(seed);


  polyParams<realT>* para = new polyParams<realT>();
  para->maxDegree = maxTotalDegree;
  polynomial<realT> * F = new polynomial<realT>(para);
  F->randomiseTestcase();

  singleRootSolver<realT> SRS(F);

  SRS.setStartPoint(F->guessStartPoint());

  int step=1;
  SRS.setPrecisionGoal(epsilon);
  solver_state state=CONTINUE;

  while (state >= CONTINUE){
    state=SRS.step();
#if DEBUG>=DETAIL
    cout << __FILE__ << " : After " << step << " iterations the solver achieved |f| = " << SRS.getAbsF() << " state = " << endl<<endl;
#endif
    step++;
  }
#if DEBUG>=STATUS
  cout << __FILE__ << " : After " << step << " iterations the solver achieved |f| = " << SRS.getAbsF() <<endl<<endl;
#endif

  if (SRS.getState()==SUCCESS){
    cout << __FILE__ << " : The solver found an approximate root f(z) =" << SRS.getLastValue() << " \nat z =" << SRS.getLastPoint()<<endl;
  }else{
    cout << __FILE__ << " : The solver got stuck after " << step << " iterations at z = " << SRS.getLastPoint() << "\nwhere f(z) =" << SRS.getLastValue()<<endl;
  }

  cout << __FILE__ << " : The solver used " << F->getFunctionCallsCounter() << " function calls and computed the jacobian " << F->getJacobianCallsCounter() << " times."<<endl;


  // complex<realT> z=SRS.getLastPoint();
  // int closest=0;
  // realT dist =abs(z-zeros[0]);

  // for(int i=0;i<degree;i++){
  //   if(abs(zeros[i]-z)<dist){
  //     dist = abs(zeros[i]-z);
  //     closest=i;
  //   }
  // }
  // cout << "\nThe closest actual zero is " << zeros[closest] << " with distance " << dist <<endl;


  return 0;
}
