/**
 * @file batchSolver.hpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0.2014-04-30
 *
 * @section DESCRIPTION
 *
 * This is a higher level interface for the rootSolver for the following scenario:
 * We are looking for solutions of the (set of) equations \f$ f(z)=0\f$ for some
 * \f$ f_p \mathbb{C}^d \to \mathbb{C}^d \f$ where $p$ denotes an set of "parameters".
 * We want to get the solution for all parameters iterating through some user defined parameter container.
 * The resulting implizit functions \f$ z(p) \f$ is formatted such as to be written to a dat file.
 *
 * This wrapper may be compiled to use multi threading. No interdependence of \f$ z(p_i)\f$
 * and \f$z(p_j)\f$ is used, hence parallelisation is trivial.
 *
 *
 *
 * @section LICENSE
 *
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


#ifndef BATCHSOLVER_HPP
#define BATCHSOLVER_HPP

#include "preProDebugFlags.h"


#ifndef MULTITHREADED
#define MULTITHREADED 1
#endif

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>

#if MULTITHREADED
#include <future>
#endif

#include <chrono>

#include "rootSolver.hpp"


using namespace std;



//------------------------------------------------------------------------------------------
// Declarations:
//------------------------------------------------------------------------------------------


template <typename rootSolverT, class parameterSetIterator>
class batch_functions : public virtual rootSolverT::functions_type{//be amazed how templates avoid inheritance diamonds here ;)
 public:
  virtual void setParameters(parameterSetIterator It)=0;
};


/** /brief static class to wrap the batch solver.
 *
 * (Not that I like static classes... Tell me if you know a nice alternative compatible with threading!)
 *
 *
 * You need to implement/provide the parameter iterator (likely iterating through your container).
 *
 * The parameterSetIterator should move through your desired points in parameter space upon increment.
 *
 * Please make sure that *(parameterSetIterator) implements the << operator to
 * produced a .dat format string, i.e. "p[1]\tp[2]\t...\tp[nParameters]"
 *
 *
 * rootSolverT shall be any of the solvers implementing the rootSolver<valueT,derivativeT> interface. 
 * (I.e. currently Single or Multi Root Solver.)
 *
 *
 */
template <typename rootSolverT, class parameterSetIterator>
class batchSolver{

private:
  //=>do not create instances!
  batchSolver();
  ~batchSolver();
  batchSolver(const batchSolver&);
  batchSolver& operator=(const batchSolver&);


  static int childProc(batch_functions<rootSolverT,parameterSetIterator> *F, parameterSetIterator It, string outFileName, double precisionGoal);

public:

  static void run(int NThreads, batch_functions<rootSolverT,parameterSetIterator> *F, parameterSetIterator It, parameterSetIterator End, string outFilePrefix, double precisionGoal);

};




//------------------------------------------------------------------------------------------
// Implementation:
//------------------------------------------------------------------------------------------

template<typename rootSolverT, class parameterSetIterator>
int batchSolver<rootSolverT, parameterSetIterator>::childProc(batch_functions<rootSolverT,parameterSetIterator> *F, parameterSetIterator It, string outFileName, double precisionGoal){

  ofstream out;

  try {
    out.open(outFileName.c_str());
  }
  catch (std::ios_base::failure &e) {
#if DEBUG>=ERR
    cerr << "ERR: CHILD: I could not open " << outFileName << " for writing." <<endl;
#endif
    throw e;
  }

  //run the root finder:
  F->setParameters(It);
  rootSolverT Solver(F);
  Solver.setStartPoint(F->guessStartPoint());

#if DEBUG>=STATUS
  cout << __FILE__ << " : Solving saddle point equations for z = " <<  Solver.getLastPoint()<< " at Parameters = " << *It <<endl;
#endif

  int steps = 0;
  solver_state state = CONTINUE;

  while (state == CONTINUE){
    state = Solver.step(precisionGoal);
#if DEBUG>=DETAIL
    steps++;
    cout << __FILE__ << " : After " << steps << " iterations the solver achieved |f| = " << Solver.getAbsF() << ". State = " << state <<endl;
#endif
  }

#if DEBUG>=DETAIL
  if (Solver.getState()==SUCCESS){
    cout << __FILE__ << " : The solver found an approximate root. ";
    cout << "at Parameters = " << *It <<endl<<endl;
  }else{
    cout << __FILE__ << " : The solver got stuck after " << steps << " iterations\nat z =\n" << Solver.getLastPoint() << "\nwhere f =\n" << Solver.getLastValue()<<endl;
    cout << __FILE__ << " : Parameters = " << *It <<endl<<endl;
  }
#endif

  if (Solver.getState()==MRS_SUCCESS){
    out << Solver.getLastPointInDatFormat() <<"\t"<<*It<<endl;
  }else{
    out << "#Could not find a root for Parameters = " << *It<<endl;
  }

  out.close();

  return 0;
}



//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------


template <typename rootSolverT, class parameterSetIterator>
void batchSolver<rootSolverT, parameterSetIterator>::run(int NThreads, batch_functions<rootSolverT,parameterSetIterator> *F, parameterSetIterator It, parameterSetIterator End, string outFilePrefix, double precisionGoal){

  //------------------------------------------------------------------------------------------
  // batch processing
  //------------------------------------------------------------------------------------------

#if MULTITHREADED
#if DEBUG>=STATUS
  if(NThreads>1){
    cout <<endl << __FILE__ << " : RUN: Start launching " << NThreads << " parallel threads."<<endl;
  }else{
#if DEBUG>=DETAIL
    cout << __FILE__ << " : RUN: You could compile the single threaded version than setting NThreads=1."<<endl;
#endif
  }
#endif

  future<int>* t = new future<int>[NThreads];
  int cT=0;

#else
#if DEBUG>=WARN
  if(NThreads>1){
    cout <<endl << "WARN: BATCH_SOLVER: RUN: This programm was compiled as the single threaded Version. NThreads is ignored."<<endl;
  }
#endif
#endif

  int fileNbr=0;


  while(It != End){
    stringstream ss;
    ss << outFilePrefix;
    ss << "__part_";
    ss << setfill('0')<< setw(5) << fileNbr;
    ss <<".dat";
    fileNbr++;

    string fileName(ss.str());


#if MULTITHREADED
    int wait =1000;
    t[cT] = std::async(launch::async, childProc, F, It, fileName, precisionGoal);
#if DEBUG>=STATUS
    cout << "Batch_SOLVER RUN: Launched new thread for Parameters = " << *It << " in slot " << cT <<endl;
#endif
    //search next free thread
    cT=-1;
    while(cT<0){
      for(int i=0;i<NThreads && cT<0;i++){
        if(!t[i].valid()){
          cT=i;
	  wait/=2;
        }else{
          // Taken that http://en.cppreference.com/w/cpp/thread/future/wait_for is correct, there seems to be a bug in gcc 4.6.3 where wait_for returns bool
          // workaround:
#if __GNUC__ == 4 && __GNUC_MINOR__ == 6
          bool status =true;
#else
          //this works with gcc 4.7.2. No other versions tested.
          future_status status = future_status::ready;
#endif
          if(status == t[i].wait_for(std::chrono::milliseconds(wait))){
            cT=i;
	    wait/=2;
          }
        }
      }
      if(cT<0)
	wait = min(wait*2, 1000);
    }//wend thread search

#else //not MULTITHREADED
#if DEBUG>=STATUS
    cout << "BATCH_SOLVER RUN: Start working on Parameters = " << *It << " in single threaded mode." << endl;
#endif
    childProc(F,It, fileName, precisionGoal);
#endif //MULTITHREADED

    It++;
  }//next Parameter set


#if MULTITHREADED
  //wait for remaining threads to terminate:
  for(int i=0;i<NThreads && cT<0;i++){
    if(t[i].valid()){
      t[i].wait();
    }
  }
#endif

#if DEBUG>=STATUS
  cout << "BATCH_SOLVER RUN: Collecting results using cat and deleting part files"<<endl;
#endif

  stringstream cmd;
  cmd << "cat " <<outFilePrefix << "__part_*.dat > " << outFilePrefix << "__collected.dat;";
  cmd << "rm "<<outFilePrefix << "__part_*.dat";
  system(cmd.str().c_str());


#if DEBUG>=STATUS
  cout <<endl << "BATCH_SOLVER RUN:: FINISHED"<<endl<<endl;
#endif

}




#endif

