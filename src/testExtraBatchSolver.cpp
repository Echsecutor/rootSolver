/**
 * @file testExtraBatchSolver.cpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0.2014-06-11
 *
 * @section DESCRIPTION
 *
 * Testcase for the complete batch + extrapolation solver. Testing MRS
 * and SRS.
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


#ifndef MULTITHREADED
#define MULTITHREADED 1
#endif

#include "preProDebugFlags.h"

#include <iostream>
#include <string>
#include <time.h>

#include <cstdlib>//rand

using namespace std;

//trial function
#include "polynomials.hpp"

//actual root finding:
#include "multiRootSolver.hpp"
#include "singleRootSolver.hpp"
#include "batchSolver.hpp"
#include "extrapolationSolver.hpp"

using namespace root_solver;


void help(){
  cout << "\nusage: test [-h] [-r] [-s SEED] [-m MAX_TOTAL_DEGREE] [-p PRECISION] [-b RUNS] [-t THREADS] [-o FILENAME]"<<endl;
  cout<<"-h\tthis help\n";
  cout<<"-r\trandomise using time\n";
  cout<<"-s\tseed the random number generator with SEED to reproduce a known test case (default 42)\n";
  cout<<"-m\tset the maximal degree of polynomials used to MAX_TOTAL_DEGREE (default 9)\n";
  cout<<"-p\tset the precision goal to PRECISION, i.e. ||f||<PRECISION is counted as a root (default 1e-8)\n";
  cout<<"-b\tset the number of batch runs to perform to RUNS\n";
  cout<<"-t\tset the number of (parallel) threads to use to THREADS\n";
  cout<<"-o\twrite dat output to FILENAME\n";
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

  int batchRuns=64;
  int nThreads=4;

  string outFileName="Extra-batch-test";

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
    }else if(strcmp(arg[i],"-b")==0){
      batchRuns=atoi(arg[i+1]);
      cout << __FILE__ << " : Set batchRuns "<<batchRuns<<endl;
      i++;
    }else if(strcmp(arg[i],"-t")==0){
      nThreads=atoi(arg[i+1]);
      cout << __FILE__ << " : Set nThreads "<< nThreads <<endl;
      i++;
    }else if(strcmp(arg[i],"-o")==0){
      outFileName = arg[i+1];
      cout << __FILE__ << " : Set outFileName "<< outFileName <<endl;
      i++;
    }else{
      cout << __FILE__ << " : Unrecognized option " << arg[i] <<endl;
    }
  }

  cout << __FILE__ << " : Seeding RNG with seed="<<seed<<endl;
  srand(seed);



  cout <<endl << __FILE__ << " : Starting Multi Root Solver Extrapolation Batch Test."<<endl<<endl;


  polyParams<realT> para;
  para.maxDegree = maxTotalDegree;
  para.dim=dim;
  PSI<realT> it(para);
  ++it;//randomise
  it.counter=0;//reset counter
  PSI<realT> end(para);
  end.counter=batchRuns;
  polynomials<dim,realT> *F = new polynomials<dim,realT>(*it);


  typedef multiRootSolver<root_solver::complex<realT>,dim> MRS;
  typedef extraSolver<MRS,realT> extraMRS;
  typedef batchSolver<extraMRS, PSI<realT> > batchExtraMRS;

  batchExtraMRS::run(nThreads,  F, it, end, outFileName + "-MRS","#Solutions of some random multi-dimensional polynomial equations\n", epsilon);

  cout << __FILE__ << " : Finished."<<endl<<endl;

  cout <<endl << __FILE__ << " : Starting Single Root Solver Batch Test."<<endl<<endl;

  typedef singleRootSolver<realT > SRS;
  typedef extraSolver<SRS,realT> extraSRS;
  typedef batchSolver<extraSRS, PSI< realT> > batchExtraSRS;


  para.maxDegree = maxTotalDegree;
  para.dim=1;
  PSI< realT> itOne(para);
  ++itOne;//randomise
  itOne.counter=0;//reset counter
  PSI<realT> endOne(para);
  endOne.counter=batchRuns;
  polynomial<realT> *FOne = new polynomial<realT>(*itOne);


  batchExtraSRS::run(nThreads, FOne, itOne, endOne, outFileName + "-SRS","#Solutions of some random polynomial equation\n", epsilon);



  return 0;
}
