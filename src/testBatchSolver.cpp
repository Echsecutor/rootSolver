/**
 * @file testBatchSolver.cpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0.2014-06-11
 *
 * @section DESCRIPTION
 *
 * Testcase generator for the batch solver. Testing MRS and SRS batches.
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

#ifndef MULTITHREADED
#define MULTITHREADED 1
#endif

#include "preProDebugFlags.h"

#include <iostream>
#include <string>
#include <time.h>

#include <cstdlib>//rand

using namespace std;

//actual root finding:
#include "multiRootSolver.hpp"
#include "singleRootSolver.hpp"
#include "batchSolver.hpp"

using namespace root_solver;

//trial function
#include "polynomials.hpp"


#include "complexMprealWrapper.hpp"
using namespace mpfr;
typedef mpreal realT;

//  typedef double realT;




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

  initScalarType<realT>::ini();

  const int dim=3;

  realT epsilon=1e-8;

  int seed=42;
  int maxTotalDegree=9;

  int batchRuns=64;
  int nThreads=4;

  string outFileName="batch-test";

  for (int i=1;i<args;i++){
    if(strcmp(arg[i],"-h")==0||strcmp(arg[1],"--help")==0){
      help();
    }else if(strcmp(arg[i],"-r")==0){
      cout << "TEST: Seeding RNG with time, hence producing 'unique' testcase"<<endl;
      seed=time(NULL) ^ clock();
      //      cout << "TEST: seed="<<seed<<endl;
    }else if(strcmp(arg[i],"-s")==0){
      seed=atoi(arg[i+1]);
      i++;
    }else if(strcmp(arg[i],"-m")==0){
      maxTotalDegree=atoi(arg[i+1]);
      cout << "TEST: Set maximum degree of polynomials to "<<maxTotalDegree<<endl;
      i++;
    }else if(strcmp(arg[i],"-p")==0){
      epsilon=atof(arg[i+1]);
      cout << "TEST: Set precision to "<<epsilon<<endl;
      i++;
    }else if(strcmp(arg[i],"-b")==0){
      batchRuns=atoi(arg[i+1]);
      cout << "TEST: Set batchRuns "<<batchRuns<<endl;
      i++;
    }else if(strcmp(arg[i],"-t")==0){
      nThreads=atoi(arg[i+1]);
      cout << "TEST: Set nThreads "<< nThreads <<endl;
      i++;
    }else if(strcmp(arg[i],"-o")==0){
      outFileName = arg[i+1];
      cout << "TEST: Set outFileName "<< outFileName <<endl;
      i++;
    }else{
      cout << "TEST: Unrecognized option " << arg[i] <<endl;
    }
  }

  cout << "TEST: Seeding RNG with seed="<<seed<<endl;
  srand(seed);



  cout <<endl << "Starting Multi Root Solver Batch Test."<<endl<<endl;


  polyParams<realT> para(dim,maxTotalDegree);
  PSI<realT> it(&para);
  ++it;//randomise
  it.counter=0;//reset counter
  PSI<realT> end(&para);
  end.counter=batchRuns;
  polynomials<dim,realT> *F = new polynomials<dim,realT>((*it));

  typedef batchSolver<multiRootSolver<root_solver::complex<realT>,dim >, PSI<realT> > batchMRS;

  batchMRS::run(nThreads, F, it, end, outFileName + "-MRS", "#Solutions of some random multi-dimensional polynomial equations\n", epsilon);

  cout << "Finished."<<endl<<endl;


  cout <<endl << "Starting Single Root Solver Batch Test."<<endl<<endl;

  para.dim=1;
  para.maxDegree = maxTotalDegree;
  para.randomiseTestcase();
  PSI<realT> itOne(&para);
  PSI<realT> endOne(&para);
  endOne.counter=batchRuns;
  polynomial<realT> *FOne = new polynomial<realT>(para);

  typedef batchSolver<singleRootSolver<realT>, PSI<realT> > batchSRS;

  batchSRS::run(nThreads, FOne, itOne, endOne, outFileName + "-SRS","#Solutions of some random polynomial equation\n", epsilon);



  return 0;
}
