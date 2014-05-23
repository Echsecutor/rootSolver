/**
 * @file exampleSelfconsistencyEquations.cpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
* @version 1.0.2014-05-22
 *
 *
 * @section DESCRIPTION
 *
 *
 * This is a real life use case for the rootSolver framework. We solve
 * the selfconsistency equations encountered in recent research on
 * disordered bosons at the Institute for Theoretical Physics,
 * University of Cologne, in the group of Prof. M. R. Zirnbauer.
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

#include "preProDebugFlags.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>
#include <ctime>


//root finding:
#include "singleRootSolver.hpp"
#include "batchSolver.hpp"
#include "extrapolationSolver.hpp"

//self consistency equations:
#include SCE.hpp

using namespace std;


void readFile(const char * confFile,
              double &bNuGoal,
              double &initialDb,
              double &alpha,
              double &epsilon,
              double &precisionGoal,
              int &desiredNrSteps,
              double &maxDb,
	      double &fromOmega,
	      double &toOmega,
	      double &maxDOmega,
	      int &NThreads,
	      double &extThreshold,
	      bool &onlySaveGoal,
	      double &minDb,
	      double &dim,
	      bool &logStep,
	      double &bsbm
              ){
  cout << "Reading parameters from " << confFile << endl;

  ifstream conf(confFile);
  string line;
  if (!conf.is_open()){
    cerr<< "ERR could not open conf File" <<endl;
    exit(26354);
  }

  char * para;
  char * val;

  while ( conf.good() ){
    getline(conf,line);
    if( line.size()>0 && line.at(0)!='#' ){//ignore empty/comment lines
      char * cstr = new char [line.size()+1];
      strcpy (cstr, line.c_str());
      para = strtok(cstr,"=");
      val = strtok(NULL,"=");

      if(strcmp(para,"bNuGoal")==0){
        bNuGoal=atof(val);
      }else if(strcmp(para,"initialDb")==0){
        initialDb=atof(val);
      }else if(strcmp(para,"alpha")==0){
        alpha=atof(val);
      }else if(strcmp(para,"bsbm")==0){
        bsbm=atof(val);
      }else if(strcmp(para,"epsilon")==0){
        epsilon=atof(val);
      }else if(strcmp(para,"precisionGoal")==0){
        precisionGoal=atof(val);
      }else if(strcmp(para,"desiredNrSteps")==0){
        desiredNrSteps=atoi(val);
      }else if(strcmp(para,"maxDb")==0){
        maxDb=atof(val);
      }else if(strcmp(para,"fromOmega")==0){
        fromOmega=atof(val);
      }else if(strcmp(para,"toOmega")==0){
        toOmega=atof(val);
      }else if(strcmp(para,"maxDOmega")==0){
        maxDOmega=atof(val);
      }else if(strcmp(para,"NThreads")==0){
        NThreads=atoi(val);
      }else if(strcmp(para,"extThreshold")==0){
        extThreshold=atof(val);
      }else if(strcmp(para,"onlySaveGoal")==0){
        onlySaveGoal=(bool)atoi(val);
      }else if(strcmp(para,"minDb")==0){
        minDb=atof(val);
      }else if(strcmp(para,"dim")==0){
        dim=(double)atoi(val);
      }else if(strcmp(para,"logStep")==0){
        logStep=(bool)atoi(val);
      }else{
        cerr << "ERR: CNF: Unknown parameter: " << line <<endl;
      }
      delete cstr;
    }
  }

}



int childProc(const char *inDatFile, const char * oDatFile, double alpha, double bsbm, double epsilon, double omega, double Db, int extPoints, int extOrder, double extThreshold, double precisionGoal, int desiredNrSteps, double maxDb, double bGoal, double minDb,bool onlySaveGoal, double dim, mt19937 * rng=0){

  Matrix<complex<double>,2,1> x;

  data * dat = new data(9);

  diag_params p;
  p.d=dim;
  p.alpha=alpha;
  p.z=epsilon+I*omega;
  p.bsbm=bsbm;

  ofstream oFile(oDatFile);
  if(!oFile.good()){
    cerr<< "ERR: CHILD: could not open '" << oDatFile << "' for writing."<<endl;
    exit(-42);
  }


  oFile << "#\n# This is a data file containing numerical solutions of the diagrammatic self consistency equations.\n";
  oFile << "#It contains data for omega=" << omega <<endl;
  oFile<< "\n# The format is:"<<endl;
  oFile <<"#0     1       2       3      4       5         6        7        8\n";
  oFile <<"#dim   alpha   b/nu    Re(z)  Im(z)   Re(nu g)  Im(nug)  Re(dBM)  Im(dBM)\n"<<endl;


  selfconsistency_Eqs * F;

  F = new selfconsistency_EqsApprox(&p);

  dataInit(inDatFile, dat, F, x, Db, extPoints, extOrder,extThreshold, precisionGoal);
  rootFinding(dat, F, x, Db, extPoints, extOrder,extThreshold, precisionGoal, desiredNrSteps, maxDb, bGoal, oFile,minDb,onlySaveGoal,rng);

  oFile.close();
  delete F;

  return 0;
}


void help(){
  cout << "batchSolver [-h] [-r INPUT_DATA_FILE] ([-c] CONF_FILE) [[-o] DATA_FILE_NAME_PREFIX]"<<endl;
  cout << "If -c and -o are omitted, the conf file has to be given first."<<endl;
  cout <<endl;
  exit(0);
}

////////////////////////////////////////////////////////////////////////////////
/// main ;)
///
int main(int args, char *arg[]){
  char * cnfFile=0;
  char * inDatFile=0;
  const char * outDatFile=0;


  //------------------------------------------------------------------------------------------
  // reading parameters and conf file
  //------------------------------------------------------------------------------------------

  for (int i=1;i<args;i++){
    if(strcmp(arg[i],"-h")==0||strcmp(arg[1],"--help")==0){
      help();
    }else if(strcmp(arg[i],"-c")==0){
      cnfFile=arg[i+1];
      i++;
    }else if(strcmp(arg[i],"-r")==0){
      inDatFile=arg[i+1];
      i++;
    }else if(strcmp(arg[i],"-o")==0){
      outDatFile=arg[i+1];
      i++;
    }else if (arg[i][0]=='-'){
      cout << "MAIN: Unrecognized option " << arg[i] <<endl;
    }else if (cnfFile==0){
      cnfFile=arg[i];
    }else{
      outDatFile=arg[i];
    }
  }

  if(cnfFile==0){
    help();
  }

  string outFilePrefix="T-inv-saddle-points__";

  if(outDatFile!=0){
    outFilePrefix = outDatFile;
  }

  double dim=1;
  double alpha = 1;
  double epsilon = 1e-8;
  double bsbm = 1.0;

  double fromOmega = 1e-2;
  double toOmega = 2;
  double maxDOmega = 1e-1;///< sets a minimal resolution in omega for the resulting data points

  bool logStep=false;///< toggle wether omega points are equally spaced on a linear or a log scale

  // double minDOmega = 1e-5;///< sets a maximal resolution in omega. Important to drop refinement around divergencies at some point.
  // int omegaIntPoints = 10;///< number of points for fittin
  // int omegaIntDegree = 3;///< degree of interpolating polynomials
  // double omegaMaxIntErr = 1e7;///< maximal error for the interpolation. Decrease to enforce smoother data by computing more intermediate points


  double bNuGoal=1;///< will extrapolate from b=0 or from file towards bNuGoal
  double Db=1e-8;///< initial b step size, dynamically adapted
  double maxDb=1e-2;///< sets a minimal resolution in b 

  double minDb=1e-20;///< avoid getting stuck

  double precisionGoal=1e-12;///< for MRS. A to high precuisionGoal will ruine the extrapolation and actually slow down the solver. So, be ambitious! ;)
  int desiredNrSteps=64;///< for the MRS solver, if more/less steps are needed to find the root, dB is adjusted

  double extThreshold=1;///< to determine wether extrapolation was reasonable
  int extPoints=7;///< for extrapolation in b
  int extOrder=3;///< for extrapolation in b

  int NThreads=8;

  bool onlySaveGoal=true;


  readFile(cnfFile,bNuGoal,Db,alpha,epsilon,precisionGoal,desiredNrSteps,maxDb,fromOmega,toOmega,maxDOmega,NThreads,extThreshold,onlySaveGoal,minDb,dim,logStep,bsbm);


  cout << "MAIN: Batch solving saddle point equations for"<<endl;
  cout << "dim = " << dim <<endl;
  cout << "alpha = " << alpha <<endl;
  cout << "bsbm = " << bsbm <<endl;
  cout << "epsilon = " << epsilon <<endl;
  cout << "fromOmega = " << fromOmega <<endl;
  cout << "toOmega = " << toOmega <<endl;
  cout << "maxDOmega = " << maxDOmega <<endl;
  cout << "logStep = " << logStep <<endl;
  cout << "bNuGoal = " << bNuGoal <<endl;
  cout << "Db = " << Db <<endl;
  cout << "precisionGoal = " << precisionGoal<<endl;
  cout << "desiredNrSteps = " << desiredNrSteps <<endl;
  cout << "maxDb = " << maxDb <<endl;
  cout << "NThreads = " << NThreads <<endl;
  cout << "extThreshold = " << extThreshold <<endl;
  cout << "onlySaveGoal = " << onlySaveGoal <<endl;

  cout <<endl;


  if (precisionGoal<1e-12){
    cout << "MAIN: CAUTION! Precision goal is very small. Function evaluation should not be expected to be precise to more than 1e-12 when using the exact formula involving square roots!"<<endl;
  }

  mt19937 * rng=0;


  //------------------------------------------------------------------------------------------
  // batch processing
  //------------------------------------------------------------------------------------------

#if MULTITHREADED
  if(NThreads>1){
    cout <<endl << "MAIN: Start launching parallel threads."<<endl;
  }else{
    cout << "MAIN: You should rather compile the single threaded version than just setting NThreads=1."<<endl;
  }

  future<int>* t = new future<int>[NThreads];
  int cT=0;

#else
  if(NThreads>1){
    cout <<endl << "WARN: MAIN: This programm was compiled as the single threaded Version. NThreads is ignored."<<endl;
  }
#endif

  int fileNbr=0;
  //create initial grid given by maxDOmega
  for(double omega=fromOmega;omega<=toOmega;omega+=maxDOmega){
    stringstream ss;
    ss << outFilePrefix;
    ss << "__part_";
    ss << setfill('0')<< setw(5) << fileNbr;
    ss <<".dat";
    fileNbr++;

    char * cstr = new char [ss.str().length()+1];
    strcpy (cstr, ss.str().c_str());

#if MULTITHREADED
    t[cT] = std::async(launch::async, childProc,inDatFile, cstr, alpha, bsbm, epsilon, omega, Db, extPoints, extOrder, extThreshold, precisionGoal, desiredNrSteps, maxDb, bNuGoal, minDb,onlySaveGoal,dim,rng);
#if DEBUG>=STATUS
    cout << "MAIN: Launched new thread for omega = " << omega << " in slot " << cT <<endl;
#endif
    //search next free thread
    cT=-1;
    while(cT<0){
      for(int i=0;i<NThreads && cT<0;i++){
	if(!t[i].valid()){
	  cT=i;
	}else{
	  // Taken that http://en.cppreference.com/w/cpp/thread/future/wait_for is correct, there seems to be a bug in gcc 4.6.3 where wait_for returns bool
	  // workaround:
#if __GNUC__ == 4 && __GNUC_MINOR__ == 6
	  bool status =true;
#else 
	  //this works with gcc 4.7.2. No other versions tested.
	  future_status status = future_status::ready;
#endif
	  if(status == t[i].wait_for(std::chrono::seconds(1))){ 
	    cT=i;
	  }
	}
      }
    }//end thread search

#else //not MULTITHREADED
#if DEBUG>=STATUS
    cout << "MAIN: Start working on omega = " << omega <<endl;
#endif
    childProc(inDatFile, ss.str().c_str(), alpha, bsbm, epsilon, omega, Db, extPoints, extOrder, extThreshold, precisionGoal, desiredNrSteps, maxDb, bNuGoal,minDb,onlySaveGoal,dim,rng);
#endif //MULTITHREADED

    /// for log step note that \f$ x \frac{(x+Dx)^2}{x^2} = (x + Dx) + (Dx \frac{Dx^2}{x}) \f$
    if(logStep){
      maxDOmega += maxDOmega*maxDOmega/omega;
    }

  }//next omega


#if MULTITHREADED
  //wait for remaining threads to terminate:
  for(int i=0;i<NThreads && cT<0;i++){
    if(t[i].valid()){
      t[i].wait();
    }
  }
#endif
  
  cout << "Collecting results using cat and deleting part files"<<endl;

  stringstream cmd;
  cmd << "cat " <<outFilePrefix << "__part_*.dat > " << outFilePrefix << "__collected.dat;";
  cmd << "rm "<<outFilePrefix << "__part_*.dat";
  system(cmd.str().c_str());

  delete rng;

  cout <<endl << "MAIN: FINISHED"<<endl<<endl;
  return 0;
}
