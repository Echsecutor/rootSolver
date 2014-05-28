/**
 * @file selfconsistencyEquations.cpp
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

#include<complex>
using namespace std;

#ifndef COMPLEX_I
#define COMPLEX_I
static complex<double> I(0.0,1.0);
#endif

//#include "logger.hpp"

//root finding:
#include "singleRootSolver.hpp"
#include "batchSolver.hpp"
#include "extrapolationSolver.hpp"

//self consistency equations:
#include "SCE.hpp"



void readFile(const char * confFile,
              double &bmTarget,
              double &initialDb,
              double &epsilon,
              double &precisionGoal,
              int &desiredNrSteps,
              double &maxDb,
	      double &fromOmega,
	      double &toOmega,
	      double &DOmega,
	      int &NThreads,
	      bool &onlySaveGoal,
	      double &minDb,
	      double &dim,
	      bool &logStep,
	      double &bsbm
              ){

  //  logger::write(string("Reading parameters from ") + string(confFile), STATUS, __FILE__, logger::removePath);
#if DEBUG>=STATUS
  cout << __FILE__ << " : Reading parameters from " <<confFile<<endl;
#endif

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

      if(strcmp(para,"bmTarget")==0){
        bmTarget=atof(val);
      }else if(strcmp(para,"initialDb")==0){
        initialDb=atof(val);
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
      }else if(strcmp(para,"DOmega")==0){
        DOmega=atof(val);
      }else if(strcmp(para,"NThreads")==0){
        NThreads=atoi(val);
      }else if(strcmp(para,"onlySaveGoal")==0){
        onlySaveGoal=(bool)atoi(val);
      }else if(strcmp(para,"minDb")==0){
        minDb=atof(val);
      }else if(strcmp(para,"dim")==0){
        dim=(double)atoi(val);
      }else if(strcmp(para,"logStep")==0){
        logStep=(bool)atoi(val);
      }else{
#if DEBUG>=ERR
        cerr << __FILE__ << " : ERR: Unknown parameter: " << line <<endl;
#endif
      }
      delete cstr;
    }
  }

}


void help(){
  cout << "selfconsistencyEquations [-h] [-r INPUT_DATA_FILE] [-c] CONF_FILE [[-o] DATA_FILE_NAME_PREFIX]"<<endl;
  cout << "If -o is omitted, the conf file has to be given first."<<endl;
  cout <<endl;
  exit(0);
}

////////////////////////////////////////////////////////////////////////////////
/// main ;)
///
int main(int args, char *arg[]){

  cout.precision(5);
  cout << std::scientific;

  char * cnfFile=0;
  //  char * inDatFile=0;
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
    // }else if(strcmp(arg[i],"-r")==0){
    //   inDatFile=arg[i+1];
    //   i++;
    }else if(strcmp(arg[i],"-o")==0){
      outDatFile=arg[i+1];
      i++;
    }else if (arg[i][0]=='-'){
      cout << "Unrecognised option " << arg[i] <<endl;
    }else if (cnfFile==0){
      cnfFile=arg[i];
    }else{
      outDatFile=arg[i];
    }
  }

  if(cnfFile==0){
    help();
  }

  string outFilePrefix="SCE-solutions__";

  if(outDatFile!=0){
    outFilePrefix = outDatFile;
  }

  double dim=1;
  double epsilon = 1e-9;
  double bsbm = 1.0;

  double fromOmega = 1e-3;
  double toOmega = 2;
  double DOmega = 1e-1;

  bool logStep=false;///< toggle wether omega points are equally spaced on a linear or a log scale

  double bmTarget=1e-2;///< will extrapolate from b=0 or from file towards bNuGoal
  double Db=1e-8;///< initial b step size, dynamically adapted
  double maxDb=1e-3;///< sets a minimal resolution in b 

  double minDb=1e-10;///< avoid getting stuck

  double precisionGoal=1e-12;
  int desiredNrSteps=32;

  int NThreads=8;

  bool onlySaveGoal=true;


  readFile(cnfFile,bmTarget,Db,epsilon,precisionGoal,desiredNrSteps,maxDb,fromOmega,toOmega,DOmega,NThreads,onlySaveGoal,minDb,dim,logStep,bsbm);

  cout << __FILE__ << " : Start solving saddle point equations for"<<endl;
  cout << "dim = " << dim <<endl;
  cout << "bs / bm = " << bsbm <<endl;
  cout << "epsilon = " << epsilon <<endl;
  cout << "fromOmega = " << fromOmega <<endl;
  cout << "toOmega = " << toOmega <<endl;
  cout << "DOmega = " << DOmega <<endl;
  cout << "logStep = " << logStep <<endl;
  cout << "bmTarget = " << bmTarget <<endl;
  cout << "Db = " << Db <<endl;
  cout << "precisionGoal = " << precisionGoal<<endl;
  cout << "desiredNrSteps = " << desiredNrSteps <<endl;
  cout << "maxDb = " << maxDb <<endl;
  cout << "NThreads = " << NThreads <<endl;
  cout << "onlySaveGoal = " << onlySaveGoal <<endl;
  cout << endl;

  cout << "The parameter format is\n" << SCE_parameters::getFormat() <<endl<<endl;

#if DEBUG>=WARN
  if (precisionGoal < 1e-12){
    cout << __FILE__ << " : CAUTION! Precision goal is very small. Function evaluation should not be expected to be precise to more than 1e-12!" << endl;
  }
#endif


#if DEBUG >= STATUS
  cout << __FILE__ << " : Initialising parameter list..."<<endl;
#endif
  list <SCE_parameters> params;
  double omega=fromOmega;
  while(omega < toOmega){
    params.push_back(SCE_parameters(dim,epsilon + I * omega));

    if(logStep){
      DOmega += DOmega * DOmega / omega;
    }

    omega += DOmega;
  }
  params.push_back(SCE_parameters(dim,epsilon + I * toOmega));

#if DEBUG >= STATUS
  cout << __FILE__ << " : I will work on " << params.size() << " parameter sets"<<endl;
  cout << __FILE__ << " : Preparing self consistency equation."<<endl;
#endif

  SCE * F = new SCE(bsbm, bmTarget, maxDb, minDb, Db);

#if DEBUG >= STATUS
  cout << __FILE__ << " : Starting solver..."<<endl;
#endif

  beSRS::run(NThreads, F, params.begin(), params.end(), outFilePrefix, string("#real(dBM)\timag(dBM)\t") + SCE_parameters::getFormat() + "\n", precisionGoal, !onlySaveGoal);


  cout << endl<< __FILE__ << " : SCE SOLVER FINISHED"<<endl<<endl;

  return 0;
}
