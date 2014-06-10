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

/// change the following to adapt the internal precision
//typedef long double real_type;

#include <mpreal.h>
using namespace mpfr;
typedef mpreal real_type;



#include "preProDebugFlags.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>
#include <ctime>


using namespace std;


//#include "logger.hpp"

//root finding:
#include "complex.hpp"
#include "singleRootSolver.hpp"
#include "batchSolver.hpp"
#include "extrapolationSolver.hpp"

using namespace root_solver;

//self consistency equations:
#include "SCE.hpp"





void readFile(const char * confFile,
              real_type &bmTarget,
              real_type &initialDb,
              real_type &epsilon,
              real_type &precisionGoal,
              int &desiredNrSteps,
              real_type &maxDb,
	      real_type &fromOmega,
	      real_type &toOmega,
	      real_type &DOmega,
	      int &NThreads,
	      bool &onlySaveGoal,
	      real_type &minDb,
	      real_type &dim,
	      bool &logStep,
	      real_type &bsbm
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
        dim=(real_type)atoi(val);
      }else if(strcmp(para,"logStep")==0){
        logStep=(bool)atoi(val);
      }
      else{
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

  real_type dim=1;
  real_type epsilon = 1e-9;
  real_type bsbm = 1.0;

  real_type fromOmega = 1e-3;
  real_type toOmega = 2;
  real_type DOmega = 1e-1;

  bool logStep=false;///< toggle wether omega points are equally spaced on a linear or a log scale

  real_type bmTarget=1e-2;///< will extrapolate from b=0 or from file towards bNuGoal
  real_type Db=1e-8;///< initial b step size, dynamically adapted
  real_type maxDb=1e-3;///< sets a minimal resolution in b 
  real_type minDb=1e-10;///< avoid getting stuck

  real_type precisionGoal=1e-12;
  int desiredNrSteps=64;

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

  cout << "The parameter format is\n" << SCE_parameters<real_type>::getFormat() <<endl<<endl;

#if DEBUG>=WARN
  if (precisionGoal < 1e-12){
    cout << __FILE__ << " : CAUTION! Precision goal is very small. Function evaluation should not be expected to be precise to more than 1e-12!" << endl;
  }
#endif



#if DEBUG>=WARN
  if(precisionGoal < 1e3 * numeric_limits<real_type>::epsilon()){
    cout << __FILE__ << " : CAUTION! Precision goal " << precisionGoal << " is close to the lower machine precision error bound " << numeric_limits<real_type>::epsilon() <<endl;
  }
#endif

  if(precisionGoal < numeric_limits<real_type>::epsilon()){
    throw runtime_error(string(__FILE__) + string(" : Precision goal is lower than machine precision. Recompile this solver with higher internal precision!"));
  }


#if DEBUG >= STATUS
  cout << __FILE__ << " : Initialising parameter list..."<<endl;
#endif
  list <SCE_parameters<real_type> > params;
  real_type omega=fromOmega;
  while(omega < toOmega){
    params.push_back(SCE_parameters<real_type>(dim, SCE_parameters<real_type>::complex_type(epsilon, omega)));

    omega += DOmega;

    if(logStep){
      DOmega = omega * omega / (omega - DOmega) - omega;
    }

  }
  params.push_back(SCE_parameters<real_type>(dim, SCE_parameters<real_type>::complex_type(epsilon,toOmega)));

#if DEBUG >= STATUS
  cout << __FILE__ << " : I will work on " << params.size() << " parameter sets"<<endl;
  cout << __FILE__ << " : Preparing self consistency equation."<<endl;
#endif

  SCE<real_type> * F = new SCE<real_type>(bsbm, bmTarget, maxDb, minDb, Db);

#if DEBUG >= STATUS
  cout << __FILE__ << " : Starting solver..."<<endl;
#endif

  batchSolver<extraSolver<singleRootSolver<real_type >, real_type>, list<SCE_parameters<real_type> >::iterator>::run(NThreads, F, params.begin(), params.end(), outFilePrefix, string("#real(dBM)\timag(dBM)\t") + SCE_parameters<real_type>::getFormat() + string("\n"), precisionGoal, !onlySaveGoal);


  cout << endl<< __FILE__ << " : SCE SOLVER FINISHED"<<endl<<endl;

  return 0;
}
