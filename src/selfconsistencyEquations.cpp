/**
 * @file selfconsistencyEquations.cpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0.2014-06-11
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

//#define MPREALPRECISION 256
#define MPREALPRECISION 512

#include "complexMprealWrapper.hpp"
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
#include "singleRootSolver.hpp"
#include "batchSolver.hpp"
#include "extrapolationSolver.hpp"

using namespace root_solver;

//self consistency equations:
#include "SCE.hpp"



struct confInfo{
  const char * confFile;

  const char * outDatFile;
  const char * inDatFile;

  real_type dim;

  real_type epsilon;

  real_type as;
  real_type am;

  real_type fromOmega;
  real_type toOmega;
  real_type DOmega;

  bool logStep;

  real_type initialDb;
  real_type maxDb;
  real_type minDb;

  disorder<real_type> b;

  approach_mode approach;

  real_type precisionGoal;
  int desiredNrSteps;

  int NThreads;

  bool onlySaveGoal;

  confInfo():
    confFile(0),
    outDatFile(0),
    inDatFile(0),
    dim(3),
    epsilon(1e-14),
    as(1.0),//alpha_s
    am(1.0),//alpha_m
    fromOmega(1e-5),
    toOmega (1),
    DOmega (1e-5),
    logStep(true),///< toggle wether omega points are equally spaced on a linear or a log scale
    initialDb(1e-8),///< initial b step size, dynamically adapted
    maxDb(1e-3),///< sets a minimal resolution in b
    minDb(1e-20),///< avoid getting stuck
    b(),
    approach(MANHATTANMAD),
    precisionGoal(1e-18),
    desiredNrSteps(64),
    NThreads(4),
    onlySaveGoal(true){b.mass.mix=2.0;b.mass.add=2.0;b.spring.mix=2.0;b.spring.add=2.0;}

};


void readFile(confInfo& I){

  //  logger::write(string("Reading parameters from ") + string(confFile), STATUS, __FILE__, logger::removePath);
#if DEBUG>=STATUS
  cout << __FILE__ << " : Reading parameters from " <<I.confFile<<endl;
#endif

  ifstream conf(I.confFile);
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

      if(strcmp(para,"bsMAD")==0||strcmp(para,"b_spring_mix")==0){
        I.b.spring.mix=atof(val);
      }else if(strcmp(para,"bsADD")==0 || strcmp(para,"b_spring_add")==0){
        I.b.spring.add=atof(val);
      }else if(strcmp(para,"bmADD")==0 || strcmp(para,"b_mas_add")==0){
        I.b.mass.add=atof(val);
      }else if(strcmp(para,"bmMAD")==0 || strcmp(para,"b_mas_mix")==0){
        I.b.mass.mix=atof(val);
      }else if(strcmp(para,"initialDb")==0){
        I.initialDb=atof(val);
      }else if(strcmp(para,"epsilon")==0){
        I.epsilon=atof(val);
      }else if(strcmp(para,"precisionGoal")==0){
        I.precisionGoal=atof(val);
      }else if(strcmp(para,"desiredNrSteps")==0){
        I.desiredNrSteps=atoi(val);
      }else if(strcmp(para,"maxDb")==0){
        I.maxDb=atof(val);
      }else if(strcmp(para,"fromOmega")==0){
        I.fromOmega=atof(val);
      }else if(strcmp(para,"toOmega")==0){
        I.toOmega=atof(val);
      }else if(strcmp(para,"DOmega")==0){
        I.DOmega=atof(val);
      }else if(strcmp(para,"NThreads")==0){
        I.NThreads=atoi(val);
      }else if(strcmp(para,"onlySaveGoal")==0){
        I.onlySaveGoal=(bool)atoi(val);
      }else if(strcmp(para,"minDb")==0){
        I.minDb=atof(val);
      }else if(strcmp(para,"dim")==0){
        I.dim=(real_type)atoi(val);
      }else if(strcmp(para,"logStep")==0){
        I.logStep=(bool)atoi(val);
      }else if(strcmp(para,"as")==0){
        I.as=atof(val);
      }else if(strcmp(para,"am")==0){
        I.am=atof(val);
      }else if(strcmp(para,"approach")==0){
        if(atoi(val)==0){
          I.approach=STRAIGHT;
        }else if(atoi(val)==1){
          I.approach=MANHATTANMAD;
        }else{
          I.approach=MANHATTANM;
        }
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

  initScalarType<root_solver::complex<real_type> >::ini();

  cout.precision(5);
  cout << std::scientific;

  confInfo I;

  //------------------------------------------------------------------------------------------
  // reading parameters and conf file
  //------------------------------------------------------------------------------------------

  for (int i=1;i<args;i++){
    if(strcmp(arg[i],"-h")==0||strcmp(arg[1],"--help")==0){
      help();
    }else if(strcmp(arg[i],"-c")==0){
      I.confFile = arg[i+1];
      i++;
    }else if(strcmp(arg[i],"-r")==0){
      I.inDatFile = arg[i+1];
      i++;
    }else if(strcmp(arg[i],"-o")==0){
      I.outDatFile = arg[i+1];
      i++;
    }else if (arg[i][0]=='-'){
      cout << "Unrecognised option " << arg[i] <<endl;
    }else if (I.confFile==0){
      I.confFile = arg[i];
    }else{
      I.outDatFile = arg[i];
    }
  }

  if(I.confFile==0){
    help();
  }

  string outFilePrefix="SCE-solutions__";

  if(I.outDatFile!=0){
    outFilePrefix = I.outDatFile;
  }



  readFile(I);

  cout << __FILE__ << " : Start solving saddle point equations for"<<endl;
  cout << "dim = " << I.dim <<endl;
  cout << "b_m^+ = " << I.b.mass.add <<endl;
  cout << "b_m^x = " << I.b.mass.mix <<endl;
  cout << "b_s^+ = " << I.b.spring.add <<endl;
  cout << "b_s^x = " << I.b.spring.mix <<endl;
  cout << "approach = " << I.approach <<endl;
  cout << "alpha_s = " << I.as <<endl;
  cout << "alpha_m = " << I.am <<endl;
  cout << "epsilon = " << I.epsilon <<endl;
  cout << "fromOmega = " << I.fromOmega <<endl;
  cout << "toOmega = " << I.toOmega <<endl;
  cout << "DOmega = " << I.DOmega <<endl;
  cout << "logStep = " << I.logStep <<endl;
  cout << "initialDb = " << I.initialDb <<endl;
  cout << "precisionGoal = " << I.precisionGoal<<endl;
  cout << "desiredNrSteps = " << I.desiredNrSteps <<endl;
  cout << "maxDb = " << I.maxDb <<endl;
  cout << "NThreads = " << I.NThreads <<endl;
  cout << "onlySaveGoal = " << I.onlySaveGoal <<endl;
  cout << "\n\n" << endl;

  cout << "The parameter format is\n" << SCE_parameters<real_type>::getFormat() <<endl<<endl;



#if DEBUG>=WARN
  if(I.precisionGoal < 1e3 * numeric_limits<real_type>::epsilon()){
    cout << __FILE__ << " : CAUTION! Precision goal " << I.precisionGoal << " is close to the lower machine precision error bound " << numeric_limits<real_type>::epsilon() <<endl;
  }
#endif

  if(I.precisionGoal < numeric_limits<real_type>::epsilon()){
    throw runtime_error(string(__FILE__) + string(" : Precision goal is lower than machine precision. Recompile this solver with higher internal precision!"));
  }


#if DEBUG >= STATUS
  cout << __FILE__ << " : Initialising parameter list..."<<endl;
#endif
  list <SCE_parameters<real_type> > params;
  real_type omega = I.fromOmega;
  while(omega < I.toOmega){
    params.push_back(SCE_parameters<real_type>(I.dim, SCE_parameters<real_type>::complex_type(I.epsilon, omega)));
    params.back().as = I.as;
    params.back().am = I.am;

    omega += I.DOmega;

    if(I.logStep){
      I.DOmega = omega * omega / (omega - I.DOmega) - omega;
    }

  }
  params.push_back(SCE_parameters<real_type>(I.dim, SCE_parameters<real_type>::complex_type(I.epsilon,I.toOmega)));
  params.back().as = I.as;
  params.back().am = I.am;


#if DEBUG >= STATUS
  cout << __FILE__ << " : I will work on " << params.size() << " parameter sets"<<endl;
  cout << __FILE__ << " : Preparing self consistency equation."<<endl;
#endif


  // SCE(const realT& target_bsMAD_,const realT& target_bsADD_, const realT& target_bmMAD_, const realT& target_bmADD_,const realT& max_d_b,  const realT& min_d_b, const realT& ini_d_b, const approach_mode& app):
   
  SCE<real_type> * F = new SCE<real_type>(I.b,I.maxDb,I.minDb,I.initialDb, I.approach);


#if DEBUG >= STATUS
  cout << __FILE__ << " : Starting solver..."<<endl;
#endif

  batchSolver<extraSolver<singleRootSolver<real_type >, real_type>, list<SCE_parameters<real_type> >::iterator>::run(I.NThreads, F, params.begin(), params.end(), outFilePrefix, string("#real(dBM)\timag(dBM)\t") + SCE_parameters<real_type>::getFormat() + string("\n"), I.precisionGoal, !I.onlySaveGoal);


  cout << endl<< __FILE__ << " : SCE SOLVER FINISHED"<<endl<<endl;

  return 0;
}
