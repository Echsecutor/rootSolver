/**
 * @file derivativeVerification.cpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0.2014-06-04
 *
 *
 * @section DESCRIPTION
 *
 * This is a small test if you got your derivatives for the single
 * root solver right. ;)
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




// #include <boost/multiprecision/cpp_dec_float.hpp>
// using boost::multiprecision::cpp_dec_float_50;
// typedef cpp_dec_float_50 real_type;

#include <mpreal.h>
using namespace mpfr;
typedef mpreal real_type;


#include "preProDebugFlags.h"

#include <iostream>
#include <string>
#include <limits>
#include "complex.hpp"
#include <random>
#include <chrono>

using namespace std;

using namespace root_solver;

static root_solver::complex<real_type> I(0.0,1.0);
#define COMPLEX_I


//self consistency equations:
//(Include your own functions implementation instead)
#include "SCE.hpp"


int main(int args, char *arg[]){

  unsigned int precision=256;
  mpreal::set_default_prec(precision);
  root_solver::complex<real_type>::PI = const_pi(precision, mpreal::get_default_rnd());
  
  cout.precision(3);
  cout << std::scientific;

  default_random_engine gen(chrono::system_clock::now().time_since_epoch().count());
  normal_distribution<double> nd;

  real_type omega=1e-1;
  real_type dim=1;
  real_type epsilon = 1e-9;
  real_type bsbm = 1.0;

  real_type bmTarget=1e-2;///< will extrapolate from b=0 or from file towards bNuGoal
  real_type Db=1e-8;///< initial b step size, dynamically adapted
  real_type maxDb=1e-3;///< sets a minimal resolution in b
  real_type minDb=1e-10;///< avoid getting stuck


  list <SCE_parameters<real_type> > params;

  params.push_back(SCE_parameters<real_type>(dim,epsilon + I * omega));

  SCE<real_type> F (bsbm, bmTarget, maxDb, minDb, Db);

  F.setParameters(params.begin());

  F.set_extra_parameter(2.0);

  //start somewhere
  real_type one("1.0");
  root_solver::complex<real_type> x1,direction,x2, f1, f2, DfDz, J;
  real_type dist;

  cout << "numeric_limits<real_type>::epsilon() = " << numeric_limits<real_type>::epsilon() <<endl <<endl;

  root_solver::complex<real_type> eps(1.0),check=one+eps;
  while( abs(abs((check - one) / eps) - one) < 1e-2){
    eps *=1e-1;
    check=one+eps;
  }

  cout << "complex 1 ~ epsilon > " << eps <<endl <<endl;

  eps = (real_type)nd(gen) + I * (real_type)nd(gen);
  eps /= abs(eps);
  check=one+I+eps;

  while( abs(abs((check - one-I) / eps) - one) < 1e-2){
    eps *=1e-1;
    check=one+I+eps;
    //  cout <<eps<< "(" << check - one << ")"<<endl;
  }

  cout << "complex epsilon > " << eps <<endl <<endl;

  bool good = true;

  real_type maxErr=0.0;
  real_type maxFinalErr=0.0;

  for(int i=0;i<5 && good;i++){
    x1 = root_solver::complex<real_type>(nd(gen), nd(gen));
    F.changePoint(x1);
    f1 = F.calcF();
    J = F.calcJ();

    direction = root_solver::complex<real_type>(nd(gen),nd(gen));
    direction /= root_solver::complex<real_type>(abs(direction));
    cout << "checking point " << i +1 << ". z = " << x1 << " f = " << f1 << " dz ~ " << direction << " (|dz|="<<abs(direction)<<")"<<endl;
    cout << "derivative = " << J<<endl;
    dist = max(1e15 * numeric_limits<real_type>::epsilon(),1e-10);
    //    dist = "1.0";
    while( abs(dist) > max(1.0e5 * numeric_limits<real_type>::epsilon(),1e-20)){
      x2 = x1 + direction * dist;

      if(abs(abs((x2 - x1) / dist/direction) - one) > 1e-2){
	cout << "\n\nprecision underflow at dist = " << dist << ", direction*dist = " << direction*dist << ", x2-x1 = " << x2-x1 <<endl;
	cout << "abs(x2 - x1) = " << abs(x2 - x1) << ", rel err = ";
	  cout << abs(abs((x2 - x1) / dist/direction) - one)<<endl;
	return 1;
      }

      F.changePoint(x2);
      f2 = F.calcF();
      //    DfDz=(f2 - f1)/(direction * dist);//much worse results...
      DfDz=(f2 - f1) / (x2 - x1);
      cout << "|Dz| = " << dist << " => Df/Dz = " << (f2-f1) << " / " << (x2-x1) << " = " << DfDz << " \t(off by " << 100.0 * abs(DfDz - J)/abs(J) << " %)"<<endl;
      //      dist *= 1e-1;
      dist*=1e-1;
      if(abs(DfDz-J)/abs(J)>maxErr){
        maxErr=abs(DfDz-J)/abs(J);
      }
    }
      if(abs(DfDz-J)/abs(J)>maxFinalErr){
        maxFinalErr=abs(DfDz-J)/abs(J);
      }



    if(abs(DfDz-J)/abs(J) < 1e-2){
      cout <<"\n\nseems legitimate."<<endl;
    }else{
      cout << "does not look good... \n"<<endl;
      good =false;
    }
  }

  cout <<endl;

  if(good){
    cout << "\n[ok]\t\tmaximal relative error was " << maxErr <<endl;

    cout << "\t\tmaximal error at smallest step was " << maxFinalErr << "\n\n"<<endl;
    return 0;
  }

  cout << "\n[fail]\n\n"<<endl;
  return 1;
}

