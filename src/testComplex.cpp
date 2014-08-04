/**
 * @file testComplex.cpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0.2014-07-31
 *
 * @section DESCRIPTION
 *
 * Testing complex.hpp
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

#include "complexMprealWrapper.hpp"
using namespace mpfr;

using namespace std;

#include "complex.hpp"

using namespace root_solver;

template<typename realT>
void doStuff(){

  srand(42); 
  root_solver::initScalarType<realT>::ini();

  root_solver::complex<realT> a,b(1.0),c(42),d(0.0,1.0);

  root_solver::complex<realT>* pa = &a;
  root_solver::complex<realT>* p = new root_solver::complex<realT>(2);

  *pa = b * c;

  cout << b << "*" << c << "=" << a<<endl;

  cout << c << "/" << d << "=" << c/d <<endl;

  cout << *p << "+" << *pa << "=" << *p + *pa <<endl;

  cout << c << "-" << d << "=" << c - d <<endl;


  realT x= ((realT)rand()/(realT)RAND_MAX);
  realT y= ((realT)rand()/(realT)RAND_MAX);
  a= x + d * y;

  cout << "real + imag ( a = " << a << ") = " << real(a) << " + " << imag(a) <<endl;

  cout << "conj a = " << conj(a)<<endl;

  cout << "abs2 a = " << abs2(a)<<endl;

  cout << "abs a = " << abs(a)<<endl;

  cout << "exp a = " << exp(a) <<endl;

  cout << "log a = " << log(a) <<endl;

  cout << "arg a = " << root_solver::arg(a)<<endl;

  cout << "sqrt a = " << sqrt(a)<<endl;

  cout << "atan a = " << atan(a)<<endl;


  cout << "test deleting... " <<endl;

  delete p;
}

////////////////////////////////////////////////////////////////////////////////
/// main ;)
///
int main(int args, char *arg[]){


  cout << "\ndouble precision:"<<endl;
  doStuff<double>();
  cout <<"\n zomfg precision:"<<endl;
  doStuff<mpreal>();


  cout << "\n\ndone."<<endl;


  return 0;
}
