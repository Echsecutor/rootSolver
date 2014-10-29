/**
 * @file testEigen.cpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0.2014-08-04
 *
 * @section DESCRIPTION
 *
 * Testing whether eigen works as expected
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

#include <eigen3/Eigen/Dense>
using namespace Eigen;

using namespace std;



template<typename realT>
bool doStuff(){

  srand(4277);

  Matrix<realT, 3, 3> *p =new Matrix<realT, 3, 3>;

  Matrix<realT, 3, 3> a,b;
  Matrix<realT, 3, 1> v,t;

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      a(i,j) = ((realT)rand()/(realT)RAND_MAX);
      b(i,j) = ((realT)rand()/(realT)RAND_MAX);
    }
    v(i)= ((realT)rand()/(realT)RAND_MAX);
  }

  *p = a * b;

  cout <<"(" << a << ")\n*\n(" << b << ")\n=\n(" << *p<<")\n"<<endl;

  cout <<"(" << a << ")\n*\n(" << v << ")\n=\n(" << a*v <<")\n"<<endl;

  delete p;
  p=&a;

  cout <<"(" << *p << ")\n+\n(" << b << ")\n=\n(" << *p + b <<")\n"<<endl;

  cout << 5.0 << " * (" << a << ")\n=(" << 5.0 * a << ")\n" << endl;

  ColPivHouseholderQR<Matrix<realT, 3,3>> QR(a);

  t=QR.solve(v);

  cout << "(" << a << ")^{-1}\n(" << v << ")\n = (" << t << ")\n"<<endl;

  if(! v.isApprox(a * t)){
    cerr << a * t << "\n != " << v <<"\nsolver failed!"<<endl;
    return false;
  }

  cout << "[OK]"<<endl;

  return true;
}

////////////////////////////////////////////////////////////////////////////////
/// main ;)
///
int main(int args, char *arg[]){

cout << "I am using Eigen " << EIGEN_WORLD_VERSION << "." << EIGEN_MAJOR_VERSION << "." << EIGEN_MINOR_VERSION <<endl;

#ifdef __GNUC__
cout << "Compiled with GNUC " << __GNUC__ << "." << __GNUC_MINOR__ <<endl;
#endif

  cout << "\ndouble precision:"<<endl;
  if(!doStuff<double>())
    return 1;

  cout << "\n\ndone."<<endl;


  return 0;
}
