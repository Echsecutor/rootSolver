/**
 * @file testExtraData.cpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0.2014-05-21
 *
 * @section DESCRIPTION
 *
 * Test case for the extrapolationData class.
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

#define DEBUG 9

#include <iostream>
#include <vector>
#include <memory>

using namespace std;

#include <eigen3/Eigen/Dense>
using namespace Eigen;


#include "extrapolationData.hpp"
#include "polynomials.hpp"


////////////////////////////////////////////////////////////////////////////////
/// main ;)
///
int main(int args, char *arg[]){
  cout << "Extrapolation test case" <<endl <<endl;

  const int dataDim=5;
  typedef complex<double> numericalT;
  typedef Matrix<numericalT,dataDim,1> dataT;

  vector<polynomial> P;
  P.resize(dataDim);

  extrapolationData<Matrix<numericalT,dataDim,1>, numericalT, dataDim> extra;

  cout << "Testing whether the fit is accurate:"<<endl;

  cout << "\tpreparing data..."<<endl;

  for(double x=0;x<10;x++){
    extraDat<dataT> dat;
    dat.extra=x;
    dat.dat=shared_ptr<dataT>(new dataT);
    for(int i=0;i<dataDim;i++){
      P.at(i).changePoint(x);
      (*(dat.dat))(i) = (P.at(i).calcF());
    }
    extra.push_back(dat);
  }

  cout << "\textrapolating..."<<endl;

  double to=42;

  Matrix<numericalT,dataDim,1> ext = extra.extrapolate(to);
  Matrix<numericalT,dataDim,1> exact;

  cout << "\tcalculating exact result..."<<endl;

  for(int i=0;i<dataDim;i++){
    P.at(i).changePoint(to);
    exact(i) = P.at(i).calcF();
  }

  cout << "\nThe extrapolated value (from [0;10] to " << to <<") is \n" << ext<<endl;
  cout << "\nThe exact value is \n" << exact <<endl;
  cout << "\nThe extrapolation error is " << (exact - ext).norm() <<endl;

  cout << "\nDone.\n"<<endl;

  return 0;
}
