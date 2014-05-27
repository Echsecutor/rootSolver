/**
 * @file extrapolationData.hpp
 * @author Sebastian Schmittner <sebastian@schmittner.de>
 * @version 1.0-2014-05-21
 *
 * @section DESCRIPTION
 *
 * This file defines a container class to store multidimensional data
 * of the form (valueT f, double x). It behaves like (and internally
 * uses) a std::list, but additional implements extrapolation of
 * values. I.e. using the known (f_i, x_i) a polynomial function p(x)
 * is fitted to minimise \f$(\sum_i ||f_i - p(x_i)||^2 \f$.
 *
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


#ifndef EXTRAPOLATIONDATA_HPP
#define EXTRAPOLATIONDATA_HPP

#include "preProDebugFlags.h"

#include <iostream>
#include <stdexcept>
#include <string>
#include <list>
#include <vector>
#include <memory>

#include <random>
#include <chrono>
using namespace std;

#include <eigen3/Eigen/Dense>
using namespace Eigen;




/// this is just (valueT f, double x)
template <typename dataT>
class extraDat{
  //  enum { NeedsToAlign = (sizeof(dataT)%16)==0 };
public:
  //  EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)
  ///  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  shared_ptr<dataT> dat;
  double extra;
  extraDat(){}
  extraDat(const dataT& dat_, const double& extra_):dat(shared_ptr<dataT>(new dataT(dat_))),extra(extra_){}//ZOMFG make_shared, at least in gcc 4.6.3, creates something unaligned internally... wtf...
};




/// grml... this work around for partial member specialisation took me
/// ages... and its still pretty ugly....
template <typename dataT, typename numericalT, unsigned int dataDim>
class specialised{

public:
  numericalT& access(dataT& x, const unsigned int index){
    return x[index];
  }
  numericalT& access(shared_ptr<dataT> x, const unsigned int index){
    return (*x)[index];
  }
};

template <typename dataT, typename numericalT>
class specialised<dataT,numericalT,1>{

public:
  numericalT& access(dataT& x, const unsigned int index){
    return x;
  }
  numericalT& access(shared_ptr<dataT> x, const unsigned int index){
    return *x;
  }
};




/**
 *
 * The data type dataT X should be accessible via X[i] for \$f i \in
 * [0,\ldots,dataDim-1]\f$ where X[i] should be of numericalT which in
 * turn should be some numerical type such as float, double,
 * complex<double>, etc.
 *
 * A typical example would be
 * extrapolationData<Matrix<valueT,dataDim,1>, numericalT, dataDim>
 * but a pointer to an appropriate array of numericalT would also be
 * an acceptable dataT.
 *
 * Caution, if dataDim=1 then dataT = numericalT is expected, i.e. no
 * extra (de)referencing is performed!
 *
 */
template <typename dataT, typename numericalT, unsigned int dataDim, unsigned int degree=3>
class extrapolationData : protected list<extraDat<dataT> >{//Eigen alignement... I should check wether the speed gain is worth the effort...
protected:

  //data for the extrapolated polynomial

  vector<shared_ptr<dataT > > coeffs;//now allignment gets REALLY annoying....maybe I should swith to dynamical allocation giving up some eigen optimisation...?

  //keep track of whether the coefficients were already calculated
  bool extrapolated;

  bool calculateCoefficients();
  void forget();

  static specialised<dataT, numericalT,dataDim> accessor;


public:
  typedef extraDat<dataT> value_type;
  typedef numericalT numerical_type;
  typedef extrapolationData<dataT,numericalT,dataDim> own_type;


  /// produce a random dataT vector with iid gaussian components of mean 0 and std. devaition 1
  /// caution, the return is real valued! If you want imaginary iid just call again and multiply by I
  static dataT gaussBlob(default_random_engine*& gen =0, int seed=42);

  /// This function fits an order degree polynomial to the stored data
  dataT extrapolate(double to);


  //The following are all esentially just passing to the corresponding
  //list function but also forgetting the now outdated extrapolation

  extrapolationData():list<value_type>(), extrapolated(false){}
  extrapolationData(const extrapolationData&x):list<value_type>(x), extrapolated(false){}

  void push_front (const value_type& val){list<value_type>::push_front(val);forget();}
  void push_front (const dataT& dat, const double& extra){
    push_front(shared_ptr<extraDat<dataT> >(new extraDat<dataT>(dat,extra)));
    forget();
  }
  void push_back (const value_type& val){list<value_type>::push_back(val);forget();}
  void push_back (const dataT& dat, const double& extra){
#if DEBUG >=SPAM
    cout << __FILE__ << " : allocating new data" <<endl;
#endif
    extraDat<dataT> t(dat, extra);
#if DEBUG >=SPAM
    cout << __FILE__ << " : created" <<endl;
#endif
    push_back(t);
#if DEBUG >=SPAM
    cout << __FILE__ << " : forget..." <<endl;
#endif
    forget();
#if DEBUG >=SPAM
    cout << __FILE__ << " : alocating new data done" <<endl;
#endif

  }

  void pop_back(){list<value_type>::pop_back();forget();}
  void pop_front(){list<value_type>::pop_front();forget();}
  void clear(){list<value_type>::clear();forget();}

  void sort(){list<value_type>::sort();forget();}
  template <class Compare>
  void sort (Compare comp){list<value_type>::sort(comp);forget();}

  typename list<value_type>::size_type size(){return list<value_type>::size();}

  //if you need more of the usual publics of list, just add forwards
  //as above ;)

};


template <typename dataT, typename numericalT, unsigned int dataDim, unsigned int degree>
specialised<dataT, numericalT,dataDim> extrapolationData<dataT,numericalT,dataDim,degree>::accessor;

//----------------------------------------------------------------------

template <typename dataT, typename numericalT, unsigned int dataDim, unsigned int degree>
void extrapolationData<dataT,numericalT,dataDim,degree>::forget(){
#if DEBUG > SPAM
  cout <<  __FILE__ << " : Old extrapolation invalid."<<endl;
#endif
  extrapolated=false;
  coeffs.clear();
}



template <typename dataT, typename numericalT, unsigned int dataDim, unsigned int degree>
dataT extrapolationData<dataT,numericalT,dataDim,degree>::extrapolate(double to){

#if DEBUG > SPAM
  cout <<  __FILE__ << " : Extrapolated point requested."<<endl;
#endif

  if(!extrapolated){
    if(!calculateCoefficients()){
      //use constant extrapolation instead
      value_type last = this->back();
      this->clear();
      this->push_back(last);
      calculateCoefficients();
    }
  }

  if(coeffs.size() == 0){
    throw std::runtime_error(string(__FILE__) + string(" : Error extrapolating 0 data points."));
  }
    
  

#if DEBUG > SPAM
  cout <<  __FILE__ << " : " << coeffs.size() << " coefficients known for degree " << degree << " polynomial" <<endl;
#endif

  dataT re;

  typename vector<shared_ptr<dataT > >::iterator it = coeffs.begin();

  // #if DEBUG > SPAM
  //   cout <<  __FILE__ << " : First coefficient is " << *(*it) <<endl;
  // #endif

  re = *(*it);

  double pow = to;
  while(++it != coeffs.end()){
    re += *(*it) * pow;
    pow*=to;
    // #if DEBUG > SPAM
    //   cout <<  __FILE__ << " : next coefficient is " << *(*it) <<endl;
    // #endif
  }

#if DEBUG > SPAM
  cout <<  __FILE__ << " : result: \n" << re <<endl;
#endif

  return re;

}


//-------------------------------------------------------------------------------------------


template <typename dataT, typename numericalT, unsigned int dataDim, unsigned int degree>
dataT extrapolationData<dataT,numericalT,dataDim,degree>::gaussBlob(default_random_engine*& gen,int seed){

  if(gen ==0){
#if DEBUG>=WARN
    cout << __FILE__ << " : Caution! Creating new random number generator. This should happen at most once!"<<endl;
#endif
    if (seed==0){
      seed = chrono::system_clock::now().time_since_epoch().count();
    }
    gen = new default_random_engine(seed);
  }

  normal_distribution<double> nd;
  
  dataT re;
  for(unsigned int k=0;k<dataDim;k++){
    accessor.access(re,k) = nd(*gen);
  }
  return re;
}


//-------------------------------------------------------------------------------------------

template <typename dataT, typename numericalT, unsigned int dataDim, unsigned int degree>
bool extrapolationData<dataT,numericalT,dataDim,degree>::calculateCoefficients(){

  extrapolated=false;
#if DEBUG >= SPAM
  cout <<  __FILE__ << " : (re) calculating coefficients of the extrapolated polynomial from " <<this->size() << " data points."<<endl;
#endif

  coeffs.clear();

  if (this->size() ==0){
#if DEBUG >= ERR
    cerr << __FILE__ << " : no data points to fit to."<<endl;
#endif
    return false;
  }


  for(unsigned int i=0;i<=degree;i++){
    coeffs.push_back(shared_ptr<dataT >(new dataT));
  }


  if (this->size() < degree+1){
#if DEBUG >= WARN
    cout << __FILE__ << " : " << this->size()<< " data points are insufficient to fit an order " << degree << " polynomial. Fitting a constant instead."<<endl;
#endif

    for (unsigned int k=0; k<dataDim; k++){
      accessor.access((*(*(coeffs.begin()))),k) = accessor.access((*(this->rbegin())).dat,k);
      for(typename vector<shared_ptr<dataT > >::iterator it = ++coeffs.begin();it != coeffs.end();++it){
        accessor.access((*(*it)),k) = 0.0;
      }
    }

    extrapolated=true;

    return extrapolated;
  }


  Matrix<numericalT, degree+1, degree+1> powers;

  for(unsigned int m=0;m<=degree;m++){
    for(unsigned int n=0;n<=degree;n++){
      powers(m,n)=0;
      for(typename list<value_type>::iterator it = this->begin(); it != this->end();it++){
        powers(m,n) += pow((*it).extra, m+n);
      }
    }
  }

#if DEBUG > SPAM
  cout <<  __FILE__ << " : Matrix of powers generated."<<endl;
#endif

  ColPivHouseholderQR<Matrix<numericalT, degree+1,degree+1>> QR(powers);

#if DEBUG > SPAM
  cout <<  __FILE__ << " : Computed QR decomposition."<<endl;
#endif

  for(unsigned int k=0; k<dataDim;k++){

    Matrix<numericalT, degree+1,1> target;
    for(unsigned int n=0;n<=degree;n++){
      target(n)=0;
      for(typename list<value_type>::iterator it = this->begin(); it != this->end();it++){
        target(n) += accessor.access((*it).dat,k) * pow((*it).extra, n);
      }
    }

#if DEBUG > SPAM
    cout <<  __FILE__ << " : " << k << "th target vector computed."<<endl;
#endif

    //actual computation:
    Matrix<numericalT, degree+1,1> coefficient = QR.solve(target);

#if DEBUG > SPAM
    cout <<  __FILE__ << " : " << k << "th coefficient extrapolated."<<endl;
#endif


    if(! target.isApprox(powers * coefficient)){
      cerr << __FILE__ << " : " << "Extrapolation failed."<<endl;
      return false;
    }

    for(unsigned int n=0;n<=degree;n++){
      accessor.access((*(coeffs.at(n))),k)=coefficient(n);
    }

#if DEBUG > SPAM
    cout <<  __FILE__ << " : " << k << "th coefficient saved."<<endl;
#endif

  }


  extrapolated=true;

#if DEBUG > SPAM
  cout <<  __FILE__ << " : Extrapolation finished."<<endl;
#endif


  return extrapolated;
}






#endif
