/**
 * @file SCE.hpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0.2014-05-26
 *
 * @section DESCRIPTION
 *
 * Selfconsistency equations encountered in recent research on
 * disordered bosons at the Institute for Theoretical Physics,
 * University of Cologne, in the group of Prof. M. R. Zirnbauer.
 *
 * We need to solve one complex equation
 * \f[
 * 1 = \int_{[0;2 \pi]^d} \frac{d^d k}{(2\pi)^d} \frac{B_S}{-(\delta B_S) B_S \Delta_k + g z + b_s \delta B_S(1-d)}
 * \f]
 * with
 * \f[
 * \delta B_S = \frac{B_S - 1}{b_s} = - \frac{g B_M (B_M + b_m)}{z + b_s g B_M(B_M+b_m)}
 * \f]
 * \f[
 * g = \frac{B_M (\delta B_M +1)}{z}
 * \f]
 * and \f$ B_M = 1 + b_m \delta B_M \f$ where \f$ \delta B_M \f$ is the solver parameter.
 * This parametrisation has the advantage that the $b_s,b_m \to 0$ limit is unproblematic, also numerically.
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




#ifndef SCE_HPP
#define SCE_HPP

#include <ostream>

#include <complex>
#include<cmath>//M_PI
#include <list>
#include <limits>
using namespace std;


#ifndef COMPLEX_I
#define COMPLEX_I
static complex<double> I(0.0,1.0);
#endif


#include "singleRootSolver.hpp"
#include "batchSolver.hpp"
#include "extrapolationSolver.hpp"

using namespace root_solver;



//------------------------------------------------------------------------------------------
// Declarations:
//------------------------------------------------------------------------------------------

template<typename realT>
class SCE_parameters{
public:

  realT dim;
  realT bs, bm, cutOff;
  complex<realT> z,g;

  SCE_parameters():dim(0),bs(0),bm(0),cutOff(0),z(0.0),g(0.0){}
  SCE_parameters(realT dim_,complex<realT> z_):dim(dim_),bs(0),bm(0),cutOff(0),z(z_),g(0.0){}

  SCE_parameters(const SCE_parameters<realT>& p2):dim(p2.dim),bs(p2.bs),bm(p2.bm),cutOff(p2.cutOff),z(p2.z),g(p2.g){}


  friend ostream& operator<< (ostream &out, SCE_parameters<realT> &p){
    out << p.dim << "\t" << p.bs << "\t" << p.bm << "\t" << real(p.z) << "\t" << imag(p.z) << "\t" << real(p.g) << "\t" << imag(p.g);//cutOff is computed from dim, irrelevant
    return out;
  }

  static string getFormat(){
    return string("dim\tbs\tbm\treal(z)\timag(z)\treal(g)\timag(g)");
  }

  bool operator==(const SCE_parameters<realT>& rhs) {
    bool re = true;
    re = re && this->dim == rhs.dim;
    re = re && this->bs == rhs.bs;
    re = re && this->bm == rhs.bm;
    re = re && this->cutOff == rhs.cutOff;
    re = re && this->z == rhs.z;
    re = re && this->g == rhs.g;
    return re;
  }

};


//-----------------------------------------------------------------------------------
/*doesnt work with g++ 4.6.3 ... o0
template<typename realT> using eSRS = extraSolver<singleRootSolver<complex<realT> > >;

template<typename realT> using beSRS = batchSolver<extraSolver<singleRootSolver<complex<realT> > >, list<SCE_parameters<realT> >::iterator>;
*/

template<typename realT>
class SCE : public virtual batch_functions<extraSolver<singleRootSolver<complex<realT> > >, typename list<SCE_parameters<realT> >::iterator > {
private:

  typedef batch_functions<extraSolver<singleRootSolver<complex<realT> > >, typename list<SCE_parameters<realT> >::iterator > parent;

  //init flags
  bool dBMSet;

  //parameters:
  SCE_parameters<realT>* p;
  realT n;///< dim = 2 n +1

  realT bs_bm;///< fixed ratio \f$ \frac{b_s}{b_m} \f$
  realT target_bm, max_change, min_change, initial_change;

  //temporaries read at changePoint:
  complex<realT> BS,dBS,BM,dBM,a,b,Int;


  void computeInt();

  complex<realT> arctanOverX(complex<realT> x);

public:



  SCE(realT bs_over_bm, realT final_bm, realT max_d_bm, realT min_d_bm, realT ini_d_bm):dBMSet(false),p(0), bs_bm(bs_over_bm),target_bm(final_bm),max_change(max_d_bm),min_change(min_d_bm),initial_change(ini_d_bm){}


  //Overrides from
  //batch:
  virtual void setParameters(typename SCE<realT>::iterator It);
  //  using parent::setParameters;


  //extra:
  virtual realT get_extra_parameter();
  virtual void set_extra_parameter(realT p);

  virtual realT get_initial();
  virtual realT get_final();

  virtual realT get_max_change();

  virtual realT get_min_change();
  virtual realT get_initial_change();


  //solver:
  virtual typename SCE<realT>::value_type calcF();
  virtual typename SCE<realT>::derivative_type calcJ();
  virtual void changePoint(const typename SCE<realT>::value_type x);
  virtual typename SCE<realT>::value_type guessStartPoint();

};




//------------------------------------------------------------------------------------------
// Implementations:
//------------------------------------------------------------------------------------------

/// the cut off is \f$ \Omega^{2n+1} = (2n+1) \pi^{n+1}\frac{2n!}{n!} \f$ for \f$ d = 2n+1 \f$
template<typename realT>
void SCE<realT>::setParameters(typename SCE<realT>::iterator It){
  p = &(*It);

  n = floor((p->dim -1.0)/2.0);

  //only odd d implemented yet!
  if(2.0 * n + 1.0 != p->dim){
    throw runtime_error(string(__FILE__) + string(" : Debeye approximation only implemented for odd dimensions."));
  }

  realT fac = M_PI;
  for(realT f = n + 1.0; f <= 2.0 * n; f++){
    fac *= f * M_PI;
  }

  p->cutOff = pow( (2.0 * n + 1.0) * fac, 1.0 / (2.0 * n + 1.0) );

#if DEBUG>=SPAM
  cout << __FILE__<<" : cutOff = " << p->cutOff << " in dimension " << p->dim<<endl;
#endif
}

template<typename realT>
realT SCE<realT>::get_extra_parameter(){
  return p->bm;
}

template<typename realT>
void SCE<realT>::set_extra_parameter(realT b){
  p->bm=b;
  p->bs = bs_bm * b;
  dBMSet=false;
}

template<typename realT>
realT SCE<realT>::get_initial(){
  return 0.0;
}

template<typename realT>
realT SCE<realT>::get_final(){
  return target_bm;
}

template<typename realT>
realT SCE<realT>::get_max_change(){
  return max_change;
}

template<typename realT>
realT SCE<realT>::get_min_change(){
  return min_change;
}

template<typename realT>
realT SCE<realT>::get_initial_change(){
  return initial_change;
}


template<typename realT>
void SCE<realT>::computeInt(){

  Int=0.0;

  if(p->dim==1.0){
    complex<realT> D = - 2.0 * b;
    complex<realT> C = a + 2.0 * p->dim * b;
    complex<realT> rt= sqrt(C*C - D*D);

#if DEBUG >= SPAM
    cout << __FILE__ << " : Using exact 1d integral, C=" << C << ", D=" << D <<endl;
#endif


    Int = 1.0 / rt;
    if(abs(C + rt) < abs(D)){//correct branch of sqrt
      Int*=-1.0;
    }

  }else{

    // This implements the Debeye approximation of the integral in question. In d=1 we would have:
    // \f[
    //   Int  = \int_0^\Omega \frac{d k}{a+b k^2} = \frac{1}{a}\sqrt{\frac{a}{b}}\arctan\left(\sqrt{\frac{b}{a}}\Omega\right)
    // \f]
    // Using the same sqrt twice cancels sign ambiguities. See paper for details about higher d.
#if DEBUG >= DETAIL
    cout << __FILE__ << " : Using Debeye approximation for integrals in dimension " << p->dim <<endl;
#endif


    complex<realT> rt= sqrt(b/a);
    
    Int += arctanOverX(rt * p->cutOff);


    for(int k=0;k<n;k++){
      Int -= pow(- b / a * p->cutOff * p->cutOff , k) / (2.0 * (realT) k + 1.0);
    }

    Int *= p->cutOff/a * pow(- a / b, (int)n);

  }

}


template<typename realT>
void SCE<realT>::changePoint(const typename SCE<realT>::value_type x){
  dBM = x;
  BM = 1.0 + p->bm * dBM;

  p->g = BM * (dBM + 1.0) / p->z;

  dBS = - p->g * BM * (BM + p->bm) / (p->z + p->bs * p->g * BM * (BM + p->bm));
  BS = 1.0 + p->bs * dBS;

  a = p->g * p->z + p->bs * dBS *(1.0 - p->dim);
  b = - BS * dBS / 2.0 / p->dim;

  computeInt();

  dBMSet=true;
}

template<typename realT>
complex<realT> SCE<realT>::arctanOverX(complex<realT> x){
  if (abs(x) < 0.001){
    return 1.0 - x * x / 3.0;//at abs(x)=0.001 the relative error of this approximation is about 2e-13
  }
  return atan(x) / x;

}


/// The actual computation is done at changePoint
template<typename realT>
typename SCE<realT>::value_type SCE<realT>::calcF(){
  if(!dBMSet)
    changePoint(dBM);//make sure internals are up to date

  return BS * Int - 1.0;
  //todo: idea: use log(BS*Int) instead???
}


template<typename realT>
typename SCE<realT>::derivative_type SCE<realT>::calcJ(){

  complex<realT> Dg = (1.0 + 2.0 * p->bm * dBM + p->bm) / p->z;

  complex<realT> DBM = p->bm;
  complex<realT> DdBS= BS * dBS * Dg / p->g  + dBS*(1.0 - p->bs * dBS / BM / (BM + p->bm)) * (2.0 * BM + p->bm) / BM / (BM + p->bm) * DBM;

  complex<realT> Da = p->z * Dg + p->bs * DdBS *(1.0 - p->dim);
  complex<realT> Db = (- p->bs * dBS * DdBS - BS * DdBS) / 2.0 / p->dim;

  complex<realT> DInt=0.0;

  if(p->dim == 1.0){
    complex<realT> D = - 2.0 * b;
    complex<realT> C = a + 2.0 * p->dim * b;
    complex<realT> DD = - 2.0 * Db;
    complex<realT> DC = Da + 2.0 * p->dim * Db;


    DInt = Int*Int*Int * (D * DD - C * DC);

  }else{

    complex<realT> t = p->cutOff / a / (1.0 + b / a * p->cutOff*p->cutOff);
    DInt += 0.5 * (t - Int) * (Db / b   - Da / a);


    for(int k=1;k<n;k++){
      DInt -= pow( (- b / a  * p->cutOff * p->cutOff ), k ) / (2.0 * (realT)k + 1.0) * (realT)k * ( Db / b - Da / a) ;
    }

 
    DInt *= p->cutOff/a * pow(- a / b, (int)n);

    DInt += Int *((n - 1.0) * Da / a - n * Db / b);

  }

  return Int * DdBS * p->bs + BS * DInt;

}


/// The deterministic limit is the well known harmonic crystal
template<typename realT>
typename SCE<realT>::value_type SCE<realT>::guessStartPoint(){

  a = p->z*p->z;
  b = 1.0 / 2.0 / p->dim;

  computeInt();

  p->g = Int * p->z;
  dBM = p->g * p->z -1.0;

  return dBM;

}




#endif

