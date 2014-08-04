/**
 * @file SCE.hpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0.2014-06-11
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

#ifndef USEEXACTINT
#define USEEXACTINT 0
#endif


#include <ostream>
#include <stdexcept>
#include <list>
#include <limits>
using namespace std;


#include "singleRootSolver.hpp"
#include "batchSolver.hpp"
#include "extrapolationSolver.hpp"

#include "complex.hpp"
using namespace root_solver;



//------------------------------------------------------------------------------------------
// Declarations:
//------------------------------------------------------------------------------------------

template<typename realT>
class SCE_parameters : public virtual ostream_insertable{
public:

  typedef root_solver::complex<realT> complex_type;

  realT dim;
  realT bs, bm, cutOff;
  root_solver::complex<realT> z,g;

  SCE_parameters():dim(0),bs(0),bm(0),cutOff(0),z(0.0),g(0.0){}
  SCE_parameters(realT dim_,root_solver::complex<realT> z_):dim(dim_),bs(0),bm(0),cutOff(0),z(z_),g(0.0){}

  // automatic member wise copy sufficient
  //  SCE_parameters(const SCE_parameters<realT>& p2):dim(p2.dim),bs(p2.bs),bm(p2.bm),cutOff(p2.cutOff),z(p2.z),g(p2.g){}

  virtual void insertMeInto(ostream &out) const {
    out << dim << "\t" << bs << "\t" << bm << "\t" << real(z) << "\t" << imag(z) << "\t" << real(g) << "\t" << imag(g);
    //cutOff is computed from dim, irrelevant
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

  // automatic assignment operator should be fine

};


//-----------------------------------------------------------------------------------
/*doesnt work with g++ 4.6.3 ... o0
  template<typename realT> using eSRS = extraSolver<singleRootSolver<complex<realT> > >;

  template<typename realT> using beSRS = batchSolver<extraSolver<singleRootSolver<complex<realT> > >, list<SCE_parameters<realT> >::iterator>;
*/

template<typename realT>
class SCE : public virtual batch_functions<extraSolver<singleRootSolver<realT>, realT >, typename list<SCE_parameters<realT> >::iterator > {
private:

  typedef batch_functions<extraSolver<singleRootSolver<root_solver::complex<realT> >, realT >, typename list<SCE_parameters<realT> >::iterator > parent;

  //init flags
  bool dBMSet;

  //parameters:
  SCE_parameters<realT> p;
  realT n;///< dim = 2 n +1 or dim =2n
  bool oddDim;

  realT bs_bm;///< fixed ratio \f$ \frac{b_s}{b_m} \f$
  realT target_bm, max_change, min_change, initial_change;

  //temporaries read at changePoint:
  root_solver::complex<realT> BS,dBS,BM,dBM,a,b,Int;


  void computeInt();

  root_solver::complex<realT> arctanOverX(root_solver::complex<realT> x);

public:

  typedef root_solver::complex<realT> complex_type;

  static const complex_type one;
  static const complex_type two;


  SCE(realT bs_over_bm, realT final_bm, realT max_d_bm, realT min_d_bm, realT ini_d_bm):dBMSet(false),bs_bm(bs_over_bm),target_bm(final_bm),max_change(max_d_bm),min_change(min_d_bm),initial_change(ini_d_bm){}

  // automatic member wise copy sufficient
  //  SCE(const SCE<realT>& other):dBMSet(false),p(other.p), bs_bm(other.bs_bm),target_bm(other.target_bm),max_change(other.max_change),min_change(other.min_change),initial_change(other.initial_change){}


  //Overrides from
  //batch:
  virtual void setParameters(const typename SCE<realT>::iterator& It);
  virtual const SCE_parameters<realT>& getParameters() const{return p;};
  virtual SCE<realT>* clone() const;


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

//stupidly the implicit double to complex<realT> cast doesn't work due
//to ambiguities
template<typename realT>
const typename SCE<realT>::complex_type SCE<realT>::one = complex_type(1.0);
template<typename realT>
const typename SCE<realT>::complex_type SCE<realT>::two = complex_type(2.0);



//------------------------------------------------------------------------------------------
// Implementations:
//------------------------------------------------------------------------------------------

template<typename realT>
SCE<realT>* SCE<realT>::clone() const{
  return new SCE(*this);
}

/// the cut off is \f$ \Omega^{2n+1} = (2n+1) \pi^{n+1}\frac{2n!}{n!} \f$ for \f$ d = 2n+1 \f$
template<typename realT>
void SCE<realT>::setParameters(const typename SCE<realT>::iterator& It){
  p = *It;

  n = floor(p.dim / 2.0);

  oddDim = p.dim - 2.0 *n > 0.0;

  assert (complex_type::PI != 0.0);

  realT fac = complex_type::PI;
  for(realT f = n + 1.0; f <= 2.0 * n; f++){
    fac *= f * complex_type::PI;
  }

  p.cutOff = pow( (2.0 * n + 1.0) * fac, 1.0 / (2.0 * n + 1.0) );

#if DEBUG>=SPAM
  cout << __FILE__<< " : cutOff = " << p.cutOff << " in dimension " << p.dim << " = 2 * " <<  n;
  if(oddDim){
    cout << " + 1.0";
  }
  cout  <<endl;
#endif
}

template<typename realT>
realT SCE<realT>::get_extra_parameter(){
  return p.bm;
}

template<typename realT>
void SCE<realT>::set_extra_parameter(realT b){
  p.bm=b;
  p.bs = bs_bm * b;
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
void SCE<realT>::changePoint(const typename SCE<realT>::value_type x){
  if(x != x)
    throw std::runtime_error(std::string(__FILE__) + std::string(" : Can not change Point to NaN!"));

  dBM = x;
  BM = one + p.bm * dBM;

  p.g = BM * (dBM + one) / p.z;

  dBS = - p.g * BM * (BM + p.bm) / (p.z + p.bs * p.g * BM * (BM + p.bm));
  BS = one + p.bs * dBS;


  //#ifdef DIRAC
  //  a = p.g * p.z;///<Dirac model
  //#else
  a = p.g * p.z + p.bs * dBS *(one - p.dim);///< localQ model
  //#endif


  b = - BS * dBS / two / p.dim;

  computeInt();

  dBMSet=true;
}

template<typename realT>
typename SCE<realT>::complex_type SCE<realT>::arctanOverX(complex_type x){
  if (abs(x) < 1e-6){
    return one - x * x / (two+ one);//at abs(x)=0.001 the relative error of this approximation is about 2e-13
  }
  return atan(x) / x;
}


/// The actual computation is done at changePoint
template<typename realT>
typename SCE<realT>::value_type SCE<realT>::calcF(){
  if(!dBMSet)
    changePoint(dBM);//make sure internals are up to date


  if(real(p.g)<-1e-10){
    return realT(1e10);
    //    return mpfr::const_infinity(); //inf quickly leads to nan, nan leads to errors, errors lead to suffering... ;)
  }


  //  return BM;
  //  return p.g;
  //  return dBS;
  // return a + two * p.dim * b;//C
  // return - two * b; //D
  //  return a;
  // return b;
  //   return Int;

  //enforce positive density of states. (There seems to be another solution for \f$ - \bar g \f$.)

  return BS * Int - one;
  //todo: idea: use log(BS*Int) instead???

}


template<typename realT>
void SCE<realT>::computeInt(){

  Int=0.0;

#if USEEXACTINT>0
  if(p.dim==1.0){
    complex_type D = - two * b;
    complex_type C = a + two * p.dim * b;
    complex_type rt= sqrt(C*C - D*D);

#if DEBUG >= SPAM
    cout << __FILE__ << " : Using exact 1d integral, C=" << C << ", D=" << D <<endl;
#endif


    Int = one / rt;
    if(abs(C + rt) < abs(D)){//correct branch of sqrt
      Int*=-one;
    }

    return;

  }
#endif


  // This implements the Debeye approximation of the integral in question. In d=1 we would have:
  // \f[
  //   Int  = \int_0^\Omega \frac{d k}{a+b k^2} = \frac{1}{a}\sqrt{\frac{a}{b}}\arctan\left(\sqrt{\frac{b}{a}}\Omega\right)
  // \f]
  // Using the same sqrt twice cancels sign ambiguities. See paper for details about higher d.
#if DEBUG >= DETAIL
  cout << __FILE__ << " : Using Debeye approximation for integrals in dimension " << p.dim <<endl;
#endif


  if(oddDim){
    complex_type rt= sqrt(b/a);

    Int = arctanOverX(rt * p.cutOff);

    complex_type pow(1.0);
    complex_type powAB(1.0);

    for(realT k=0.0;k<n;k++){
      Int -= pow / (two * k + one);
      pow *= - b / a * p.cutOff * p.cutOff;
      powAB *= - a / b;
    }

    Int *= p.cutOff / a * powAB;

  }else{
    complex_type X = - a / b /p.cutOff / p.cutOff;

    Int = log((X - one)/X);

    complex_type pow(1.0);
    complex_type powAB(1.0);

    for(realT k=1.0;k<=n;k++){
      pow /= X;
      Int += pow / k;

      powAB *= - a / b;
    }

    Int *= powAB / two / b;

  }


}



template<typename realT>
typename SCE<realT>::derivative_type SCE<realT>::calcJ(){

  if(real(p.g)<-1e-10){
    return realT(0.0);
  }


  complex_type Dg = (one + two * p.bm * dBM + p.bm) / p.z;

  complex_type DBM = p.bm;

  complex_type DdBS = - BS * BS / p.z *(BM * (BM + p.bm) * Dg + p.g * (two * BM + p.bm) * DBM);

  //#ifdef DIRAC
  //  complex_type Da = p.z * Dg ;///< Dirac model
  //#else
  complex_type Da = p.z * Dg + p.bs * DdBS *(one - p.dim);///< local Q model
  //#endif

  complex_type Db = - (one + two * p.bs * dBS) * DdBS / two / p.dim;

  complex_type DInt(0.0);

#if USEEXACTINT>0
  if(p.dim == 1.0){
    complex_type D = - two * b;
    complex_type C = a + two * p.dim * b;
    complex_type DD = - two * Db;
    complex_type DC = Da + two * p.dim * Db;
    //return DC;
    //    return DD;
    DInt = Int*Int*Int * (D * DD - C * DC);

  }else{
#endif


    if(oddDim){

#if DEBUG >= DETAIL
    cout << __FILE__ << " : Using Debeye approximation for derivatives of integrals in odd dimension " << p.dim <<endl;
#endif

      complex_type X = sqrt(b/a) * p.cutOff;
      complex_type DX = X / two * (Db / b - Da / a);

      DInt = one / (one + X*X) - arctanOverX(X) ;


      complex_type pow(1.0);
      complex_type powAB(1.0);
      for(realT k=0.0;k<n;k++){
        DInt +=  two * k / (two * k + one) * pow;

        powAB *= - a / b;
        pow *= X*X;
      }

      DInt *= p.cutOff / a * powAB * DX / X;

      DInt += Int * ((n - one) * Da / a - n * Db / b);

    }else{

#if DEBUG >= DETAIL
    cout << __FILE__ << " : Using Debeye approximation for derivatives of integrals in even dimension " << p.dim <<endl;
#endif


      complex_type X = - a / b / p.cutOff/p.cutOff;
      complex_type DX = X * (Da / a - Db / b);

      DInt = one / (X - one);

      complex_type pow(1.0);
      complex_type powAB(1.0);
      for(realT k=1.0;k<=n;k++){
        pow /= X;
        DInt -= pow ;

        powAB *= - a / b;
      }

      DInt *= powAB/two/b *DX/X;
      DInt += n * Int *(Da/a-Db/b) - Int * Db/b;

    }

#if USEEXACTINT>0
  }
#endif



  // return DBM;
  //  return Dg;
  //  return DdBS;
  //  return Da;
  // return Db;
  //    return DInt;

  return Int * DdBS * p.bs + BS * DInt;

}


/// The deterministic limit is the well known harmonic crystal
template<typename realT>
typename SCE<realT>::value_type SCE<realT>::guessStartPoint(){

  a = p.z*p.z;
  b = one / two / p.dim;

  computeInt();

  p.g = Int * p.z;
  dBM = p.g * p.z -one;

  return dBM;

}




#endif

