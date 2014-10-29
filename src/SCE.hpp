/**
 * @file SCE.hpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0.2014-09-24`
 *
 * @section DESCRIPTION
 *
 * Selfconsistency equations encountered in recent research on
 * disordered bosons at the Institute for Theoretical Physics,
 * University of Cologne, in the group of Prof. M. R. Zirnbauer.
 *
 * We need to solve one complex equation
 * (see phd thesis)
 * and \f$ B_M = 1 + b_m \delta B_M \f$ where \f$ \delta B_M \f$ is the solver parameter.
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



namespace root_solver{



  //------------------------------------------------------------------------------------------
  // Declarations:
  //------------------------------------------------------------------------------------------

  template<typename realT>
  struct disorder_imp{
    realT add;
    realT mix;
    disorder_imp(){add=0.0;mix=0.0;}
  };

  template<typename realT>
  struct disorder{
    disorder_imp<realT> mass;
    disorder_imp<realT> spring;
    disorder():mass(),spring(){}
  };


  template<typename realT>
  class SCE_parameters : public virtual ostream_insertable{
  public:

    typedef root_solver::complex<realT> complex_type;

    realT dim;
    realT cutOff, cutOffD, as,am;
    disorder<realT> b;
    root_solver::complex<realT> z,g;

    SCE_parameters():dim(0.0),cutOff(0.0),as(1.0),am(1.0),b(),z(0.0),g(0.0){}
    SCE_parameters(realT dim_,root_solver::complex<realT> z_):dim(dim_),cutOff(0.0),as(1.0),am(1.0),b(),z(z_),g(0.0){}

    // automatic member wise copy sufficient
    //  SCE_parameters(const SCE_parameters<realT>& p2):dim(p2.dim),bs(p2.bs),bm(p2.bm),cutOff(p2.cutOff),z(p2.z),g(p2.g){}

    virtual void insertMeInto(ostream &out) const {
      out << dim << "\t" << b.spring.mix << "\t" << b.mass.mix << "\t" << real(z) << "\t" << imag(z) << "\t" << real(g) << "\t" << imag(g) << "\t"<<  b.spring.add << "\t" << b.mass.add << "\t" <<as <<"\t"<<am;
      //cutOff is computed from dim, irrelevant
    }

    static string getFormat(){
      return string("dim\tb_s^mixed\tb_m^mix\treal(z)\timag(z)\treal(g)\timag(g)\tb_s^add\tb_m^add\talpha_s\talpha_m");
    }

    bool operator==(const SCE_parameters<realT>& rhs) {
      bool re = true;
      re = re && this->dim == rhs.dim;
      re = re && this->b.spring.add == rhs.b.spring.add;
      re = re && this->b.mass.add == rhs.b.mass.add;
      re = re && this->b.spring.add == rhs.b.spring.add;
      re = re && this->b.mass.add == rhs.b.mass.add;
      re = re && this->cutOff == rhs.cutOff;
      re = re && this->cutOffD == rhs.cutOffD;
      re = re && this->z == rhs.z;
      re = re && this->g == rhs.g;
      re = re && this->as == rhs.as;
      re = re && this->am == rhs.am;
      return re;
    }

    // automatic assignment operator should be fine

  };


  //-----------------------------------------------------------------------------------
  /*doesnt work with g++ 4.6.3 ... o0
    template<typename realT> using eSRS = extraSolver<singleRootSolver<complex<realT> > >;

    template<typename realT> using beSRS = batchSolver<extraSolver<singleRootSolver<complex<realT> > >, list<SCE_parameters<realT> >::iterator>;
  */



  typedef root_solver::solver_state approach_mode; //reuse custom enum
  static const approach_mode STRAIGHT(0);
  static const approach_mode MANHATTANMAD(1);
  static const approach_mode MANHATTANM(2);


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

    disorder<realT> target_b;
    approach_mode approach;

    int nonZero;//<for MANHATTAN approach
    realT** valueOrder;//<for MANHATTAN approach
    realT** targetOrder;//<for MANHATTAN approach


    realT b_rel, max_change, min_change, initial_change;

    //temporaries read at changePoint:
    root_solver::complex<realT> BS,BM,dBM,ZM,a,b,Int, BMADD,BSADD;


    void computeInt();

    root_solver::complex<realT> arctanOverX(root_solver::complex<realT> x);

  public:

    //  interpolate_in interpolate;

    typedef root_solver::complex<realT> complex_type;

    static const complex_type one;
    static const complex_type two;

    SCE(const disorder<realT>& target_b_,
        const realT& max_d_b, //< relative scale, i.e. for 0=<b<=1
        const realT& min_d_b,
        const realT& ini_d_b,
        const approach_mode& app):
      dBMSet(false),
      target_b(target_b_),
      approach(app),
      valueOrder(new realT*[4]),
      targetOrder(new realT*[4]),
      b_rel(0.0),
      max_change(max_d_b),
      min_change(min_d_b),
      initial_change(ini_d_b){}

    SCE():
      dBMSet(false),
      target_b(),
      approach(STRAIGHT),
      valueOrder(new realT*[4]),
      targetOrder(new realT*[4]),
      b_rel(0.0),
      max_change(1.0),
      min_change(1.0),
      initial_change(1.0){}

    ~SCE(){
      delete[] valueOrder;
      delete[] targetOrder;
    }

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


  //todo: quick n dirty, should go into apropriate trait
  template<typename realT>
  int toDouble(realT x){
    return (double) x;
  }

  template<>
  int toDouble<mpfr::mpreal>(mpfr::mpreal x){
    return x.toDouble();
  }


  /// the cut off is \f$ \Omega^{2n+1} = (2n+1) \pi^{n+1}\frac{2n!}{n!} \f$ for \f$ d = 2n+1 \f$
  template<typename realT>
  void SCE<realT>::setParameters(const typename SCE<realT>::iterator& It){
    p = *It;

    int iDim = floor(toDouble<realT>(p.dim));

    if((realT)((double)iDim) != p.dim){
      throw std::runtime_error(string(__FILE__)+string(" : I cant cast the real data type to int properly... o0"));
    }

    n = (double) floor((iDim-1) / 2);
    oddDim = (iDim % 2 == 1);

    assert (complex_type::PI != 0.0);

    realT VolS;
    if(oddDim){
      VolS = 2.0;
      for(realT i=0.0; i < p.dim - 1.0; i+=2.0){
        VolS *= 2.0 * complex_type::PI / (i+1.0);
      }
    }else{
      VolS = 2.0 * complex_type::PI;
      for(realT i=1.0; i < p.dim - 1.0; i+=2.0){
        VolS *= 2.0 * complex_type::PI / (i+1.0);
      }
    }
    //now VolS = Vol(S^{d-1})
#if DEBUG>=SPAM
    cout << __FILE__<< " : Vol(S^{d-1}) = " <<  VolS <<endl;
#endif

    p.cutOffD = p.dim * pow(2.0*complex_type::PI, iDim) / VolS;

    p.cutOff = pow( p.cutOffD, real(one /p.dim));


#if DEBUG>=SPAM
    cout << __FILE__<< " : cutOff ^ d = " << p.cutOff << " ^ " << p.dim << " = " << p.cutOffD << " in dimension d = " << p.dim <<  " = 2 * " <<  n;
    if(!oddDim){
      cout << " + 2";
    }else{
      cout << " + 1";
    }
    cout  <<endl;
#endif



    nonZero=4;
    if(approach==MANHATTANMAD){
      valueOrder[0]=&(p.b.mass.mix);
      valueOrder[1]=&(p.b.spring.mix);
      valueOrder[2]=&(p.b.mass.add);
      valueOrder[3]=&(p.b.spring.add);

      targetOrder[0]=&(target_b.mass.mix);
      targetOrder[1]=&(target_b.spring.mix);
      targetOrder[2]=&(target_b.mass.add);
      targetOrder[3]=&(target_b.spring.add);

    }else{
      valueOrder[0]=&(p.b.mass.mix);
      valueOrder[1]=&(p.b.mass.add);
      valueOrder[2]=&(p.b.spring.mix);
      valueOrder[3]=&(p.b.spring.add);

      targetOrder[0]=&(target_b.mass.mix);
      targetOrder[1]=&(target_b.mass.add);
      targetOrder[2]=&(target_b.spring.mix);
      targetOrder[3]=&(target_b.spring.add);

    }

    for(int i=0;i<nonZero;i++){
      if(*(targetOrder[i])==0.0){
        for(int j=i;j<nonZero -1;j++){
          valueOrder[j] = valueOrder[j+1];
          targetOrder[j] = targetOrder[j+1];
        }
        i--;
        nonZero--;
      }
    }

#if DEBUG>=DETAIL
    cout << __FILE__ << " : approach is " << approach <<endl;
    cout << __FILE__ << " : non-zero parameters in target " << nonZero <<endl;
    cout << "order " << *(targetOrder[0]) << ", "<< *(targetOrder[1]) << ", "<< *(targetOrder[2]) << ", "<< *(targetOrder[3])<<endl;

#endif

  }//end set Parameters


  template<typename realT>
  realT SCE<realT>::get_extra_parameter(){
    return b_rel;
  }

  template<typename realT>
  void SCE<realT>::set_extra_parameter(realT b){
    b_rel=b;

    if(approach==STRAIGHT){
      p.b.spring.add= b_rel *  target_b.spring.add;
      p.b.mass.add= b_rel *  target_b.mass.add;
      p.b.spring.mix= b_rel *  target_b.spring.mix;
      p.b.mass.mix= b_rel *  target_b.mass.mix;
    }else{

      //MANHATTAN path
      realT partLength = 1.0/(double)nonZero;

      for(int i=0;i<nonZero;i++){
        realT iR((double)i);
        if( b_rel <= (iR+1.0) * partLength && b_rel >= iR * partLength){
          *(valueOrder[i]) = (b_rel / partLength - iR) * *(targetOrder[i]);
#if DEBUG >=SPAM
          cout << __FILE__ << " : setting parameter " << i << " to " << *(valueOrder[i]) << " (target = " << *(targetOrder[i]) << ")"  <<endl;
#endif
        }else if( b_rel >= (iR+1.0) * partLength){//todo: kA y needed...
          *(valueOrder[i])=  *(targetOrder[i]);
        }else{
          *(valueOrder[i])=0.0;
        }

      }//next
    }//end if manhattan

    dBMSet=false;
  }//end set extra parameter

  template<typename realT>
  realT SCE<realT>::get_initial(){
    return 0.0;
  }

  template<typename realT>
  realT SCE<realT>::get_final(){
    return 1.0;
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
    BM = one + p.b.mass.mix * dBM;

    if(p.b.mass.add>0.0){//piece of cake for branching prediction =>speed up + less numerical errors
      BMADD =  BM / (BM - p.b.mass.add * dBM / p.am);

      ZM = BM * (BM + p.b.mass.mix) + p.b.mass.add * BMADD;

      p.g = (BM * (one + dBM) + p.b.mass.add * BMADD *dBM/BM )/p.z;

    }else{
      BMADD=one;
      ZM = BM * (BM + p.b.mass.mix);
      p.g = BM * (one + dBM) / p.z;
    }

    BS = p.z / ( p.z + p.b.spring.mix * p.g * ZM);

    if(p.b.spring.add>0.0){//s.above
      BSADD = p.z / ( p.z + p.b.spring.add / p.as * p.g * ZM);
      a = p.z*p.z + (p.dim - one + p.g * p.z) * p.b.spring.mix * ZM + p.b.spring.add * p.dim * ZM * BSADD/BS;
    }else{
      BSADD=one;
      a = p.z*p.z + (p.dim -one + p.g * p.z) * p.b.spring.mix *ZM;
    }

    b = ZM*BS;

    computeInt();

#if DEBUG >= CALCULATION
    cout << __FILE__ << " : changed point to " << x <<endl;
    cout << __FILE__ << " : with parameters b.mass.mix = " << p.b.mass.mix << ", b.mass.add = " << p.b.mass.add << ", alpha_m = " << p.am<<endl;
    cout << __FILE__ << " : b.spring.mix = " << p.b.spring.mix << ", b.spring.add = " << p.b.spring.add << ", alpha_s = " << p.as<<endl;

    cout << __FILE__ << " : dBM = " << dBM <<endl;
    cout << __FILE__ << " : BM = " << BM <<endl;
    cout << __FILE__ << " : DBMADD = " << DBMADD <<endl;
    cout << __FILE__ << " : ZM = " << ZM <<endl;
    cout << __FILE__ << " : g = " << p.g <<endl;
    cout << __FILE__ << " : DBSADD = " << DBSADD <<endl;
    cout << __FILE__ << " : BS = " << BS <<endl;
    cout << __FILE__ << " : a = " << a << ", b=" << b <<endl;
    cout << __FILE__ << " : Int = " << Int <<endl;
    cout <<endl;
#endif


    dBMSet=true;
  }

  template<typename realT>
  typename SCE<realT>::complex_type SCE<realT>::arctanOverX(complex_type x){

    return atan(x)/x;

    complex_type re=0.0;
    if (abs(x) < 1e-6){
      re= one - x * x / (two+ one);//at abs(x)=0.001 the relative error of this approximation is about 2e-13
    }else{
      re=atan(x);


#if DEBUG >= CALCULATION
      cout << __FILE__ << " : Plain arctan is " << re <<endl;
      cout << __FILE__ << " : real(arctan) / pi *2 + 1.0 is " << real(re) / complex_type::PI * 2.0 +1.0<<endl;
#endif

      //choose 'right' branch: ... ?
      /*
        while(real(re)>complex_type::PI/2.0){
        re-= complex_type::PI;
        }
        while(real(re)<=-complex_type::PI/2.0 + NumTraits<complex_type>::epsilon()){
        re+= complex_type::PI;
        }
      */

      //  re+= complex_type::PI;

#if DEBUG >=  CALCULATION
      cout << __FILE__ << " : Corrected arctan is " << re <<endl;
#endif

      re/=x;
    }

#if DEBUG >=  CALCULATION
    cout << __FILE__ << " : arctan(x)/x =  " << re <<endl;
#endif


    return re;
  }



  /// The actual computation is done at changePoint
  template<typename realT>
  typename SCE<realT>::value_type SCE<realT>::calcF(){
    if(!dBMSet)
      changePoint(dBM);//make sure internals are up to date


    if(real(p.g)<-1e-20){
      return realT(1e10);
      //    return mpfr::const_infinity(); //inf quickly leads to nan, nan leads to errors, errors lead to suffering... ;)
    }

    //  return BM;
    //    return BMADD;
    //return p.g;
    //  return dBS;
    // return a + two * p.dim * b;//C
    // return - two * b; //D
    //  return a;
    // return b;
    //    return Int;

    //enforce positive density of states. (There seems to be another solution for \f$ - \bar g \f$.)

    return Int * p.z - p.g;

  }



  /**
****************************************************************************
* (Approximation of)
* \f[
*    Int  = \int_{[0,2\pi]^d} \frac{d^{d} k}{(2 \pi)^d} \frac{1}{a + b 2(d - \sum_{i=1}^d \cos(k_i))} \approx \frac{d}{\Omega^d} \int_0^\Omega \frac{k^{d-1} d k}{a + b k^2}
*  \f]
*
****************************************************************************
*/
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

    complex_type X = - p.cutOff * p.cutOff / a * b;

    if(oddDim){

#if DEBUG >= DETAIL
      cout << __FILE__ << " : Using Debeye approximation for integrals in odd dimension d = " << p.dim;
      cout << " = 2 *" << n << " + 1"<<endl;
#endif

      // \f[
      //   Int  = \int_0^\Omega \frac{d^{2n} k}{a + b k^2}
      // \f]

      complex_type rt= sqrt(-X);

#if DEBUG >= CALCULATION
      cout << __FILE__ << " : X = " << X << " => sqrt(-X) = " << rt << endl;
#endif

      Int = -X * arctanOverX(rt);

      complex_type pow(1.0);
      complex_type powAB(1.0);

      for(realT k=1.0; k<=n;k++){
        pow*=X;
        Int +=  pow / (two * k - one);
        powAB *= - a / b;
      }

      Int *= powAB / p.cutOff;

    }else{

#if DEBUG >= DETAIL
      cout << __FILE__ << " : Using Debeye approximation for integrals in even dimension " << p.dim <<endl;
      cout << " = 2 *" << n << " + 2"<<endl;
#endif

      // \f[
      //   Int  = \int_0^\Omega \frac{d^{2n + 1} k}{a + b k^2}
      // \f]

      Int = log(one - X);

      complex_type pow(1.0);
      complex_type powAB(1.0);

      if((int)(toDouble<realT>(n)) % 2 ==1){
#if DEBUG >= DETAIL
        cout << __FILE__ << " : n = " << n << " is odd." <<endl;
#endif
        pow=-1.0;
      }

      for(realT k=1.0;k<=n;k++){
        pow *= X;
        Int += pow / k;
        powAB *= a / b;
      }
      Int *= powAB / two;
    }


    Int *= p.dim / p.cutOffD / b;

#if DEBUG >= DETAIL
    cout << __FILE__ << " : The integral was computed to be " << Int <<endl;
#endif

  }





  template<typename realT>
  typename SCE<realT>::derivative_type SCE<realT>::calcJ(){

    if(real(p.g)<-1e-20){
      return realT(0.0);
    }


    complex_type DBM = p.b.mass.mix;

    complex_type DBMADD, DZM, Dg;

    if(p.b.mass.add>0.0){
      //      BMADD =  BM / (BM - p.b.mass.add * dBM / p.am);
      DBMADD = (DBM * (BM - p.b.mass.add * dBM / p.am) - (DBM - p.b.mass.add / p.am) * BM) * BMADD * BMADD / BM / BM;

      //      ZM = BM * (BM + p.b.mass.mix) + p.b.mass.add * BMADD;
      DZM = DBM * (two * BM + p.b.mass.mix) + p.b.mass.add * DBMADD;

      //      p.g = (BM * (one + dBM) + p.am * (BMADD - one) )/p.z;
      Dg = (BM + DBM*(one+dBM) + p.am * DBMADD ) / p.z ;
    }else{
      DZM = DBM * (two * BM + p.b.mass.mix);
      Dg = (BM + DBM*(one+dBM)) / p.z ;
    }

    complex_type DBS = - BS * BS * p.b.spring.mix / p.z * (Dg * ZM + DZM * p.g);

    complex_type DBSADD,Da;

    if(p.b.spring.add>0.0){
      DBSADD = - BSADD * BSADD * p.b.spring.add / p.as / p.z * (Dg * ZM + DZM * p.g);
      Da = Dg * p.z * p.b.spring.mix *ZM + DZM*p.b.spring.mix * (p.dim - one + p.g * p.z) + p.b.spring.add * ZM * p.dim * (DBSADD * BS - DBS * BSADD) / BS / BS + p.dim * p.b.spring.add * DZM * BSADD/BS;
    }else{
      Da = Dg * p.z * p.b.spring.mix * ZM + DZM*p.b.spring.mix * (p.dim - one + p.g * p.z);
    }

    complex_type Db = DBS*ZM + BS*DZM;

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

      complex_type X = - p.cutOff * p.cutOff / a * b;
      complex_type DX = X * (Db / b - Da / a);

      if(oddDim){

#if DEBUG >= DETAIL
        cout << __FILE__ << " : Using Debeye approximation for derivatives of integrals in odd dimension " << p.dim <<endl;
#endif

        complex_type rt = sqrt(-X);

        DInt = -X*((one / (one - X) - arctanOverX(rt)) * DX / X / two) - arctanOverX(rt)*DX ;

        complex_type pow(1.0);
        complex_type powAB(1.0);
        for(realT k=1.0;k<=n;k++){
          DInt +=  k * pow / (two * k - one) *DX;
          powAB *= - a / b;
          pow *= X;
        }

        DInt *= powAB /p.cutOff;

        DInt *= p.dim / p.cutOffD /b;

        DInt += Int * (n * Da / a - (n+one) * Db / b);

      }else{

#if DEBUG >= DETAIL
        cout << __FILE__ << " : Using Debeye approximation for derivatives of integrals in even dimension " << p.dim <<endl;
#endif


        DInt = one / (X - one);

        complex_type pow(1.0);
        complex_type powAB(1.0);

        if((int)(toDouble<realT>(n)) % 2 == 1){
          pow=-1.0;
        }

        for(realT k=1.0;k<=n;k++){
          DInt += pow ;
          pow *= X;

          powAB *= a / b;
        }

        DInt *= powAB/two *DX;
        DInt *= p.dim / p.cutOffD /b;
        DInt +=  Int *( n * Da/a - (n+one)* Db/b);
      }


#if USEEXACTINT>0
    }
#endif



    //   return DBM;
    //    return DBMADD;
    //  return Dg;
    //  return DdBS;
    //  return Da;
    // return Db;
    //    return DInt;

     return DInt * p.z - Dg;

  }


  /// The deterministic limit is the well known harmonic crystal
  template<typename realT>
  typename SCE<realT>::value_type SCE<realT>::guessStartPoint(){


    a = p.z*p.z;
    b = one;

    computeInt();

    p.g = Int * p.z;

    dBM = p.g * p.z -one;

    changePoint(dBM);


#if DEBUG >= DETAIL
    cout << __FILE__ << " : Deterministic DOS at z = " << p.z << " is " <<p.g <<endl;
#endif


    if(real(p.g)<0.0)
      throw std::runtime_error(string(__FILE__)+string(" : exact starting point must not have real(p.g) < 0 !"));





    if(abs(calcF())>1e-8)
      throw std::runtime_error(string(__FILE__)+string(" : exact starting point must be a solution!"));


    return dBM;
  }



}//end namespace

#endif
