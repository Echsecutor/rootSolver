/**
 * @file complex.hpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0.2014-06-11
 *
 * @section DESCRIPTION
 *
 *
 * The trouble with with std::complex<T> is
 * (according to http://www.cplusplus.com/reference/complex/complex/ ):
 *
 * The effect of instantiating a complex with a T other than
 * float, double or long double is undefined (certain library
 * implementations may support it, but the resulting code is
 * non-portable).
 *
 * What actually happens within my setup is that accuracy is lost
 * unpredictably at various arithmetic operations... o0
 *
 * Googling for a while did not produce a complex high precision
 * solution for Eigen, hence here is my own implementation.
 *
 * (Notice: https://stackoverflow.com/questions/1759300/when-should-i-write-the-keyword-inline-for-a-function-method )
 *
 *
 * This class can be used like
 * #include <mpreal.h>
 * using namespace mpfr;
 * typedef mpreal real_type;
 * typedef root_solver::complex<real_type> complex_type;
 * int main(int args, char *arg[]){
 * unsigned int precision=256;
 * mpreal::set_default_prec(precision);
 * root_solver::complex<real_type>::PI = const_pi(precision, mpreal::get_default_rnd());
 * complex_type I(0.0,1.0);
 * ...
 *
 * See e.g. derivativeVerification.cpp for a working example.
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


#ifndef COMPLEX_HPP
#define COMPLEX_HPP

#include <stdexcept>
#include <eigen3/Eigen/Core>
#include "util.hpp"

#include <cmath>//for std::sqrt and log


namespace root_solver{
  /** Use MatrixBase as a trait to obtain
   *
   * typedef Eigen::MatrixBase<EigenStuff<root_solver::complex<_real_type>>>::RealScalar real_type;
   *
   * and also
   *
   * typedef Eigen::MatrixBase<root_solver::complex<_real_type>>::RealScalar real_type;
   *
   * Notice that NumTraits<...>::Real only works in the
   * second case but produces a matrix with real entries in the first.
   *
   */
  using Eigen::MatrixBase;

  //only for std types... hopefully ;)
  using std::sqrt;
  using std::log;


  template<typename realT>
  class complex{
  private :
    realT RE;
    realT IM;

  public:

    typedef realT real_type;

    static real_type PI;
    static const complex<real_type> I;

    complex():RE(0.0),IM(0.0){}
    complex(realT re):RE(re),IM(0.0){}
    complex(int re):RE(double(re)),IM(0.0){} // e.g. mpreal can not be constructed from int by default...? o0
    complex(realT re, realT im):RE(re),IM(im){}

    realT real() const{
      return RE;
    }

    realT imag() const{
      return IM;
    }


    complex<realT> conjugate() const{
      return complex<realT>(real(), - imag());
    }

    realT abs2() const{
      return real()*real() + imag()*imag();
    }

    realT abs() const{
      return sqrt(abs2());
    }

    realT arg() const{
      if(IM==0.0 && RE==0.0)
        return 0.0;
      return atan2(IM,RE);
    }


    complex<realT> operator-() const{
      return complex<realT>(-RE, -IM);
    }

    complex<realT> operator +(const complex<realT>& rhs) const{
      return complex<realT>(RE + rhs.RE, IM + rhs.IM);
    }

    complex<realT> operator +(const realT& rhs) const{
      return complex<realT>(RE + rhs , IM);
    }

    complex<realT> operator -(const complex<realT>& rhs) const{
      return complex<realT>(RE - rhs.RE, IM - rhs.IM);
    }

    complex<realT> operator -(const realT& rhs) const{
      return complex<realT>(RE - rhs , IM);
    }

    complex<realT> operator *(const complex<realT>& rhs) const{
      return complex<realT>(RE * rhs.RE - IM * rhs.IM,  RE * rhs.IM +  IM * rhs.RE );
    }

    complex<realT> operator *(const realT& rhs) const{
      return complex<realT>(RE * rhs,  IM * rhs);
    }

    complex<realT> operator /(const complex<realT>& rhs) const{
      return (*this) * rhs.conjugate() / rhs.abs2();
    }

    complex<realT> operator /(const realT& rhs) const{
      return complex<realT>(RE / rhs,  IM / rhs);
    }


    complex<realT>& operator+=(const complex<realT>& rhs){
      RE += rhs.RE;
      IM += rhs.IM;
      return *this;
    }

    complex<realT>& operator+=(const realT& rhs){
      RE += rhs;
      return *this;
    }

    complex<realT>& operator-=(const complex<realT>& rhs){
      RE -= rhs.RE;
      IM -= rhs.IM;
      return *this;
    }

    complex<realT>& operator-=(const realT& rhs){
      RE -= rhs;
      return *this;
    }

    complex<realT>& operator*=(const complex<realT>& rhs){
      *this = *this * rhs;
      return *this;
    }

    complex<realT>& operator*=(const realT& rhs){
      RE = RE * rhs;
      IM = IM * rhs;
      return *this;
    }

    complex<realT>& operator/=(const complex<realT>& rhs){
      *this = *this / rhs;
      return *this;
    }

    complex<realT>& operator/=(const realT& rhs){
      RE = RE / rhs;
      IM = IM / rhs;
      return *this;
    }


    template<typename castable>
    complex<realT>& operator=(const castable& other)
    {
      if ((void *)this != (void *)&other) {
        *this = complex<realT>(other);
      }
      return *this;
    }



    complex<realT>& operator=(const complex<realT>& other) // copy assignment
    {
      if (this != &other) {
        RE=other.RE;
        IM=other.IM;
      }
      return *this;
    }


    bool operator==(const complex<realT>& rhs) {
      return RE == rhs.RE && IM == rhs.IM;
    }

    bool operator!=(const complex<realT>& rhs) {
      return *this != rhs;
    }


  };



  //---------------------------------------------------------------------------------------------------------




  //to be set at run time, like for mpfr
  template<typename realT>
  realT complex<realT>::PI=0.0;


  template<typename realT>
  const complex<realT> complex<realT>::I = complex<real_type>(0.0,1.0);

  template<typename realT>
  std::ostream& operator<<(std::ostream& os, const complex<realT>& obj)
  {
    //    os << "(" << obj.real() << ", " << obj.imag() << ")"; //std format
    os << obj.real() << "\t" << obj.imag(); //dat format
    return os;
  }

  template<typename realT>
  complex<realT> operator +(const realT& lhs, const complex<realT>& rhs){
    return rhs + lhs;
  }

  template<typename realT>
  complex<realT> operator -(const realT& lhs, const complex<realT>& rhs){
    return -rhs + lhs;
  }

  template<typename realT>
  complex<realT> operator *(const realT& lhs, const complex<realT>& rhs){
    return rhs * lhs;
  }

  template<typename realT>
  complex<realT> operator /(const realT& lhs, const complex<realT>& rhs){
    return rhs.conjugate() * lhs / rhs.abs2();
  }


  //c style:

  template<typename realT>
  realT real(complex<realT> z){
    return z.real();
  }

  template<typename realT>
  realT imag(complex<realT> z){
    return z.imag();
  }

  template<typename realT>
  realT abs(complex<realT> z){
    return z.abs();
  }

  template<typename realT>
  realT abs2(complex<realT> z){
    return z.abs2();
  }

  template<typename realT>
  realT norm(complex<realT> z){
    return z.abs();
  }

  template<typename realT>
  realT arg(complex<realT> x){
    return x.arg();
  }

  template<typename realT>
  complex<realT> conj (const complex<realT>& x){
    return x.conjugate();
  }


  //some important functions
  //todo: increase precision internally for mpfr types

  /// computes exp(i x)
  template<typename realT>
  complex<realT> fromAngle(realT x){
    return complex<realT> (cos(x), sin(x));
  }

  template<typename realT>
  complex<realT> exp(complex<realT> x){
    return exp(x.real()) * fromAngle(x.imag());
  }


  /// branch cut inherited from atan implementation for realT
  template<typename realT>
  complex<realT> sqrt(complex<realT> x){
    return fromAngle(arg(x) / realT(2.0)) * sqrt(x.abs());
  }


  template<typename realT>
  complex<realT> log(complex<realT> x){
    realT re = log(x.abs());
    realT im = arg(x);
    return complex<realT>(re,im);
  }

  //using log
  template<typename realT>
  complex<realT> atan(complex<realT> x){
    return complex<realT>::I * (log(realT(1.0) - complex<realT>::I * x) - log(realT(1.0) + complex<realT>::I * x)) / realT(2.0);
  }



}//end namespace root_solver



// Explaining to Eigen that this is a complex type.
// I hope, following the std conventions, i.e.
// adding specialisations for my own types is ok. ;)
namespace Eigen{
  using root_solver::complex;

  template<typename _Real> struct NumTraits<root_solver::complex<_Real> >
  : GenericNumTraits<root_solver::complex<_Real> >
  {
    typedef _Real Real;
    enum {
      IsComplex = 1,
      RequireInitialization = NumTraits<_Real>::RequireInitialization,
      ReadCost = 2 * NumTraits<_Real>::ReadCost,
      AddCost = 2 * NumTraits<Real>::AddCost,
      MulCost = 4 * NumTraits<Real>::MulCost + 2 * NumTraits<Real>::AddCost
    };

    inline static Real epsilon() { return NumTraits<Real>::epsilon(); }
    inline static Real dummy_precision() { return NumTraits<Real>::dummy_precision(); }
  };



  namespace internal{

    template<typename _Real>
    struct traits<root_solver::complex<_Real> > : traits<Matrix<root_solver::complex<_Real>, 1, 1> >{};


    template<typename T> struct scalar_product_traits<T,root_solver::complex<T> >
    {
      //enum { Cost = 2*NumTraits<T>::MulCost };
      typedef root_solver::complex<T> ReturnType;
    };

    template<typename T> struct scalar_product_traits<root_solver::complex<T>, T>
    {
      //enum { Cost = 2*NumTraits<T>::MulCost  };
      typedef root_solver::complex<T> ReturnType;
    };

    /****************************************************************************
     * Implementation of abs
     * todo:needed?
     ****************************************************************************/


    template<typename RealScalar>
    struct abs_impl<root_solver::complex<RealScalar> >
    {
      static inline RealScalar run(const root_solver::complex<RealScalar>& x)
      {
        return x.abs();
      }
    };

    /****************************************************************************
     * Implementation of abs2
     ****************************************************************************/


    template<typename RealScalar>
    struct abs2_impl<root_solver::complex<RealScalar> >
    {
      static inline RealScalar run(const root_solver::complex<RealScalar>& x)
      {
        return x.abs2();
      }
    };



    /****************************************************************************
     * Implementation of real                                                 *
     ****************************************************************************/

    template<typename RealScalar>
    struct real_impl<root_solver::complex<RealScalar> >
    {
      static inline RealScalar run(const root_solver::complex<RealScalar>& x)
      {
        return x.real();
      }
    };



    /****************************************************************************
     * Implementation of imag                                                 *
     ****************************************************************************/

    template<typename RealScalar>
    struct imag_impl<root_solver::complex<RealScalar> >
    {
      static inline RealScalar run(const root_solver::complex<RealScalar>& x)
      {
        return x.imag();
      }
    };


    /****************************************************************************
     * Implementation of conj                                             *
     ****************************************************************************/

    template<typename RealScalar>
    struct conj_impl<root_solver::complex<RealScalar> >
    {
      static inline root_solver::complex<RealScalar> run(const root_solver::complex<RealScalar>& x)
      {
        return x.conjugate();
      }
    };

  }
}


#endif

