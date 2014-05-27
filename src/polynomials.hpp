/**
 * @file polynomials.hpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0.2014-05-21
 *
 * @section DESCRIPTION
 *
 * The polynomials class generates multidimensional polynomial functions
 * as trial functions for the solvers.
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




#ifndef POLYNOMIALS_HPP
#define POLYNOMIALS_HPP

#include <ostream>

#include <complex>
#include <iterator>
#include <memory>
using namespace std;


#ifndef COMPLEX_I
#define COMPLEX_I
static complex<double> I(0.0,1.0);
#endif


#include <eigen3/Eigen/Dense>
using namespace Eigen;


#include "singleRootSolver.hpp"
#include "multiRootSolver.hpp"
#include "batchSolver.hpp"
#include "extrapolationSolver.hpp"



class polyParams{
public:
  int dim;
  int maxDegree;
  int * m;
  complex<double> ** s;
  int *** degree;
  complex<double> **** zeros;

  friend ostream& operator<< (ostream &out, polyParams &p){
    out << p.dim <<"\t" << p.maxDegree;
    return out;
  }

  
};




template <int dim>
class polynomials : public virtual MRS_functions<complex<double>,dim>{
protected:

  int functionCalls;
  int jacobianCalls;

  // using a pointer instead of aligning this whole class with  EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)
  shared_ptr<Matrix<complex<double>, dim, 1>> x;

  polyParams * para;

  complex<double> p(complex<double> z, int nZeros, complex<double> * zeros){
    complex<double> re=1.0;
    for(int i=0;i<nZeros;i++){
      re *= (z-zeros[i]);
    }
    return re;
  }

  complex<double> P(const Matrix<complex<double>,dim,1> z, int m,complex<double> *s,int ** degree, complex<double> *** zeros){
    complex<double> re=0.0;

    for(int i=0;i<m;i++){
      complex<double> prod=s[i];
      for(int k=0;k<dim;k++){
        prod*=p(z(k),degree[i][k],zeros[i][k]);
      }
      re +=prod;
    }
    return re;
  }

  Matrix<complex<double>,dim,1> Poly(const Matrix<complex<double>,dim,1> z, polyParams * p){
    Matrix<complex<double>,dim,1> re;
    for(int j=0;j<dim;j++){
      complex<double> Pj=P(z,p->m[j],p->s[j],p->degree[j],p->zeros[j]);
      re(j)=Pj;
    }
    return re;
  }


  // derivatives of the polynomial function:

  complex<double> dp(complex<double> z, int nZeros, complex<double> * zeros){
    complex<double> re=0.0;
    for(int j=0;j<nZeros;j++){
      complex<double> prod=1.0;
      for(int i=0;i<nZeros;i++){
        if(i!=j){
          prod *= (z-zeros[i]);
        }
      }
      re+=prod;
    }

    return re;
  }

  complex<double> dPdzj(int j,const Matrix<complex<double>,dim,1> z, int m,complex<double> *s,int ** degree, complex<double> *** zeros){
    complex<double> re=0.0;

    for(int i=0;i<m;i++){
      complex<double> prod=s[i];
      for(int k=0;k<dim;k++){
        if(k==j){
          prod*=dp(z(k),degree[i][k],zeros[i][k]);
        }else{
          prod*=p(z(k),degree[i][k],zeros[i][k]);
        }
      }
      re +=prod;
    }
    return re;
  }

  Matrix<complex<double>,dim,dim> dPolyDz(const Matrix<complex<double>,dim,1> z, polyParams * p){
    Matrix<complex<double>,dim,dim> re;
    for(int j=0;j<dim;j++){
      for(int i=0; i<dim; i++){
        complex<double> Pj=dPdzj(i,z,p->m[j],p->s[j],p->degree[j],p->zeros[j]);
        re(j,i)=Pj;
      }
    }
    return re;
  }

  double absF0;

public:

  polynomials():functionCalls(0),jacobianCalls(0),x(new Matrix<complex<double>, dim, 1>()),para(new polyParams()),absF0(0){
    para->maxDegree=dim*3;
    randomiseTestcase();
}

  polynomials(polyParams * p_):functionCalls(0),jacobianCalls(0),x(new   Matrix<complex<double>, dim, 1>()),para(p_),absF0(0){}

  polynomials(const polynomials& p2):functionCalls(p2.functionCalls),jacobianCalls(p2.jacobianCalls),x(p2.x),para(p2.para),absF0(p2.absF0){}


  polyParams * getParameters(){return para;}
  void setParameters(polyParams * p){para=p;}

  Matrix<complex<double>, dim, 1> calcF(){
    if(absF0==0){absF0=Poly(*x,para).norm();}
    functionCalls++;
    return Poly(*x, para)/absF0;
  }

  Matrix<complex<double>, dim, dim> calcJ(){
    if(absF0==0){absF0=Poly(*x,para).norm();}
    jacobianCalls++;
    return dPolyDz(*x,para)/absF0;
  }

  void changePoint(Matrix<complex<double>, dim, 1> x_){
    *x = x_;
  }

  int getFunctionCallsCounter(){return functionCalls;}
  int getJacobianCallsCounter(){return jacobianCalls;}


  /// generates random zeros and prefactors defining the polynomials
  void randomiseTestcase(){
    int maxPols=dim;
    int maxDegree=para->maxDegree;

    complex<double> d;//dummy
    para->dim=dim;
    para->m = new int[dim];
    para->s= new complex<double>*[dim];
    para->degree= new int **[dim];
    para->zeros = new complex<double> ***[dim];

    for(int j=0;j<dim;j++){
      para->m[j]= (rand() % (maxPols))+1;
      para->s[j]=new complex<double>[para->m[j]];
      para->degree[j]=new int *[para->m[j]];
      para->zeros[j] = new complex<double> **[para->m[j]];
      for(int i=0;i<para->m[j];i++){
        int totalDegreeLeft=maxDegree;
        double x= ((double)rand()/(double)RAND_MAX);
        d = x;
        para->s[j][i]=d;
        para->degree[j][i] = new int [dim];
        para->zeros[j][i] = new complex<double> *[dim];
        for(int k=0;k<dim;k++){
          para->degree[j][i][k]=rand() % (totalDegreeLeft+1);
          totalDegreeLeft -= para->degree[j][i][k];
          para->zeros[j][i][k] = new complex<double>[para->degree[j][i][k]];
          for(int n=0;n<para->degree[j][i][k];n++){
            double x= 2.0*((double)rand()/(double)RAND_MAX)-1.0;
            double y= 2.0*((double)rand()/(double)RAND_MAX)-1.0;
            d= x + I * y;
            para->zeros[j][i][k][n]=d;
          }
        }
      }
    }
  }

};



/// wraps the above to a SRS_function.
// Yeah I know, implementing it the other way around would have made more sense... clean it up if you wish. ;)
class polynomial : public virtual SRS_functions<complex<double>>{
 protected:
  polynomials<1> ps;

 public:

  polynomial():ps(){}
  polynomial(polyParams * p_):ps(p_){}

  virtual complex<double>calcF(){
    return ps.calcF()(0);
  }

  virtual complex<double> calcJ(){
    return ps.calcJ()(0,0);
  }

  virtual void changePoint(const complex<double> x){
    Matrix<complex<double>, 1, 1> z;
    z << x;
    ps.changePoint(z);
  }

  int getFunctionCallsCounter(){return ps.getFunctionCallsCounter();}
  int getJacobianCallsCounter(){return ps.getJacobianCallsCounter();}

  void randomiseTestcase(){
    ps.randomiseTestcase();
  }

  polyParams * getParameters(){return ps.getParameters();}



};



//---------------------------------------------------------------------------------------------
// Batch polynomials:

///parameterSetIterator

template <int dim>
class PSI : public iterator<std::input_iterator_tag, polyParams>{
private:
  shared_ptr<polynomials<dim>> p;
public:
  int counter;
  PSI(polyParams* x):p(new polynomials<dim>(x)),counter(0){}
  PSI(const PSI& p2):p(new polynomials<dim>(*(p2.p))),counter(p2.counter){}
  PSI& operator++() {this->p->randomiseTestcase();this->counter++;return *this;}
  PSI operator++(int) {PSI tmp(*this); operator++(); return tmp;}
  bool operator==(const PSI& rhs) {return this->counter == rhs.counter;}
  bool operator!=(const PSI& rhs) {return !operator==(rhs);}
  polyParams& operator*() {return *(this->p->getParameters());}
};



template <int dim>
class batch_polynomials :public virtual batch_functions<multiRootSolver<complex<double>,dim>, PSI<dim>>, public virtual polynomials<dim>{
 public:

  batch_polynomials():polynomials<dim>(){}
  batch_polynomials(polyParams * p_):polynomials<dim>(p_){}
  batch_polynomials(const batch_polynomials& p2):polynomials<dim>(p2){}

  virtual void setParameters(PSI<dim> It){
    this->para = &(*It);
  }
};



class batch_polynomial :public virtual batch_functions<singleRootSolver<complex<double>>, PSI<1>>, public virtual polynomial{
 public:

  batch_polynomial():polynomial(){}
  batch_polynomial(polyParams * p_):polynomial(p_){}
  batch_polynomial(const batch_polynomial& p2):polynomial(p2){}

  virtual void setParameters(PSI<1> It){
    this->ps.setParameters(&(*It));
  }
};


//---------------------------------------------------------------------------------------------
// trivial polynomial shift to test extrapolation:

class extra_polynomial : public virtual extra_functions<singleRootSolver<complex<double> > >, public virtual polynomial{//you can avoid diamonds if you are not using the same function for dozens of solvers or are less lazy ;)

private:
  double shift;

public:

  extra_polynomial():polynomial(),shift(0){}
  extra_polynomial(polyParams * p):polynomial(p),shift(0){}
  extra_polynomial(const extra_polynomial& p2):polynomial(p2), shift(0){}

  virtual double get_extra_parameter(){return shift;}
  virtual void set_extra_parameter(double p){shift = p;}

  virtual double get_initial(){return 0;}
  virtual double get_final(){return 1;}

  virtual double get_max_change(){return 1.0/8.0;}


  virtual void changePoint(complex<double> x){
    int shiftOrder=4;
    double c[] = {0.1,-1.0,1.0,-1.0,1.0};

    //    polynomial::changePoint(x + shift);
    complex<double>shifted = 0.0;
    for(int i=0;i<=shiftOrder;i++){
      shifted *=shift;
      shifted+=c[i];
    }
    shifted += x;
    polynomial::changePoint(shifted);
  }

};

template <int dim>
class extra_polynomials : public virtual extra_functions<multiRootSolver<complex<double>,dim> >, public virtual polynomials<dim>{

private:
  double shift;

public:

  extra_polynomials():polynomials<dim>(),shift(0){}
  extra_polynomials(polyParams * p):polynomials<dim>(p),shift(0){}
  extra_polynomials(const extra_polynomial& p2):polynomials<dim>(p2), shift(0){}

  virtual double get_extra_parameter(){return shift;}
  virtual void set_extra_parameter(double p){shift = p;}

  virtual double get_initial(){return 0;}
  virtual double get_final(){return 1;}

  virtual double get_max_change(){return 1.0/8.0;}

  virtual void changePoint(Matrix<complex<double>, dim,1> x){
    int shiftOrder=4;
    double c[] = {0.1,-1.0,1.0,-1.0,1.0};

    //    polynomial::changePoint(x + shift);
    complex<double>shifted = 0.0;
    for(int i=0;i<=shiftOrder;i++){
      shifted *=shift;
      shifted+=c[i];
    }
    for(int i=0; i<dim;i++){
      x[i]+=shifted;
    }
    polynomials<dim>::changePoint(x);
  }

};

//---------------------------------------------------------------------------------------------
// combine shift and batch


class batch_extra_polynomial : public virtual extra_polynomial, public virtual batch_functions<extraSolver<singleRootSolver<complex<double> > >, PSI<1> >
{
 public:

  batch_extra_polynomial():extra_polynomial(){}
  batch_extra_polynomial(polyParams * p_):extra_polynomial(p_){}
  batch_extra_polynomial(const batch_extra_polynomial& p2):extra_polynomial(p2){}

  virtual void setParameters(PSI<1> It){
    this->ps.setParameters(&(*It));
  }

};


template<int dim>
  class batch_extra_polynomials :public virtual extra_polynomials<dim>, public virtual batch_functions<extraSolver<multiRootSolver<complex<double>,dim> >, PSI<dim> >{
 public:

  batch_extra_polynomials():extra_polynomials<dim>(){}
  batch_extra_polynomials(polyParams * p_):extra_polynomials<dim>(p_){}
  batch_extra_polynomials(const batch_extra_polynomial& p2):extra_polynomials<dim>(p2){}

  virtual void setParameters(PSI<dim> It){
    this->para = &(*It);
  }

};



#endif

