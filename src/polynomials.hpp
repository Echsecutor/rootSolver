/**
 * @file polynomials.hpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0.2014-06-11
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

#include <iterator>
#include <memory>
using namespace std;

#include <eigen3/Eigen/Dense>
using namespace Eigen;

#include "complex.hpp"
#include "singleRootSolver.hpp"
#include "multiRootSolver.hpp"
#include "batchSolver.hpp"
#include "extrapolationSolver.hpp"

using namespace root_solver;

template<typename realT=double>
class polyParams{
public:
  int dim;
  int maxDegree;
  int * m;
  root_solver::complex<realT> ** s;
  int *** degree;
  root_solver::complex<realT> **** zeros;

  friend ostream& operator<< (ostream &out, polyParams<realT> &p){
    out << p.dim <<"\t" << p.maxDegree;
    return out;
  }


};




template <int dim, typename realT=double>
class polynomials : public virtual MRS_functions<root_solver::complex<realT>,dim>{
protected:

  int functionCalls;
  int jacobianCalls;

  // using a pointer instead of aligning this whole class with  EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)
  shared_ptr<Matrix<root_solver::complex<realT>, dim, 1>> x;

  polyParams<realT> * para;

  root_solver::complex<realT> p(root_solver::complex<realT> z, int nZeros, root_solver::complex<realT> * zeros){
    root_solver::complex<realT> re=1.0;
    for(int i=0;i<nZeros;i++){
      re *= (z-zeros[i]);
    }
    return re;
  }

  root_solver::complex<realT> P(const Matrix<root_solver::complex<realT>,dim,1> z, int m,root_solver::complex<realT> *s,int ** degree, root_solver::complex<realT> *** zeros){
    root_solver::complex<realT> re=0.0;

    for(int i=0;i<m;i++){
      root_solver::complex<realT> prod=s[i];
      for(int k=0;k<dim;k++){
        prod*=p(z(k),degree[i][k],zeros[i][k]);
      }
      re +=prod;
    }
    return re;
  }

  Matrix<root_solver::complex<realT>,dim,1> Poly(const Matrix<root_solver::complex<realT>,dim,1> z, polyParams<realT> * p){
    Matrix<root_solver::complex<realT>,dim,1> re;
    for(int j=0;j<dim;j++){
      root_solver::complex<realT> Pj=P(z,p->m[j],p->s[j],p->degree[j],p->zeros[j]);
      re(j)=Pj;
    }
    return re;
  }


  // derivatives of the polynomial function:

  root_solver::complex<realT> dp(root_solver::complex<realT> z, int nZeros, root_solver::complex<realT> * zeros){
    root_solver::complex<realT> re=0.0;
    for(int j=0;j<nZeros;j++){
      root_solver::complex<realT> prod=1.0;
      for(int i=0;i<nZeros;i++){
        if(i!=j){
          prod *= (z-zeros[i]);
        }
      }
      re+=prod;
    }

    return re;
  }

  root_solver::complex<realT> dPdzj(int j,const Matrix<root_solver::complex<realT>,dim,1> z, int m,root_solver::complex<realT> *s,int ** degree, root_solver::complex<realT> *** zeros){
    root_solver::complex<realT> re=0.0;

    for(int i=0;i<m;i++){
      root_solver::complex<realT> prod=s[i];
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

  Matrix<root_solver::complex<realT>,dim,dim> dPolyDz(const Matrix<root_solver::complex<realT>,dim,1> z, polyParams<realT> * p){
    Matrix<root_solver::complex<realT>,dim,dim> re;
    for(int j=0;j<dim;j++){
      for(int i=0; i<dim; i++){
        root_solver::complex<realT> Pj = dPdzj(i,z,p->m[j],p->s[j],p->degree[j],p->zeros[j]);
        re(j,i)=Pj;
      }
    }
    return re;
  }

  realT absF0;

public:

  polynomials():functionCalls(0),jacobianCalls(0),x(new Matrix<root_solver::complex<realT>, dim, 1>()),para(new polyParams<realT>()),absF0(0){
    para->maxDegree=dim*3;
    randomiseTestcase();
  }

  polynomials(polyParams<realT> * p_):functionCalls(0),jacobianCalls(0),x(new   Matrix<root_solver::complex<realT>, dim, 1>()),para(p_),absF0(0){}

  polynomials(const polynomials& p2):functionCalls(p2.functionCalls),jacobianCalls(p2.jacobianCalls),x(p2.x),para(p2.para),absF0(p2.absF0){}


  polyParams<realT> * getParameters(){return para;}
  void setParameters(polyParams<realT> * p){para=p;}

  Matrix<root_solver::complex<realT>, dim, 1> calcF(){
    if(absF0 == 0){absF0 = Poly(*x,para).norm();}
    functionCalls++;
    return Poly(*x, para)/absF0;
  }

  Matrix<root_solver::complex<realT>, dim, dim> calcJ(){
    if(absF0==0){absF0 = abs(Poly(*x,para).norm());}
    jacobianCalls++;
    return dPolyDz(*x,para)/absF0;
  }

  void changePoint(Matrix<root_solver::complex<realT>, dim, 1> x_){
    *x = x_;
  }

  int getFunctionCallsCounter(){return functionCalls;}
  int getJacobianCallsCounter(){return jacobianCalls;}


  /// generates random zeros and prefactors defining the polynomials
  void randomiseTestcase(){
    int maxPols=dim;
    int maxDegree=para->maxDegree;

    root_solver::complex<realT> d;//dummy
    para->dim=dim;
    para->m = new int[dim];
    para->s= new root_solver::complex<realT>*[dim];
    para->degree= new int **[dim];
    para->zeros = new root_solver::complex<realT> ***[dim];

    for(int j=0;j<dim;j++){
      para->m[j]= (rand() % (maxPols))+1;
      para->s[j]=new root_solver::complex<realT>[para->m[j]];
      para->degree[j]=new int *[para->m[j]];
      para->zeros[j] = new root_solver::complex<realT> **[para->m[j]];
      for(int i=0;i<para->m[j];i++){
        int totalDegreeLeft=maxDegree;
        realT x= ((realT)rand()/(realT)RAND_MAX);
        d = x;
        para->s[j][i]=d;
        para->degree[j][i] = new int [dim];
        para->zeros[j][i] = new root_solver::complex<realT> *[dim];
        for(int k=0;k<dim;k++){
          para->degree[j][i][k]=rand() % (totalDegreeLeft+1);
          totalDegreeLeft -= para->degree[j][i][k];
          para->zeros[j][i][k] = new root_solver::complex<realT>[para->degree[j][i][k]];
          for(int n=0;n<para->degree[j][i][k];n++){
            realT x= 2.0*((realT)rand()/(realT)RAND_MAX)-1.0;
            realT y= 2.0*((realT)rand()/(realT)RAND_MAX)-1.0;
            d = root_solver::complex<realT>(x,y);
            para->zeros[j][i][k][n]=d;
          }
        }
      }
    }
  }

};



/// wraps the above to a SRS_function.
// Yeah I know, implementing it the other way around would have made more sense... clean it up if you wish. ;)
template <typename realT=double>
class polynomial : public virtual SRS_functions<root_solver::complex<realT>>{
 protected:
  polynomials<1,realT> ps;

 public:

  polynomial():ps(){}
  polynomial(polyParams<realT> * p_):ps(p_){}

  virtual root_solver::complex<realT>calcF(){
    return ps.calcF()(0);
  }

  virtual root_solver::complex<realT> calcJ(){
    return ps.calcJ()(0,0);
  }

  virtual void changePoint(const root_solver::complex<realT> x){
    Matrix<root_solver::complex<realT>, 1, 1> z;
    z << x;
    ps.changePoint(z);
  }

  int getFunctionCallsCounter(){return ps.getFunctionCallsCounter();}
  int getJacobianCallsCounter(){return ps.getJacobianCallsCounter();}

  void randomiseTestcase(){
    ps.randomiseTestcase();
  }

  polyParams<realT> * getParameters(){return ps.getParameters();}



};



//---------------------------------------------------------------------------------------------
// Batch polynomials:

///parameterSetIterator

template <int dim,typename realT=double>
class PSI : public iterator<std::input_iterator_tag, polyParams<realT>>{
                            private:
  shared_ptr<polynomials<dim,realT>> p;
                            public:
  int counter;
  PSI(polyParams<realT>* x):p(new polynomials<dim,realT>(x)),counter(0){}
  PSI(const PSI& p2):p(new polynomials<dim,realT>(*(p2.p))),counter(p2.counter){}
  PSI& operator++() {this->p->randomiseTestcase();this->counter++;return *this;}
  PSI operator++(int) {PSI tmp(*this); operator++(); return tmp;}
  bool operator==(const PSI& rhs) {return this->counter == rhs.counter;}
  bool operator!=(const PSI& rhs) {return !operator==(rhs);}
  polyParams<realT>& operator*() {return *(this->p->getParameters());}
                            };



template <int dim,typename realT=double>
class batch_polynomials :public virtual batch_functions<multiRootSolver<root_solver::complex<realT>, dim >, PSI<dim,realT> >, public virtual polynomials<dim,realT>{
 public:

  batch_polynomials():polynomials<dim,realT>(){}
  batch_polynomials(polyParams<realT> * p_):polynomials<dim,realT>(p_){}
  batch_polynomials(const batch_polynomials<dim, realT>& p2):polynomials<dim,realT>(p2){}

  virtual void setParameters(PSI<dim,realT> It){
    this->para = &(*It);
  }
};


template<typename realT=double>
class batch_polynomial :public virtual batch_functions<singleRootSolver<realT>, PSI<1,realT>>, public virtual polynomial<realT>{
 public:

  batch_polynomial():polynomial<realT>(){}
  batch_polynomial(polyParams<realT> * p_):polynomial<realT>(p_){}
  batch_polynomial(const batch_polynomial<realT>& p2):polynomial<realT>(p2){}

  virtual void setParameters(PSI<1,realT> It){
    this->ps.setParameters(&(*It));
  }
};


//---------------------------------------------------------------------------------------------
// trivial polynomial shift to test extrapolation:

template<typename realT=double>
class extra_polynomial : public virtual extra_functions<singleRootSolver<realT>, realT>, public virtual polynomial<realT>{//you can avoid diamonds if you are not using the same function for dozens of solvers or are less lazy ;)

private:
  realT shift;

public:

  extra_polynomial():polynomial<realT>(),shift(0){}
  extra_polynomial(polyParams<realT> * p):polynomial<realT>(p),shift(0){}
  extra_polynomial(const extra_polynomial<realT>& p2):polynomial<realT>(p2), shift(0){}

  virtual realT get_extra_parameter(){return shift;}
  virtual void set_extra_parameter(realT p){shift = p;}

  virtual realT get_initial(){return 0;}
  virtual realT get_final(){return 1;}

  virtual realT get_max_change(){return 1.0/8.0;}


  virtual void changePoint(root_solver::complex<realT> x){
    int shiftOrder=4;
    realT c[] = {0.1,-1.0,1.0,-1.0,1.0};

    //    polynomial::changePoint(x + shift);
    root_solver::complex<realT>shifted = 0.0;
    for(int i=0;i<=shiftOrder;i++){
      shifted *=shift;
      shifted+=c[i];
    }
    shifted += x;
    polynomial<realT>::changePoint(shifted);
  }

};

template <int dim,typename realT=double>
class extra_polynomials : public virtual extra_functions<multiRootSolver<root_solver::complex<realT>,dim>,realT >, public virtual polynomials<dim>{

private:
  realT shift;

public:

  extra_polynomials():polynomials<dim>(),shift(0){}
  extra_polynomials(polyParams<realT> * p):polynomials<dim>(p),shift(0){}
  extra_polynomials(const extra_polynomial<realT>& p2):polynomials<dim,realT>(p2), shift(0){}

  virtual realT get_extra_parameter(){return shift;}
  virtual void set_extra_parameter(realT p){shift = p;}

  virtual realT get_initial(){return 0;}
  virtual realT get_final(){return 1;}

  virtual realT get_max_change(){return 1.0/8.0;}

  virtual void changePoint(Matrix<root_solver::complex<realT>, dim,1> x){
    int shiftOrder=4;
    realT c[] = {0.1,-1.0,1.0,-1.0,1.0};

    //    polynomial::changePoint(x + shift);
    root_solver::complex<realT>shifted = 0.0;
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

template<typename realT=double>
class batch_extra_polynomial : public virtual extra_polynomial<realT>, public virtual batch_functions<extraSolver<singleRootSolver<realT>,realT >, PSI<1,realT> >
{
public:

  batch_extra_polynomial():extra_polynomial<realT>(){}
  batch_extra_polynomial(polyParams<realT> * p_):extra_polynomial<realT>(p_){}
  batch_extra_polynomial(const batch_extra_polynomial<realT>& p2):extra_polynomial<realT>(p2){}

  virtual void setParameters(PSI<1,realT> It){
    this->ps.setParameters(&(*It));
  }

};


template<int dim,typename realT=double>
class batch_extra_polynomials :public virtual extra_polynomials<dim,realT>, public virtual batch_functions<extraSolver<multiRootSolver<root_solver::complex<realT>,dim>,realT >, PSI<dim,realT> >{
public:

  batch_extra_polynomials():extra_polynomials<dim,realT>(){}
  batch_extra_polynomials(polyParams<realT> * p_):extra_polynomials<dim,realT>(p_){}
  batch_extra_polynomials(const batch_extra_polynomial<realT>& p2):extra_polynomials<dim,realT>(p2){}

  virtual void setParameters(PSI<dim,realT> It){
    this->para = &(*It);
  }

};



#endif

