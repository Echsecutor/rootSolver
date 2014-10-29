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


//#include <algorithm>    // wrong std::swap ?
#include <utility>    // C++11 std::swap for arrays


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
class polyParams;

template <typename realT=double>
class PSI;

template <int dim,typename realT=double>
class polynomials;




template<typename realT>
class polyParams : virtual public root_solver::ostream_insertable{
public:
  int dim;
  int maxDegree;
  int * m;
  root_solver::complex<realT> ** s;
  int *** degree;
  root_solver::complex<realT> **** zeros;

  void swap(polyParams&p){
    std::swap(dim,p.dim);
    std::swap(maxDegree,p.maxDegree);

    cout << __FILE__ << " : before swap : m=" << m << ", p.m = " << p.m <<endl;
    std::swap(m,p.m);
    cout << __FILE__ << " : after swap : m=" << m << ", p.m = " << p.m <<endl;

    std::swap(s,p.s);
    std::swap(degree,p.degree);
    std::swap(zeros,p.zeros);
  }

  polyParams(){
    dim=0;
    maxDegree=0;
    m=0;
    s=0;
    degree=0;
    zeros=0;
  }

  polyParams(int dim_, int maxDegree_){
    dim=dim_;
    maxDegree=maxDegree_;
    m=0;
    s=0;
    degree=0;
    zeros=0;
    randomiseTestcase();
  }


  polyParams(const polyParams&p){
    //    cout << __FILE__ << " : calling polyParams const copy constructor"<<endl;
    dim = p.dim;
    maxDegree = p.maxDegree;
    m=0;
    s=0;
    degree=0;
    zeros=0;

    if(p.m!=0){
      m = new int[dim];
      s = new root_solver::complex<realT>*[dim];
      degree= new int **[dim];
      zeros = new root_solver::complex<realT> ***[dim];
      for(int j=0;j<dim;j++){
        m[j]= p.m[j];
        s[j]=new root_solver::complex<realT>[m[j]];
        degree[j]=new int *[m[j]];
        zeros[j] = new root_solver::complex<realT> **[m[j]];
        for(int i=0;i<m[j];i++){
          s[j][i]=p.s[j][i];
          degree[j][i] = new int [dim];
          zeros[j][i] = new root_solver::complex<realT> *[dim];
          for(int k=0;k<dim;k++){
            degree[j][i][k]=p.degree[j][i][k];
            if(degree[j][i][k]>0){
              zeros[j][i][k] = new root_solver::complex<realT>[degree[j][i][k]];
              for(int n=0; n<degree[j][i][k]; n++){
                zeros[j][i][k][n]=p.zeros[j][i][k][n];
              }
            }else{
              zeros[j][i][k]=0;
            }
          }
        }
      }
    }
    //    cout << __FILE__ << " : new " << *this << " copied from " << p << endl;
  }

  ~polyParams(){
    clear();
  }

  polyParams& operator=(const polyParams&p){
    polyParams tmp(p);
    swap(tmp);
    return *this;
  }

  virtual void insertMeInto(ostream &out) const{
    if(m==0){
      out << "polyParams: uninitialised!";
    }else{
      out << "polyParams: dim = " << dim << ", maxDegree = " << maxDegree;
      //      out << ", m = " << m <<endl;
      if(zeros[0][0][0] !=0){
        out << ", first zero: *" << zeros[0][0][0] << " = " << zeros[0][0][0][0];
      }else{
        cout << ", zeros[0][0][0] == 0 i.e. uninitialised.";
        cout << ", degree[0][0][0] = " << degree[0][0][0]<<endl;
      }
    }
  }

  void clear(){
    if(m==0)
      return;

    //    cout << __FILE__ << " : clearing polyParams " << *this <<endl;
    for(int j=0;j<dim;j++){
      for(int i=0;i<m[j];i++){
        for(int k=0;k<dim;k++){
	  //	  cout << "deleting zeros[" << j <<"][" << i <<"][" << k <<"]" <<endl;
	  //	  cout << zeros[j][i][k]<<endl;
          delete[] zeros[j][i][k];
        }
	//	cout << "deleting zeros and degree[" << j <<"][" << i <<"]" <<endl;
        delete[] degree[j][i];
        delete[] zeros[j][i];
      }
      //      cout << "deleting zeros and degree and s[" << j <<"]" <<endl;

      delete[] s[j];
      delete[] degree[j];
      delete[] zeros[j];
    }
    //    cout << "deleting top level" <<endl;

    delete[] m;
    m=0;
    delete[] s;
    s=0;
    delete[] degree;
    degree=0;
    delete[] zeros;
    zeros=0;
    //    cout << __FILE__ << " : cleared: " << *this << endl;
  }

  /// generates random zeros and prefactors defining the polynomials
  void randomiseTestcase(){
    if(dim==0){
      throw std::runtime_error("dim=0 can not be randomised!");
    }
    int maxPols=dim;
    root_solver::complex<realT> d;//dummy
    clear();
    m = new int[dim];
    s= new root_solver::complex<realT>*[dim];
    degree= new int **[dim];
    zeros = new root_solver::complex<realT> ***[dim];

    for(int j=0;j<dim;j++){
      m[j]= (rand() % (maxPols))+1;
      s[j]=new root_solver::complex<realT>[m[j]];
      degree[j]=new int *[m[j]];
      zeros[j] = new root_solver::complex<realT> **[m[j]];
      for(int i=0;i<m[j];i++){
        int totalDegreeLeft=maxDegree;
        realT x= ((realT)rand()/(realT)RAND_MAX);
        d = x;
        s[j][i]=d;
        degree[j][i] = new int [dim];
        zeros[j][i] = new root_solver::complex<realT> *[dim];
        for(int k=0;k<dim;k++){
          degree[j][i][k]=rand() % (totalDegreeLeft+1);
	  if(k==dim-1 && i == m[j]-1 && j == dim-1)
	    degree[j][i][k]=totalDegreeLeft;

          totalDegreeLeft -= degree[j][i][k];
          if(degree[j][i][k]>0){
            zeros[j][i][k] = new root_solver::complex<realT>[degree[j][i][k]];
            for(int n=0;n<degree[j][i][k];n++){
              realT x= 2.0*((realT)rand()/(realT)RAND_MAX)-1.0;
              realT y= 2.0*((realT)rand()/(realT)RAND_MAX)-1.0;
              d = root_solver::complex<realT>(x,y);
              zeros[j][i][k][n]=d;
            }
          }else{
            zeros[j][i][k]=0;
          }
        }
      }
    }

    cout << __FILE__ << " : after randomising : " << *this << endl;

  }




};



///parameterSetIterator

template <typename realT>
class PSI : public iterator<std::input_iterator_tag, polyParams<realT>>{
 private:
  polyParams<realT>* p;
 public:
  int counter;
  PSI(polyParams<realT>* x):p(x),counter(0){}
  PSI(const PSI& p2):p(p2.p),counter(p2.counter){}
  PSI& operator++() {this->p->randomiseTestcase();this->counter++;return *this;}
  PSI operator++(int) {PSI tmp(*this); operator++(); return tmp;}
  bool operator==(const PSI& rhs) {return this->counter == rhs.counter;}
  bool operator!=(const PSI& rhs) {return !operator==(rhs);}
  polyParams<realT> const& operator*() const{return *(this->p);}
};




template <int dim, typename realT>
class polyCore
{
protected:

  realT shift;

  int functionCalls;
  int jacobianCalls;

  // using a pointer instead of aligning this whole class with  EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)
  shared_ptr<Matrix<root_solver::complex<realT>, dim, 1>> x;

  polyParams<realT> para;

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

  Matrix<root_solver::complex<realT>,dim,1> Poly(const Matrix<root_solver::complex<realT>,dim,1> z, polyParams<realT>& p){
    Matrix<root_solver::complex<realT>,dim,1> re;
    for(int j=0;j<dim;j++){
      root_solver::complex<realT> Pj=P(z,p.m[j],p.s[j],p.degree[j],p.zeros[j]);
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

  Matrix<root_solver::complex<realT>,dim,dim> dPolyDz(const Matrix<root_solver::complex<realT>,dim,1> z, polyParams<realT>& p){
    Matrix<root_solver::complex<realT>,dim,dim> re;
    for(int j=0;j<dim;j++){
      for(int i=0; i<dim; i++){
        root_solver::complex<realT> Pj = dPdzj(i,z,p.m[j],p.s[j],p.degree[j],p.zeros[j]);
        re(j,i)=Pj;
      }
    }
    return re;
  }

  realT absF0;

public:

  polyCore():shift(0.0),functionCalls(0),jacobianCalls(0),x(new Matrix<root_solver::complex<realT>, dim, 1>()),para(),absF0(0){
    para.maxDegree=dim*3;
    randomiseTestcase();
  }

  polyCore(const polyParams<realT>& p_):shift(0.0),functionCalls(0),jacobianCalls(0),x(new   Matrix<root_solver::complex<realT>, dim, 1>()),para(p_),absF0(0){}

  polyCore(const polyCore& p2):shift(0.0),functionCalls(p2.functionCalls),jacobianCalls(p2.jacobianCalls),x(p2.x),para(p2.para),absF0(p2.absF0){}


  virtual void setParameters(const PSI<realT>& It){
    this->para = *It;
  }

  virtual polyCore* clone() const{
    return new polyCore(*this);
  }

  const polyParams<realT>& getParameters() const{return para;}

  //  const polyParams<realT> getConstParameters() const{return para;}

  void setParameters(polyParams<realT> p){para=p;}

  //////////////////////

  int getFunctionCallsCounter(){return functionCalls;}
  int getJacobianCallsCounter(){return jacobianCalls;}


  virtual realT get_extra_parameter(){return shift;}
  virtual void set_extra_parameter(realT p){shift = p;}

  virtual realT get_initial(){return 0;}
  virtual realT get_final(){return 1;}

  virtual realT get_max_change(){return 1.0/8.0;}


  Matrix<root_solver::complex<realT>, dim, 1> McalcF(){
    cout << __FILE__ << " : polyCore calculating F for " << para << ", at x= " << *x << endl;

    if(absF0 == 0){absF0 = Poly(*x,para).norm()+1e-9;
      cout << __FILE__ << " : polyCore setting absF0 = " << absF0 <<endl;
    }
    functionCalls++;

    cout << __FILE__ << " : polyCore calculated " << Poly(*x, para)/absF0 <<endl;

    return Poly(*x, para)/absF0;
  }

  Matrix<root_solver::complex<realT>, dim, dim> McalcJ(){
    if(absF0==0){absF0 = abs(Poly(*x,para).norm());}
    jacobianCalls++;
    return dPolyDz(*x,para)/absF0;
  }

  void MchangePoint(Matrix<root_solver::complex<realT>, dim, 1> x_){
    int shiftOrder=4;
    realT c[] = {0.1,-1.0,1.0,-1.0,1.0};

    //    polynomial::changePoint(x + shift);
    root_solver::complex<realT>shifted = 0.0;
    for(int i=0;i<=shiftOrder;i++){
      shifted *=shift;
      shifted+=c[i];
    }

    for(int i=0; i<dim;i++){
      (*x)[i]=shifted + x_[i];
    }

    cout << __FILE__ << " : polyCore changed point to " << *x<< ", with " << para << endl;
  }

  void randomiseTestcase(){
    para.dim=dim;
    para.randomiseTestcase();
  }

};




template <int dim, typename realT>
class polynomials : public virtual root_solver::batch_functions<root_solver::extraSolver<root_solver::multiRootSolver<root_solver::complex<realT>, dim>, realT>, PSI<realT> >,
                    public virtual root_solver::batch_functions<root_solver::multiRootSolver<root_solver::complex<realT>, dim>, PSI<realT> >, //thought that this is implicit...?
                    protected polyCore<dim,realT>
{

public:
  polynomials():polyCore<dim,realT>(){}
  polynomials(const polyParams<realT>& p_):polyCore<dim,realT>(p_){}
  polynomials(const polynomials& p2):polyCore<dim,realT>(p2){}

  virtual void setParameters(const PSI<realT>& It){
    polyCore<dim,realT>::setParameters(It);
  }

  virtual polynomials* clone() const{
    return new polynomials(*this);
  }

  const polyParams<realT>& getParameters() const{return polyCore<dim,realT>::getParameters();}

  void setParameters(polyParams<realT> p){polyCore<dim,realT>::setParameters(p);}

  Matrix<root_solver::complex<realT>, dim, 1> calcF(){return polyCore<dim,realT>::McalcF();}

  Matrix<root_solver::complex<realT>, dim, dim> calcJ(){return polyCore<dim,realT>::McalcJ();}

  void changePoint(Matrix<root_solver::complex<realT>, dim, 1> x_){polyCore<dim,realT>::MchangePoint(x_);}

  int getFunctionCallsCounter(){return polyCore<dim,realT>::getFunctionCallsCounter();}
  int getJacobianCallsCounter(){return polyCore<dim,realT>::getJacobianCallsCounter();}


  void randomiseTestcase(){polyCore<dim,realT>::randomiseTestcase();}

  virtual realT get_extra_parameter(){return polyCore<dim,realT>::get_extra_parameter();}
  virtual void set_extra_parameter(realT p){polyCore<dim,realT>::set_extra_parameter(p);}

  virtual realT get_initial(){return polyCore<dim,realT>::get_initial();}
  virtual realT get_final(){return polyCore<dim,realT>::get_final();}

  virtual realT get_max_change(){return polyCore<dim,realT>::get_max_change();}

};



template <typename realT>
class polynomial : public virtual root_solver::batch_functions<root_solver::extraSolver<root_solver::singleRootSolver<realT>, realT>, PSI<realT> >,
                   public virtual root_solver::batch_functions<root_solver::singleRootSolver<realT>, PSI<realT> >,
                   protected polyCore<1,realT>
{

public:
  polynomial():polyCore<1,realT>(){}
  polynomial(const polyParams<realT>& p_):polyCore<1,realT>(p_){}
  polynomial(const polynomial& p2):polyCore<1,realT>(p2){}

  virtual void setParameters(const PSI<realT>& It){
    polyCore<1,realT>::setParameters(It);
  }

  virtual polynomial* clone() const{
    return new polynomial(*this);
  }

  const polyParams<realT>& getParameters() const{return polyCore<1,realT>::getParameters();}

  void setParameters(polyParams<realT> p){polyCore<1,realT>::setParameters(p);}

  root_solver::complex<realT> calcF(){return polyCore<1,realT>::McalcF()(0);}

  root_solver::complex<realT> calcJ(){return polyCore<1,realT>::McalcJ()(0,0);}

  void changePoint(root_solver::complex<realT> x_){
    Matrix<root_solver::complex<realT>, 1, 1> x;
    x(0)=x_;
    polyCore<1,realT>::MchangePoint(x);
  }

  int getFunctionCallsCounter(){return polyCore<1,realT>::getFunctionCallsCounter();}
  int getJacobianCallsCounter(){return polyCore<1,realT>::getJacobianCallsCounter();}


  void randomiseTestcase(){polyCore<1,realT>::randomiseTestcase();}

  virtual realT get_extra_parameter(){return polyCore<1,realT>::get_extra_parameter();}
  virtual void set_extra_parameter(realT p){polyCore<1,realT>::set_extra_parameter(p);}

  virtual realT get_initial(){return polyCore<1,realT>::get_initial();}
  virtual realT get_final(){return polyCore<1,realT>::get_final();}

  virtual realT get_max_change(){return polyCore<1,realT>::get_max_change();}

};






#endif

