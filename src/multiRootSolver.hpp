/**
 * @file multiRootSolver.hpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0.2014-05-21
 *
 * @section DESCRIPTION
 *
 * The multiRootSolver class is a solver for multi dimensional non-linear
 * equations using type casting, in particular it can easily be used for complex
 * functions and variables, other than
 * GSL http://www.gnu.org/software/gsl/manual/html_node/Multidimensional-Root_002dFinding.html
 * from which this project is completely independend. Instead we include
 * Eigen http://eigen.tuxfamily.org
 * for solving linear equation systems internally.
 *
 * Currently only Newtons mthod works properly.
 *
 * last changes: 
 * - uses the new mother class rootSolver to also include its sister singleRootSolver
 * into the family.
 *
 *
 * @section LICENSE
 *
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


#ifndef MULTIROOTSOLVER_H
#define MULTIROOTSOLVER_H

#include "preProDebugFlags.h"


#include <iostream>
#include <string>
#include <sstream>

using namespace std;

#include <eigen3/Eigen/Dense>
using namespace Eigen;

#include "rootSolver.hpp"

/** /brief This is an interface which shall be inherited by a user defined class implementing the relevant functions.
 *
 * - The solver will call changePoint explicitly, so that in between temporaryly computed values can be stored
 *
 * - Note that given that processors are much faster than memory this might not pay off, depending
 * on how expensive the evaluation of certain terms is.)
 *
 * - The solver will not ask for calculating F or J at the same point twice, so you will likely not want
 * to store the output.
 *
 * - Please throw an exception if the calculation fails.
 *
 */


template <typename T,  int dim>
  class MRS_functions : public virtual functions<Matrix<T, dim, 1>, Matrix<T, dim, dim>>{
 public:
  virtual Matrix<T, dim, 1> guessStartPoint(){
    Matrix<T, dim, 1> re;
    for (int i=0;i<dim;i++){
      re(i) = 0;
    }
    return re;
  };
};


/// MRS prefix for backward compatibility
static const solver_state MRS_NOINIT = solver_state::NOINIT;
static const solver_state MRS_STUCK = solver_state::STUCK;
static const solver_state MRS_SUCCESS = solver_state::SUCCESS;
static const solver_state MRS_INIT = solver_state::INIT;
static const solver_state MRS_CONTINUE = solver_state::CONTINUE;

typedef solver_state MRS_state;


enum MRS_mode
  {
    MRS_NEWTON,
    MRS_GRADIENT,
    MRS_DOGLEG
  };

/** \brief The multiRootSolver class is instanced for a given type of variables (float, double, complex<double>,...).
 *
 * The main steps for solving your equation are:
 *
 * 1. Derive a class from MRS_function and implement the function for which you search \f$ f(x)=0 \f$
 * as well as the jacobian \f$ J_{i,j} = \frac{\partial f_i}{\partial x_j} \f$ accordingly
 *
 * 2. Construct an instance multiRootSolver MRS(MRS_function new yourFunction());
 *
 * 3. Call MRS.setStartPoint(Matrix<T, dim, 1> z);
 *
 * 4. Run the solver like while (MRS.step(epsilonF,epsilonZ) == MRS_CONTINUE)
 *
 *
 * You can set a preferred mode of operation at construction
 * MRS(MRS_function new yourFunction(),MRS_mode preferredMode) or change it later MRS.setPrefferedMode(MRS_mode)
 * to be one of the following:
 *
 * - MRS_NEWTON (default) for primarily using the Newton-Raphson method with backtracking.
 * This converges quadratically if you are close to a root and backtracking improves global convergence. The down side is that
 * a linear system of equations has to be solved in each step, which might be slow in high dimension.
 *
 * - MRS_GRADIENT for primarily following the gradient flow of \f$ ||f||^2 \f$ which can be faster than Newton
 * if the calculation of the Jacobian matrix \f$ J \f$ is fast compared to solving the linear system \f$ Jx=f \f$.
 * Be aware that \f$ ||f||^2 \$f will likely have many local minima \f$ >0 \f$ where this method alone gets stuck.
 *
 * - MRS_DOGLEG searches the triangle with corners \f$\{ z, z-J^{-1} , z -\frac{||f||^2}{||J^T \bar{f}||^2} J^T \bar{f} \}\f$,
 * given by a Newton and a gradient flow step, for a point with significant lower value of \f$ ||f||^2 \f$.
 *
 * You can enforce using only the preffered mode by the corresponding flag at construction
 * MRS(MRS_function new yourFunction(),MRS_mode preferredMode=MRS_GRADIENT, bool onlyUsePrefferedMode=false) or later
 * MRS.setOnlyUsePrefferedMode(bool YN). Unless this flag is true, the solver will fallback to dogleg if one of the other two gets stuck.
 *
 */
template <typename T,  int dim>
  class multiRootSolver : public virtual rootSolver<Matrix<T, dim, 1>, Matrix<T, dim, dim>>{
 protected:

  static const int minImproveFactor = 8;///< a step is only counted as successfull if it decreases \f$ ||f|| \f$ by at least the precision goal devided by this factor
  static const int defaultEpsFByEpsZ = 16;///< default value for \f$ \frac{\epsilon_f}{\epsilon_z}\f$

  bool calculatedJ;
  double absLastStep;
  double lastAbsValueChange;

  void update();

  //variables controlling the actual solver
  MRS_mode currentMode;
  MRS_mode preferredMode;
  bool onlyUsePreferredMode;

  //temporary variables for the actual algorithms
  double epsF;
  double epsZ;
  Matrix<T, dim, 1> NewtonPoint;
  Matrix<T, dim, 1> NewtonF;
  void calcNewtonStep();
  Matrix<T, dim, 1> gradientPoint;
  Matrix<T, dim, 1> gradientF;
  void calcGradientStep();

  void simplexSearch();
  void backTrack(Matrix<T, dim, 1> farPoint, Matrix<T, dim, 1> farValue);


 public:

  typedef MRS_functions<T,dim> functions_type;
  typedef T numerical_type;
  static const int value_dimension = dim;


 multiRootSolver(functions_type * f_, MRS_mode preferredMode_=MRS_NEWTON, bool onlyUsePreferredMode_=true) : rootSolver<Matrix<T, dim, 1>, Matrix<T, dim, dim>>(f_),calculatedJ(false),absLastStep(0),lastAbsValueChange(0),currentMode(preferredMode_),preferredMode(preferredMode_),onlyUsePreferredMode(onlyUsePreferredMode_){}
  ~multiRootSolver(){}

  double getAbsLastStep(){return absLastStep;}
  double getLastAbsValueChange(){return lastAbsValueChange;}

  MRS_mode getPrefferedMode(){return preferredMode;}
  void setPrefferedMode(MRS_mode m){preferredMode=m;}

  /// Setting a new starting point will reset the solver
  void setStartPoint(Matrix<T, dim, 1> z_);

  /// This is the main part of the solver, actually computing a step towards the solution.
  /// Should be run while (step(epsilonF,epsilonZ) == MRS_CONTINUE)
  /// The first parameter sets a scale for f, defining the minimal improvement per step and maximal allowed final value for \f$ ||f||^2\f$.
  /// The first parameter similarly sets a scale for z, steps smaller than epsilonZ will be counted as stuck.
  MRS_state step(double epsilonF,double epsilonZ);
  MRS_state step(double epsilonF);

  string getLastPointInDatFormat();

};


//------------------------------------------------------------------------------------------------------


template <typename T, int dim>
  string multiRootSolver<T,dim>::getLastPointInDatFormat(){
  stringstream re;
  for (int i=0;i<dim;i++){
    if(i>0){re << "\t";}
    re << real(this->z(i)) << "\t" <<imag(this->z(i));
  }
  return re.str();
}


template <typename T,  int dim>
  void multiRootSolver<T,dim>::update(){
  this->calc->changePoint(this->z);
  this->f = this->calc->calcF();
  this->absF = this->f.norm();
#if DEBUG > 0
  cout << __FILE__ << " : Changed point to \n"<< this->z <<"\nwith f=\n"<< this->f <<"\ni.e. ||f|| = "<< this->absF <<endl;
#endif
}


template <typename T,  int dim>
  void multiRootSolver<T,dim>::setStartPoint(Matrix<T, dim, 1> z_){
  this->z=z_;
  this->calculatedJ=false;
  this->state=MRS_INIT;
  update();
}

/// tries to guess a good step size
template <typename T,  int dim>
  MRS_state multiRootSolver<T,dim>::step(double epsilonF){
  double epsilonZ= epsilonF / (double)defaultEpsFByEpsZ;
  if(this->state==MRS_CONTINUE && lastAbsValueChange>0){//should be equivalent
    double epsilonZ2 = absLastStep / lastAbsValueChange * this->absF;
    if(epsilonZ2<epsilonZ){
      epsilonZ=epsilonZ2;
    }
  }
  return step(epsilonF,epsilonZ);
}

template <typename T,  int dim>
  MRS_state multiRootSolver<T,dim>::step(double epsF_,double epsZ_){
  epsF=epsF_;
  epsZ=epsZ_;

  if(this->absF<epsF){
    this->state=MRS_SUCCESS;
#if DEBUG>=STATUS
    cout << __FILE__ << " : Current (start?) point is a root."<<endl;
#endif
  }
  if(this->state==MRS_SUCCESS && this->absF>epsF){//likely the user changed epsF on the fly
    this->state=MRS_CONTINUE;
  }

  if(this->state < MRS_INIT){
#if DEBUG>=WARN
    cout << __FILE__ << " : step() was called although the solver is not in a usable state. I won't do aything."<<endl;
#endif
    return this->state;
  }

  //------------------------------------
  // Choosing the solver algorithm:
  //------------------------------------
  switch(currentMode){
  case MRS_NEWTON:
    calcNewtonStep();
    backTrack(NewtonPoint, NewtonF);
    if(this->state==MRS_STUCK && !onlyUsePreferredMode){
      calcGradientStep();
      simplexSearch();
    }
    break;
  case MRS_GRADIENT:
    calcGradientStep();
    backTrack(gradientPoint, gradientF);
    if(this->state==MRS_STUCK && !onlyUsePreferredMode){
      calcNewtonStep();
      simplexSearch();
    }
    break;
  case MRS_DOGLEG:
    calcNewtonStep();
    calcGradientStep();
    simplexSearch();
    break;
  }
  //------------------------------------
#if DEBUG>=SPAM
  cout << __FILE__ << " : returning " << this->state <<endl;
#endif
  return this->state;
}





//------------------------------------
// Actual algorithms:
//------------------------------------

/// NewtonStep \f$ = -J^{-1}f \f$ using Eigen
template <typename T,  int dim>
  void multiRootSolver<T,dim>::calcNewtonStep(){
  if(!calculatedJ){
    this->J=this->calc->calcJ();
    calculatedJ=true;
  }
  NewtonPoint = this->z -1.0* this->J.colPivHouseholderQr().solve(this->f); //todo: optimize: only partial updates of J!
  this->calc->changePoint(NewtonPoint);
  NewtonF=this->calc->calcF();
#if DEBUG>DETAIL
  cout << __FILE__ << " : Newtons method proposed new point \n"<<NewtonPoint<<endl;
  cout << "with f(z)=\n"<<NewtonF<<endl;
#endif
}

/// gradientStep \f$ = -this->J^T \bar{f} \f$
template <typename T,  int dim>
  void multiRootSolver<T,dim>::calcGradientStep(){
  if(!calculatedJ){
    this->J=this->calc->calcJ();
    calculatedJ=true;
  }
  gradientPoint = this->J.transpose() * this->f.conjugate();
  double scale= this->absF*this->absF/gradientPoint.squaredNorm();
  gradientPoint = this->z - scale*gradientPoint; // aliasing ok as in i=a+b*i
  this->calc->changePoint(gradientPoint);
  gradientF=this->calc->calcF();
#if DEBUG>DETAIL
  cout << __FILE__ << " : Gradient method proposed new point \n"<<gradientPoint<<endl;
  cout << "with f(z)=\n"<<gradientF<<endl;
#endif

}

/// assuming that calc->changePoint(farPoint) was performed
template <typename T,  int dim>
  void multiRootSolver<T,dim>::backTrack( Matrix<T, dim, 1> farPoint, Matrix<T, dim, 1> farValue){
  double stepSize=(this->z-farPoint).norm();
#if DEBUG>=STATUS
  cout << __FILE__ << " : back tracking"<<endl;
#endif

  double farAbsF = farValue.norm();
  double change = this->absF - farAbsF;

  if(epsZ < stepSize/1024){
    epsZ=stepSize/2048;
#if DEBUG>=DETAIL
    cout << __FILE__ << " : epsZ was increased to " << epsZ <<" in order to prevent to deep recursion."<<endl;
#endif
  }

  if(change!=change){//nan
#if DEBUG>=ERR
    cerr << "ERR: MRS: Error in function evaluation!"<<endl;
#endif
    this->state=MRS_STUCK;
    return;
  }

#if DEBUG>=DETAIL
  cout << __FILE__ << " : step size " << stepSize << " leads to absF improvement: " <<change <<endl;
#endif

  if(farAbsF>epsF && stepSize < epsZ && change<epsF){
    this->state=MRS_STUCK;
#if DEBUG>=STATUS
    cout << __FILE__ << " : Step size smaller than epsZ = " << epsZ <<". Got Stuck."<<endl;
#endif
    return;//reject tiny steps
  }


  if((double)minImproveFactor*change > epsF || farAbsF<epsF){
#if DEBUG>=STATUS
    cout << __FILE__ << " : Back tracking successfull."<<endl;
#endif

    lastAbsValueChange=change;
    this->f=farValue;
    this->absF=farAbsF;
    absLastStep=stepSize;
    this->z=farPoint;
    calculatedJ=false;
    if(this->absF<epsF){
      this->state= MRS_SUCCESS;
    }else{
      this->state=MRS_CONTINUE;
    }
    return; //step accepted
  }

  //else: simple backtracking
  Matrix<T, dim, 1> newFarPoint=0.5*(farPoint+this->z);
  this->calc->changePoint(newFarPoint);

  if((farPoint-newFarPoint).norm()<epsZ){//this differenz might be zero due to rounding errors, i.e. machine precision, although (z-newFarPoint)>epsZ in next iteration... -.-
    this->state = MRS_STUCK;
#if DEBUG>=STATUS
    cout << __FILE__ << " : NewFarPoint to close to old one. state = STUCK = " << this->state <<endl;
#endif
    return;//reject tiny steps
  }

  backTrack(newFarPoint, this->calc->calcF());
}


void swap(void *A, void *B){
  void *C=A;
  A=B;
  B=C;
}

void swap3(void* A1,void* A2,void *A3,void * B1, void * B2, void * B3){
  swap(A1,B1);
  swap(A2,B2);
  swap(A3,B3);
}

template <typename T,  int dim>
  void multiRootSolver<T,dim>::simplexSearch(){
#if DEBUG>0
  cout << __FILE__ << " : Starting simplex search"<<endl;
#endif

  Matrix<T, dim, 1>  oldZ =this->z;

  //order triangle (for n=3, bubble sort is optimal):
  Matrix<T, dim, 1>  *minP =&this->z;
  Matrix<T, dim, 1>  *minF =&this->f;
  double min=this->absF;
  Matrix<T, dim, 1> * midP =&gradientPoint;
  Matrix<T, dim, 1> * midF =&gradientF;
  double mid=(*midF).norm();
  Matrix<T, dim, 1> * maxP =&NewtonPoint;
  Matrix<T, dim, 1> * maxF =&NewtonF;
  double max=(*maxF).norm();

  if(min!=min || mid!=mid || max!=max){//nan
    cerr << __FILE__ << " : Error in function evaluation!"<<endl;
    this->state=MRS_STUCK;
    return;
  }

  if(min>mid){
    swap3(minP,minF,&min,midP,midF,&mid);
  }
  if(mid>max){
    swap3(midP,midF,&mid,maxP,maxF,&max);
  }
  if(min>mid){
    swap3(minP,minF,&min,midP,midF,&mid);
  }

  Matrix<T, dim, 1> centerP;
  Matrix<T, dim, 1> centerF;
  double center;

  double size;
  double change;

  while (true){//actually this loop is always terminated by return
    centerP=(*minP+*midP+*maxP)/3.0;
    this->calc->changePoint(centerP);
    centerF=this->calc->calcF();
    center = centerF.norm();
    size = (*maxP - centerP).norm();
    change = this->absF - center;
#if DEBUG>0
    cout << __FILE__ << " : New center point \n"<<centerP<<endl;
    cout << __FILE__ << " : New center value \n"<<centerF<<endl;
#endif

    if(size < epsZ){
      this->state=MRS_STUCK;
#if DEBUG>0
      cout << __FILE__ << " : Step size " << size << " smaller than epsZ = " << epsZ <<". Got Stuck."<<endl;
#endif
      return;//reject tiny steps
    }

    if((double)minImproveFactor*change > epsF || center<epsF){
#if DEBUG>0
      cout << __FILE__ << " : Triangle search successfull."<<endl;
#endif
      lastAbsValueChange=change;
      this->f=centerF;
      this->absF=center;
      absLastStep=(oldZ-centerP).norm();
      this->z=centerP;
      calculatedJ=false;
      if(center<epsF){
        this->state= MRS_SUCCESS;
      }else{
        this->state=MRS_CONTINUE;
      }
      return; //step accepted
    }

#if DEBUG>0
    cout << __FILE__ << " : Replacing \n"<<*maxP<<"\n("<<max<<") with\n"<<centerP<<"\n("<<center<<")"<<endl;
#endif

    //update triangle, replace max:
    maxP=&centerP;
    maxF=&centerF;
    max=center;

    if(mid>max){
      swap3(midP,midF,&mid,maxP,maxF,&max);
    }
    if(min>mid){
      swap3(minP,minF,&min,midP,midF,&mid);
    }


  }//wend center

}



//------------------------------------


#endif
