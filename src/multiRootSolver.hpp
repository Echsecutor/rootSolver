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
#include <stdexcept>
using namespace std;

#include "eigenWrapper.hpp"
using namespace Eigen;

#include "rootSolver.hpp"
#include "complex.hpp"

namespace root_solver{


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


  template <typename scalarT,  int dim>
  class MRS_functions : public virtual functions<Matrix<scalarT, dim, 1>, Matrix<scalarT, dim, dim> >{
  public:
    typedef functions<Matrix<scalarT, dim, 1>, Matrix<scalarT, dim, dim> > functions_type;

    virtual typename functions_type::value_type guessStartPoint(){
      typename functions_type::value_type re;
      for (int i=0;i<dim;i++){
        re(i) = 0;
      }
      return re;
    };
  };


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
   * 3. Set appropriate thresh-holds to determine success/failure of the solver with
   * MRS.setPrecisionGoal(goal);
   * and at least one of the following:
   * MRS.setMinStepSize(stepSize);
   * MRS.setMinImprovePerStep(impro);
   * MRS.setMinRelativeImprovePerStep(impro);
   *
   * 4. Call MRS.setStartPoint(Matrix<T, dim, 1> z);
   *
   * 5. Run the solver like while (MRS.step() >= CONTINUE)
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
   *
   * If you want to use this for real numbers (better implementations are available and) you need to specify scalarT=realT
   */
  template <typename scalarT, int dim>
  class multiRootSolver : public virtual rootSolver<Matrix<scalarT, dim, 1>, Matrix<scalarT, dim, dim> >{
  public:
    //import typenames
    typedef rootSolver<Matrix<scalarT, dim, 1>, Matrix<scalarT, dim, dim> > root_solver_type;
    typedef typename root_solver_type::real_type real_type;
    typedef typename root_solver_type::scalar_type scalar_type;
    typedef typename root_solver_type::value_type value_type;

  protected:
    
    bool calculatedJ;
    real_type absLastStep;
    real_type lastAbsValueChange;

    void update();

    //variables controlling the actual solver
    MRS_mode currentMode;
    MRS_mode preferredMode;
    bool onlyUsePreferredMode;

    //temporary variables for the actual algorithms
    value_type NewtonPoint;
    value_type NewtonF;
    void calcNewtonStep();

    /*
      Matrix<T, dim, 1> gradientPoint;
      Matrix<T, dim, 1> gradientF;
      void calcGradientStep();

      void simplexSearch();
    */

    void backTrack(value_type farPoint,value_type farValue);


  public:

    typedef MRS_functions<scalar_type,dim> functions_type;
    static const int value_dimension = dim;


    multiRootSolver(functions_type * f_, MRS_mode preferredMode_=MRS_NEWTON, bool onlyUsePreferredMode_=true) : root_solver_type(f_),calculatedJ(false),absLastStep(0),lastAbsValueChange(0),currentMode(preferredMode_),preferredMode(preferredMode_),onlyUsePreferredMode(onlyUsePreferredMode_){}
    ~multiRootSolver(){}

    real_type getAbsLastStep(){return absLastStep;}
    real_type getLastAbsValueChange(){return lastAbsValueChange;}

    MRS_mode getPrefferedMode(){return preferredMode;}
    void setPrefferedMode(MRS_mode m){preferredMode=m;}


    /// Setting a new starting point will reset the solver
    void setStartPoint(value_type z_);

    using root_solver_type::step;

    /// This is the main part of the solver, actually computing a step towards the solution.
    /// Should be run while (step() >= CONTINUE)
    solver_state step();

    //no clue why this is not automatic.....
    // virtual solver_state step(double epsF){return rootSolver<Matrix<T, dim, 1>, Matrix<T, dim, dim> >::step(epsF);}
    // virtual solver_state step(double epsF, double epsZ){return rootSolver<Matrix<T, dim, 1>, Matrix<T, dim, dim> >::step(epsF,epsZ);}


  };


  //------------------------------------------------------------------------------------------------------


  template <typename scalarT,  int dim>
  void multiRootSolver<scalarT,dim>::update(){
    this->calc->changePoint(this->z);
    this->f = this->calc->calcF();
    this->absF = this->f.norm();
#if DEBUG > 0
    cout << __FILE__ << " : Changed point to \n"<< this->z <<"\nwith f=\n"<< this->f <<"\ni.e. ||f|| = "<< this->absF <<endl;
#endif
    if(!(this->absF == this->absF)){//nan
      throw runtime_error(string(__FILE__) + string(" : function evaluation failed."));
    }
  }


  template <typename scalarT,  int dim>
  void multiRootSolver<scalarT,dim>::setStartPoint(value_type z_){
    this->z=z_;
    this->calculatedJ=false;
    this->state=INIT;
    update();
  }

  template <typename scalarT,  int dim>
  solver_state multiRootSolver<scalarT,dim>::step(){

    if(this->absF < this->getPrecisionGoal()){
      this->state=SUCCESS;
#if DEBUG>=STATUS
      cout << __FILE__ << " : Current (start?) point is a root."<<endl;
#endif
    }

    if(this->state < INIT){
#if DEBUG>=WARN
      cout << __FILE__ << " : step() was called although the solver is not in a usable state. I won't do aything."<<endl;
#endif
      return this->state;
    }

    this->state=CONTINUE;

    //------------------------------------
    // Choosing the solver algorithm:
    //------------------------------------
    switch(currentMode){
    case MRS_NEWTON:
      calcNewtonStep();
      backTrack(NewtonPoint, NewtonF);
      /*
        if(this->state==STUCK && !onlyUsePreferredMode){
        calcGradientStep();
        simplexSearch();
        }
      */
      break;

      /*
        case MRS_GRADIENT:
        calcGradientStep();
        backTrack(gradientPoint, gradientF);
        if(this->state==STUCK && !onlyUsePreferredMode){
        calcNewtonStep();
        simplexSearch();
        }
        break;
        case MRS_DOGLEG:
        calcNewtonStep();
        calcGradientStep();
        simplexSearch();
        break;
      */
    default:
      throw runtime_error(string(__FILE__) + string(" : Currently only Newton's method is implemented!"));
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
  template <typename scalarT,  int dim>
  void multiRootSolver<scalarT,dim>::calcNewtonStep(){
    if(!calculatedJ){
      this->J=this->calc->calcJ();
      calculatedJ=true;
    }
    NewtonPoint = this->z - this->J.colPivHouseholderQr().solve(this->f); //todo: optimize: only partial updates of J!
    this->calc->changePoint(NewtonPoint);
    NewtonF=this->calc->calcF();
#if DEBUG>DETAIL
    cout << __FILE__ << " : Newtons method proposed new point \n"<<NewtonPoint<<endl;
    cout << "with f(z)=\n"<<NewtonF<<endl;
#endif
  }

  /// assuming that calc->changePoint(farPoint) was performed
  template <typename scalarT,  int dim>
  void multiRootSolver<scalarT,dim>::backTrack(value_type farPoint, value_type farValue){
    real_type stepSize=(this->z-farPoint).norm();
#if DEBUG>=STATUS
    cout << __FILE__ << " : back tracking"<<endl;
#endif

    real_type farAbsF = farValue.norm();
    real_type change = this->absF - farAbsF;

    if(change!=change){//nan
      throw runtime_error(string(__FILE__) + string(" :  Error in function evaluation!"));
    }

#if DEBUG>=DETAIL
    cout << __FILE__ << " : step size " << stepSize << " leads to absF improvement: " <<change <<endl;
#endif

    this->state = this->checkStatus(farAbsF, this->absF, stepSize);

    // simple backtracking
    if(this->state==REJECT){
      value_type newFarPoint = 0.5 * (farPoint + this->z);
      this->calc->changePoint(newFarPoint);
      backTrack(newFarPoint, this->calc->calcF());
    }else{
      lastAbsValueChange=change;
      this->f=farValue;
      this->absF=farAbsF;
      absLastStep=stepSize;
      this->z=farPoint;
      calculatedJ=false;
    }
  }


  /*
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
  */

  //todo:rewrite for new checkStatus() etc
  /*
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
    this->state=STUCK;
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
    this->state=STUCK;
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
    this->state= SUCCESS;
    }else{
    this->state=CONTINUE;
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
  */


  /*
 /// gradientStep \f$ = -this->J^T \bar{f} \f$
 template <typename T,  int dim>
 void multiRootSolver<T,dim>::calcGradientStep(){
 if(!calculatedJ){
 this->J=this->calc->calcJ();
 calculatedJ=true;
 }
 gradientPoint = this->J.transpose() * this->f.conjugate();
 double scale= this->absF*this->absF/gradientPoint.squaredNorm();
 gradientPoint = this->z - scale*gradientPoint;
 this->calc->changePoint(gradientPoint);
 gradientF=this->calc->calcF();
 #if DEBUG>DETAIL
 cout << __FILE__ << " : Gradient method proposed new point \n"<<gradientPoint<<endl;
 cout << "with f(z)=\n"<<gradientF<<endl;
 #endif

 }
  */


  //------------------------------------


}//end namespace

#endif
