/**
 * @file singleRootSolver.hpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0-2014.05.21
 *
 * @section DESCRIPTION
 *
 * This is a simple implementation of Newtons method for one complex variable/equation using back traking.
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License and a copy of the GNU General Public Lic
 ense along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 *
 */


#ifndef SINGLEROOTSOLVER_HPP
#define SINGLEROOTSOLVER_HPP


#include "preProDebugFlags.h"


#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>


using namespace std;

#include "rootSolver.hpp"

namespace root_solver{

  template <typename T>
  class SRS_functions : public virtual functions<T,T >{
  public:
    virtual T guessStartPoint(){return 0.0;}
  };



  /** The singleRootSolver class is instanced for a given type T of variables (complex<double>, complex<...>, ...).
   * Notice that for one dimensional root finding much more efficient braketing algorithms exist,
   * hence you should rather use e.g. the gsl routines instead of singleRootSolver<double> or so.
   *
   * The main steps for solving your equation are:
   *
   * 1. Derive a class from SRS_function and implement the function for which you search \f$ f(z)=0 \f$
   * as well as the derivative \f$ J = \frac{\partial f}{\partial z} \f$ accordingly.
   *
   * 2. Construct an instance singleRootSolver SRS(SRS_function new yourFunction());
   *
   * 3. Set appropriate thresh-holds to determine success/failure of the solver with
   * SRS.setPrecisionGoal(goal);
   * and at least one of the following:
   * SRS.setMinStepSize(stepSize);
   * SRS.setMinImprovePerStep(impro);
   * SRS.setMinRelativeImprovePerStep(impro);
   *
   * 4. Call SRS.setStartPoint(T z);
   *
   * 4. Run the solver like while (SRS.step() >= CONTINUE)
   *
   */
  template <typename T>
  class singleRootSolver : public virtual rootSolver<T,T>{
  protected:

    double epsZ;
    double epsF;

    void update();

  public:

    typedef SRS_functions<T> functions_type;
    typedef T numerical_type;
    static const int value_dimension = 1;


    singleRootSolver(functions_type * f_) : rootSolver<T,T >(f_){}
    ~singleRootSolver(){}


    /// Setting a new starting point will reset the solver
    virtual void setStartPoint(T z_);

    virtual solver_state step();

    //no clue why this is not automatic.....
    virtual solver_state step(double epsF){return rootSolver<T,T>::step(epsF);}
    virtual solver_state step(double epsF, double epsZ){return rootSolver<T,T>::step(epsF,epsZ);}


  };


  //---------------------------------------------------------------------------

  template <typename T>
  void singleRootSolver<T>::update(){
    this->calc->changePoint(this->z);
    this->f = this->calc->calcF();
    this->absF = abs((this->f));
#if DEBUG >= SPAM
    cout << __FILE__ << " : Changed point to "<< (this->z) <<"\nwith f="<< (this->f) <<"\ni.e. |f| = "<< this->absF <<endl;
#endif
    if(!(this->absF == this->absF)){//nan
      throw runtime_error(string(__FILE__) + string(" : function evaluation failed."));
    }
  }

  template <typename T>
  void singleRootSolver<T>::setStartPoint(T z_){
    this->z=z_;
    this->state=INIT;
    update();
  }


  template <typename T>
  solver_state singleRootSolver<T>::step(){

    if(this->absF < this->getPrecisionGoal()){
      this->state=SUCCESS;
#if DEBUG >= DETAIL
      cout << __FILE__ << " : Current (start?) point is a root."<<endl;
#endif
      return this->state;
    }else if(this->state == SUCCESS){//likely the user changed epsF on the fly
#if DEBUG >= DETAIL
      cout << __FILE__ << " : |f| to large, solver will continue."<<endl;
#endif
      this->state=CONTINUE;
    }

    if(this->state < INIT){
#if DEBUG >= WARN
      cout << __FILE__ << " : step() was called although the solver is not in a usable state. I won't do aything."<<endl;
#endif
      return this->state;
    }

    this->state = CONTINUE;

    this->J = this->calc->calcJ();
    if(!(this->J == this->J)){//nan
      throw runtime_error(string(__FILE__) + string(" : derivative evaluation failed."));
    }

#if DEBUG >= SPAM
    cout << __FILE__ << " : Calculated J=" << (this->J)<<endl;
#endif

    T direction = (this->f) / (this->J);
    double stepSize=abs(direction);

    direction /= stepSize; // direction \in U(1)


    T newZ = (this->z) - stepSize * direction;
    this->calc->changePoint(newZ);
    T newF = this->calc->calcF();
    double newAbsF =abs(newF);

    solver_state newState = checkStatus(newAbsF, this->absF, stepSize);

    while(newState==REJECT){
#if DEBUG >= SPAM
      cout << __FILE__ << " : Backtracking at "<< newZ <<"\nwith f="<< newF <<"\ni.e. |f| = "<< newAbsF <<endl;
#endif
      stepSize /= 2.0;
      newZ = (this->z) - stepSize * direction;
      this->calc->changePoint(newZ);
      newF = this->calc->calcF();
      newAbsF = abs(newF);
      newState = checkStatus(newAbsF, this->absF, stepSize);
    }

    this->f = newF;
    this->z = newZ;
    this->absF=newAbsF;

#if DEBUG >= DETAIL
    cout << __FILE__ << " : Changed point to "<< (this->z) <<"\nwith f="<< (this->f) <<"\ni.e. |f| = "<< this->absF <<endl;
#endif

    this->state=newState;
    return this->state;

  }



}//end namespace

#endif
