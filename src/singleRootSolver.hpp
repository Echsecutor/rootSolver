/**
 * @file singleRootSolver.hpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0.2014-06-11
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

#include <random>
#include <chrono>


using namespace std;

#include "rootSolver.hpp"
#include "complex.hpp"


namespace root_solver{

  template <typename valueT>
  class SRS_functions : public virtual functions<valueT,valueT>{
  public:
    virtual valueT guessStartPoint(){return valueT(0.0);}
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
  template <typename realT>
  class singleRootSolver : public virtual rootSolver<root_solver::complex<realT>,root_solver::complex<realT> >{
  public:
    //import typenames
    typedef rootSolver<root_solver::complex<realT>,root_solver::complex<realT> > root_solver_type;
    typedef typename root_solver_type::real_type real_type;
    typedef typename root_solver_type::scalar_type scalar_type;
    typedef typename root_solver_type::value_type value_type;

  protected:

    void update();

  public:

    typedef SRS_functions<scalar_type> functions_type;
    static const int value_dimension = 1;


    singleRootSolver(functions_type * f_) : root_solver_type(f_){}
    ~singleRootSolver(){}


    /// Setting a new starting point will reset the solver
    virtual void setStartPoint(value_type z_);

    using root_solver_type::step;

    virtual solver_state step();

  };


  //---------------------------------------------------------------------------

  template <typename realT>
  void singleRootSolver<realT>::update(){
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

  template <typename realT>
  void singleRootSolver<realT>::setStartPoint(value_type z_){
    this->z=z_;
    this->state=INIT;
    update();
  }


  template <typename realT>
  solver_state singleRootSolver<realT>::step(){

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


    if(this->J == 0){
      real_type c=1.0;
      default_random_engine* gen;
      normal_distribution<double> nd;
      int seed = chrono::system_clock::now().time_since_epoch().count();
      gen = new default_random_engine(seed);
      scalar_type actualZ = this->z;
      real_type dx = this->getMinStepSize();
      if(dx==0){
	dx= abs(this->z)*1e-2;
      }
      if(dx<1e-10){
	dx=1e-10;
      }
      while(this->J == 0.0){
#if DEBUG >= WARN
        cout << __FILE__ << " : I encountered a point with vanishing derivative. I will try to wiggle a little (scale = " << dx << ") away, but the solution found this way is likely very unstable."<<endl;
#endif
        this->z = actualZ + dx * nd(*gen);
        update();
        this->J = this->calc->calcJ();
        if(!(this->J == this->J)){//nan
          throw runtime_error(string(__FILE__) + string(" : derivative evaluation failed."));
        }
        if(c>100.0){
          throw runtime_error(string(__FILE__) + string(" : I can not find a point with non zero derivative."));
        }
        c++;
	dx*=2;
      }
    }

#if DEBUG >= SPAM
    cout << __FILE__ << " : Calculated J=" << (this->J)<<endl;
#endif

    value_type direction = (this->f) / (this->J);
    real_type stepSize = abs(direction);

    direction /= stepSize; // direction \in U(1)


    value_type newZ = (this->z) - stepSize * direction;
    this->calc->changePoint(newZ);
    value_type newF = this->calc->calcF();
    real_type newAbsF =abs(newF);

#if DEBUG>=ERR
    real_type checkError = abs(((newF - this->f)/(newZ - this->z) - this->J))/abs(this->J);
#if DEBUG>=SPAM
    cout << __FILE__ << " : relative deviation of finite differenze from derivative " << checkError << " for step size " << stepSize <<endl;
#endif
    if(checkError > 1e-2){
      cerr <<endl << __FILE__ << " : CAUTION! There is most likely an error in your derivatives at " << newZ <<endl;
    }
#endif

    solver_state newState = this->checkStatus(newAbsF, this->absF, stepSize);

    while(newState==REJECT){
#if DEBUG >= SPAM
      cout << __FILE__ << " : Backtracking at "<< newZ <<"\nwith f="<< newF <<"\ni.e. |f| = "<< newAbsF <<endl;
#endif
      stepSize /= 2.0;
      newZ = (this->z) - stepSize * direction;
      this->calc->changePoint(newZ);
      newF = this->calc->calcF();
      newAbsF = abs(newF);
      newState = this->checkStatus(newAbsF, this->absF, stepSize);
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
