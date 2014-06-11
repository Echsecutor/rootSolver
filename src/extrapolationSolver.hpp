/**
 * @file extrapolationSolver.hpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0.2014-06-11
 *
 * @section Version number format
 *
 * The Version number is formatted as "M.S.D" where M is the major
 * release branch (backward compatibility to all non-alpha releases of
 * the same branch is guaranteed), S is the state of this release (0
 * for alpha, 1 for beta, 2 for stable), and D is the date formatted
 * as yyyy-mm-dd.)
 *
 *
 * @section DESCRIPTION
 *
 * This is a higher level interface for the rootSolver for the following scenario:
 * We are looking for solutions of the (set of) equations \f$ f(z)=0\f$ for some
 * \f$ f_p \mathbb{C}^d \to \mathbb{C}^d \f$ where $p$ denotes a real parameter.
 * For some initial value \f$ f_{p_i}(z_i)=0\f$ an (approximate) solution is known. The solution
 * for some other final value \f$ f_{p_f}(z_f)=0\f$ is desired and it is assumed that there is a continuous
 * function \f$ z(p)\f$ such that \f$ f_{p}(z(p))=0\f$ for all \f$p\in[p_i,p_f]\f$ with \f$z(p_i)=z_i \f$.
 *
 * The purpose of this wrapper is to slowly tune \f$p\f$ from the initial to the target value,
 * solving the equations along the way. The speed of this tuning is adapted such as to reach the target
 * quickly while staying within the bassin of attraction of the right zero in each step.
 * The new starting value for the solver will be extrapolated from the known ones.
 *
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


#ifndef EXTRAPOLATIONSOLVER_HPP
#define EXTRAPOLATIONSOLVER_HPP

#include "preProDebugFlags.h"

#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "rootSolver.hpp"
#include "extrapolationData.hpp"


//------------------------------------------------------------------------------------------
// Declarations:
//------------------------------------------------------------------------------------------

namespace root_solver{

  // static const solver_state EXTRAPOLATION_SUCCESS(10);
  // static const solver_state EXTRAPOLATION_STUCK(-11);


  /// extra_functions interface
  /// only one dimensional extrapolation is (currently) implemented,
  /// hence the type of \f$p\f$ is fixed to double
  template <typename rootSolverT, typename parameterT>
  class extra_functions : public virtual rootSolverT::functions_type{
  public:

    virtual parameterT get_extra_parameter()=0;
    virtual void set_extra_parameter(parameterT p)=0;

    virtual parameterT get_initial()=0;
    virtual parameterT get_final()=0;

    virtual parameterT get_max_change()=0;


    virtual parameterT get_direction(){ //< unlikely to be overwritten
      if (get_final() < get_initial())
        return -1.0;
      return 1.0;
    }

    // overwrite these if you want something more sophisticated. ;)
    virtual parameterT get_min_change(){
      return get_max_change() / pow(2,20);
    }

    virtual parameterT get_initial_change(){
      return get_min_change() * pow(2,10);
    }


  };


  /** /brief The extraSolver is a wrapper class.
   *
   * The actual equation solving is done by rootSolverT internally,
   * hence an appropriate implementation of the rootSolver
   * (e.g. multiRootSolver or singleRootSolver) has to be specified.
   *
   * numericalT should be the same as the numerical type of the underlying rootSolver (e.g. complex<double>) and dataDim should be the dimension of value_type. so extraSolver<multiRootSolver<complex<double >, dim >, dim>
   *
   */
  template <typename rootSolverT, typename parameterT>
  class extraSolver : public virtual rootSolverT{
  public:
    //import typenames
    typedef rootSolverT root_solver_type;
    typedef typename root_solver_type::real_type real_type;
    typedef typename root_solver_type::scalar_type scalar_type;
    typedef typename root_solver_type::value_type value_type;
    typedef typename root_solver_type::derivative_type derivative_type;

  protected:

    default_random_engine* gen;

    extra_functions<rootSolverT, parameterT>* calc;

    parameterT dP;

    bool randomiseExtrapolation;

    typedef extrapolationData<value_type, scalar_type, real_type, root_solver_type::value_dimension> data_type;

    data_type dat;

    void speedUp(const unsigned int steps);
    bool slowDown(parameterT factor=1.5);


  public:

    typedef extra_functions<rootSolverT, parameterT> functions_type;

    unsigned int usePoints;///< maximal number of points to use for extrapolation
    unsigned int desiredSteps;///< per internal solver run


    //explicit call to rootSolver constructor needed?
    extraSolver(functions_type * f_, bool randomise = true, int seed = 42, unsigned int desiredNumberOfSteps=50, unsigned int usePoints_=8) : rootSolver<value_type, derivative_type>(f_), rootSolverT(f_), gen(new default_random_engine (seed)), calc(f_), randomiseExtrapolation(randomise), usePoints(usePoints_), desiredSteps(desiredNumberOfSteps){
      dP = f_->get_initial_change();
    }


    virtual solver_state step();

  };



  //------------------------------------------------------------------------------------------
  // Implementation:
  //------------------------------------------------------------------------------------------


  /// One "step" of the extrapolation solver actually consists of a full
  /// run of the underlying solver to produce a new data point for the
  /// next extrapolation.
  template <typename rootSolverT, typename parameterT>
  solver_state extraSolver<rootSolverT,parameterT>::step(){

    if(dat.size()==0){
#if DEBUG>=SPAM
      cout << __FILE__ << " : Guessing start point." <<endl;
#endif
      rootSolverT::setStartPoint(calc->guessStartPoint());// for p_initial a good starting point should be known
#if DEBUG>=WARN
      if(rootSolverT::getAbsF() > this->getPrecisionGoal()){
        cout << __FILE__ << " : Bad starting point guess! |f|=" << rootSolverT::getAbsF()  <<endl;
      }
#endif

    }else{

    calc->set_extra_parameter(calc->get_extra_parameter() + dP);
    if(calc->get_direction() * calc->get_extra_parameter() > calc->get_direction() * calc->get_final())
      calc->set_extra_parameter(calc->get_final());

#if DEBUG>=SPAM
    cout << __FILE__ << " : Changed to P = " << calc->get_extra_parameter() <<endl;
#endif

      bool goodPoint=false;
      real_type randomness=1e-12;
      int tries=0;
      while(!goodPoint){
        tries++;
#if DEBUG>=DETAIL
        cout << __FILE__ << " : Extrapolating start point." <<endl;
#endif
        if(this->randomiseExtrapolation){
          value_type ex = dat.extrapolate(calc->get_extra_parameter(), &randomness);
	  randomness += 1e-12;
	  scalar_type I(0,1.0);
          rootSolverT::setStartPoint(ex + randomness * (dat.gaussBlob(gen) + I * dat.gaussBlob(gen)));
        }else{
          rootSolverT::setStartPoint(dat.extrapolate(calc->get_extra_parameter()));
        }

        if(rootSolverT::getAbsF() > 1e-2){ //todo:hardcoded...
#if DEBUG>=DETAIL
          cout << __FILE__ << " : Bad starting point. " <<endl;
#endif
          //we are likely running to fast
          if(!slowDown(1.3)){
            goodPoint=true;//at least give this one a try
          }
          //or extrapolation is bad -> fall back towards constant extrapolation
          if(tries > 3 && dat.size() > dat.degree){
            dat.pop_front();
          }
	  calc->set_extra_parameter(calc->get_extra_parameter() + dP);
#if DEBUG>=SPAM
	  cout << __FILE__ << " : Changed to P = " << calc->get_extra_parameter() <<endl;
#endif
        }else{
          goodPoint=true;
        }

      }
    }



#if DEBUG>=DETAIL
    cout << __FILE__ << " : Start solving for parameter " << calc->get_extra_parameter() <<endl;
#endif

    unsigned int step=1;
    while (rootSolverT::step() >= CONTINUE){
#if DEBUG>=SPAM
      cout << __FILE__ << " : after " << step << " steps, the solver achieved |f| = " << rootSolverT::getAbsF() <<endl;
#endif
      step++;
    }

    if (rootSolverT::getState()==SUCCESS){
#if DEBUG>=DETAIL
      cout << __FILE__ << " : The solver found an approximate root for parameter " << calc->get_extra_parameter() << " at " << rootSolverT::getLastPoint() << " with |f|= " << rootSolverT::getAbsF() << endl;
#endif

      if (calc->get_extra_parameter() == calc->get_final()){
        this->state=SUCCESS;
        return this->state;
      }

      dat.push_back(rootSolverT::getLastPoint(), calc->get_extra_parameter());

#if DEBUG>=DETAIL
      cout << __FILE__ << " : Saving root for extrapolation. Now we have "<< dat.size() << " in store." << endl;
#endif

      speedUp(step);

      while(dat.size() > usePoints){
        dat.pop_front();
      }

    }else{
#if DEBUG>=DETAIL
      cout << __FILE__ << " : The solver got stuck. Reducing step size and retrying."<<endl;
#endif

      if(dat.size()==0){
#if DEBUG>=WARN
        cout << __FILE__ << " : Could not find an initial root. Check your starting point / reduce the precision goal!"<<endl;
#endif
        this->state=STUCK;
        return this->state;
      }

      if(!slowDown()){
#if DEBUG>=WARN
        cout << __FILE__ << " : Step size to small. Extrapolation got stuck."<<endl;
#endif
	if(this->getPrecisionGoal() > 1e-5){//todo:hardcoded
	  this->state=STUCK;
	  return this->state;
	}

	this->setPrecisionGoal(this->getPrecisionGoal()*10.0);

#if DEBUG>=WARN
        cout << __FILE__ << " : Reducing precision Goal to " << this->getPrecisionGoal() << " at parameter " << calc->get_extra_parameter() <<endl;
#endif
	speedUp(0);
      }


    }

    this->state=CONTINUE;
    return this->state;

  }


  template <typename rootSolverT, typename parameterT>
  void extraSolver<rootSolverT,parameterT>::speedUp(const unsigned int steps){
    if(steps < desiredSteps/2){
      dP*=2.0;
      if(dP>calc->get_max_change()){
        dP=calc->get_max_change();
      }
#if DEBUG>=SPAM
      cout << __FILE__ << " : Speed up to dP = " << dP <<endl;
#endif
    }else if(steps > desiredSteps){
      dP/=2.0;
#if DEBUG>=SPAM
      cout << __FILE__ << " : Slowdown to dP = " << dP <<endl;
#endif
    }else{
      desiredSteps++;//avoids wasting time with tiny dPs
    }
  }

  template <typename rootSolverT, typename parameterT>
  bool extraSolver<rootSolverT,parameterT>::slowDown(parameterT factor){
    parameterT oldDp=dP;

    dP /= factor;
#if DEBUG>=SPAM
    cout << __FILE__ << " : Slowdown to dP = " << dP << " ( minDP = " << calc->get_min_change() << " )" << endl;
#endif
    if(dP < calc->get_min_change()){
#if DEBUG>=DETAIL
      cout << __FILE__ << " : Step size to small!" <<endl;
#endif
      return false;
    }

    calc->set_extra_parameter(calc->get_extra_parameter() - oldDp);

#if DEBUG>=SPAM
    cout << __FILE__ << " : slowDown() reset to P = " << calc->get_extra_parameter() <<endl;
#endif

    return true;
  }

}//end namespace

#endif
