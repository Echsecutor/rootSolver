/**
 * @file extrapolationSolver.hpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0.2014-05-15
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


/// if you want complex random shifts on top of the extrapolation, use something like
/// #ifndef COMPLEX_I
/// #define COMPLEX_I
/// static complex<double> I(0.0,1.0);
/// #endif
/// BEFORE including this header

#ifndef COMPLEX_I
static double I(1.0);
#endif

//------------------------------------------------------------------------------------------
// Declarations:
//------------------------------------------------------------------------------------------

/// extra_functions interface
/// only one dimensional extrapolation is (currently) implemented,
/// hence the type of \f$p\f$ is fixed to double
template <typename rootSolverT>
class extra_functions : public virtual rootSolverT::functions_type{
public:

  virtual double get_extra_parameter()=0;
  virtual void set_extra_parameter(double p)=0;

  virtual double get_initial()=0;
  virtual double get_final()=0;

  virtual double get_max_change()=0;


  virtual double get_direction(){ //< unlikely to be overwritten
    if (get_final() < get_initial())
      return -1.0;
    return 1.0;
  }

  // overwrite these if you want something more sophisticated. ;)
  virtual double get_min_change(){
    return get_max_change() / pow(2,20);
  }

  virtual double get_initial_change(){
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
template <typename rootSolverT>
class extraSolver : public virtual rootSolverT{

protected:

  default_random_engine* gen;

  extra_functions<rootSolverT>* calc;

  double dP;

  bool randomiseExtrapolation;

  typedef extrapolationData<typename rootSolverT::value_type, typename rootSolverT::numerical_type, rootSolverT::value_dimension> data_type;

  data_type dat;

  void speedUp(const double& steps);
  bool slowDown(double factor=1.5);


public:

  typedef extra_functions<rootSolverT> functions_type;

  unsigned int usePoints;///< maximal number of points to use for extrapolation
  unsigned int desiredSteps;///< per internal solver run


  //explicit call to rootSolver constructor needed?
  extraSolver(functions_type * f_, bool randomise = true, int seed = 42, unsigned int desiredNumberOfSteps=50, unsigned int usePoints_=8) : rootSolver<typename rootSolverT::value_type, typename rootSolverT::derivative_type>(f_), rootSolverT(f_), gen(new default_random_engine (seed)), calc(f_), randomiseExtrapolation(randomise), usePoints(usePoints_), desiredSteps(desiredNumberOfSteps){
    dP = f_->get_initial_change();
  }


  virtual solver_state step(double epsilonF, double epsilonZ);
  virtual solver_state step(double epsilonF){return rootSolverT::step(epsilonF);};//no clue why this is not automatic....

};



//------------------------------------------------------------------------------------------
// Implementation:
//------------------------------------------------------------------------------------------

/// One "step" of the extrapolation solver actually consists of a full
/// run of the underlying solver to produce a new data point for the
/// next extrapolation.
template <typename rootSolverT>
solver_state extraSolver<rootSolverT>::step(double epsilonF,double epsilonZ){



  if(calc->get_extra_parameter()==calc->get_initial() || dat.size()==0){//should be equivalent
#if DEBUG>=SPAM
    cout << __FILE__ << " : Guessing start point." <<endl;
#endif
    rootSolverT::setStartPoint(calc->guessStartPoint());// for p_initial a good starting point should be known
  }else{

    bool goodPoint=false;

    while(!goodPoint){
#if DEBUG>=DETAIL
      cout << __FILE__ << " : Extrapolating start point." <<endl;
#endif
      if(this->randomiseExtrapolation){
        rootSolverT::setStartPoint(dat.extrapolate(calc->get_extra_parameter()) + dP*1e-2 * (dat.gaussBlob(gen)+ I*dat.gaussBlob(gen)) );
      }else{
        rootSolverT::setStartPoint(dat.extrapolate(calc->get_extra_parameter()));
      }

      if(rootSolverT::getAbsF() > 1e-2){ //todo:hardcoded...
        //we are likely running to fast
#if DEBUG>=DETAIL
        cout << __FILE__ << " : Bad starting point. " <<endl;
#endif
        if(!slowDown(1.1)){
          goodPoint=true;//at least give this one a try
        }
      }else{
        goodPoint=true;
      }
    }

  }



#if DEBUG>=DETAIL
  cout << __FILE__ << " : Start solving for parameter " << calc->get_extra_parameter() <<endl;
#endif

  int step=1;
  while (rootSolverT::step(epsilonF,epsilonZ) == CONTINUE){
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

    if(calc->get_direction() * calc->get_extra_parameter() > calc->get_direction() * calc->get_final())
      calc->set_extra_parameter(calc->get_final());

    while(dat.size() > usePoints){
      dat.pop_front();
    }

  }else{
#if DEBUG>=DETAIL
    cout << __FILE__ << " : The solver got stuck. Reducing step size and retrying."<<endl;
#endif
    if(!slowDown()){
#if DEBUG>=WARN
      cout << __FILE__ << " : Step size to small. Extrapolation got stuck."<<endl;
#endif
      this->state=STUCK;
      return this->state;
    }

  }

  this->state=CONTINUE;
  return this->state;

}


template <typename rootSolverT>
void extraSolver<rootSolverT>::speedUp(const double& steps){
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
  }

  calc->set_extra_parameter(calc->get_extra_parameter() + dP);
}

template <typename rootSolverT>
bool extraSolver<rootSolverT>::slowDown(double factor){
  double oldDp=dP;

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

  calc->set_extra_parameter(calc->get_extra_parameter() - oldDp + dP);

  return true;
}



#endif
