/**
 * @file rootSolver.hpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0-2014-06-11
 *
 * @section DESCRIPTION
 *
 * This file hosts the base class rootSolver from which multiRootSolver and singleRoot solver are derived.
 * MultirootSolver is a solver for multi dimensional and singleRoot Solver for one dimensional (non-linear)
 * equations using type casting, in particular it can easily be used for complex
 * functions and variables, other than
 * GSL http://www.gnu.org/software/gsl/manual/html_node/Multidimensional-Root_002dFinding.html
 * from which this project is completely independend.
 *
 * The functions interface is derived in a single and multi version to interfaces which are
 * ultimately implemented by the user to provide the set of equations to be solved.
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


#ifndef ROOTSOLVER_HPP
#define ROOTSOLVER_HPP

#include "preProDebugFlags.h"

#include <iostream>
#include <string>
#include<limits>

//for traits:
#include <type_traits>
#include <eigen3/Eigen/Core>


using namespace std;


//todo: do we need to align this calss if some types are aligned eigen classes...?

namespace root_solver{

  /// interface to implement operator<< behaviour
  class ostream_insertable{
  public:
    virtual void insertMeInto(ostream &out) const = 0;		 
  };

  ostream& operator<< (ostream &out, const ostream_insertable& p){
    p.insertMeInto(out);
    return out;
  }


  /**
   * This is the general interface for the set of equations \f$ f \mathbb{C}^d \to \mathbb{C}^d$ to be solved.
   * Single and multi root solver specify it further.
   *
   */
  template <typename valueT, typename derivativeT>
  class functions{
  public:
    static_assert (std::is_same<typename Eigen::MatrixBase<valueT>::Scalar, typename Eigen::MatrixBase<derivativeT>::Scalar>::value, "Value and derivate type have different underlying scalar types.");
    static_assert (std::is_same<typename Eigen::MatrixBase<valueT>::RealScalar, typename Eigen::MatrixBase<derivativeT>::RealScalar>::value, "Value and derivate type have different underlying real types.");
    typedef typename Eigen::MatrixBase<valueT>::RealScalar real_type;
    typedef typename Eigen::MatrixBase<valueT>::Scalar scalar_type;

    typedef valueT value_type;
    typedef derivativeT  derivative_type;
    typedef functions<valueT,derivativeT>  functions_type;

    virtual valueT calcF()=0;
    virtual derivativeT calcJ()=0;
    virtual void changePoint(const valueT x)=0;
    virtual valueT guessStartPoint()=0;

  };


  /// States of the root solvers. Basically an extensible enum
  class solver_state
  {
  private:
    int state;
  public:
    solver_state(int s):state(s){}

    solver_state& operator =(const solver_state& rhs){
      state=rhs.state;
      return *this;
    }

    friend bool operator< (const solver_state& lhs, const solver_state& rhs){return lhs.state < rhs.state;}
    friend bool operator==(const solver_state& lhs, const solver_state& rhs){return lhs.state == rhs.state;}

    friend bool operator> (const solver_state& lhs, const solver_state& rhs){return rhs < lhs;}
    friend bool operator<=(const solver_state& lhs, const solver_state& rhs){return !(lhs > rhs);}
    friend bool operator>=(const solver_state& lhs, const solver_state& rhs){return !(lhs < rhs);}
    friend bool operator!=(const solver_state& lhs, const solver_state& rhs){return !(lhs == rhs);}

  friend ostream& operator<< (ostream &out, solver_state &s){
    out << s.state;
    return out;
  }


  };

  static const solver_state ERR(-3); ///< ohoh...
  static const solver_state NOINIT(-2); ///< no starting point given
  static const solver_state STUCK(-1); ///< solver not making progress, aborting this run
  static const solver_state SUCCESS(0); ///< root was found within the given accuracy
  static const solver_state INIT(1); ///< fully initialised, no step made
  static const solver_state CONTINUE(2); ///< Solving in progress, required accuracy not yet reached
  static const solver_state REJECT(3); ///< Solving in progress, required accuracy not yet reached



  /**
   * This interface will be specified in the decendend multiRootSolver and singleRootSolver class.
   * See there for a detailed desribtion of how to use them.
   */
  template <typename valueT, typename derivativeT>
  class rootSolver{
  public:
    //typedefs:
    static_assert (std::is_same<typename Eigen::MatrixBase<valueT>::Scalar, typename Eigen::MatrixBase<derivativeT>::Scalar>::value, "Value and derivate type have different underlying scalar types.");
    static_assert (std::is_same<typename Eigen::MatrixBase<valueT>::RealScalar, typename Eigen::MatrixBase<derivativeT>::RealScalar>::value, "Value and derivate type have different underlying real types.");

    typedef typename Eigen::MatrixBase<valueT>::RealScalar real_type;
    typedef typename Eigen::MatrixBase<valueT>::Scalar scalar_type;
    typedef valueT value_type;
    typedef derivativeT  derivative_type;
    typedef functions<valueT,derivativeT> functions_type;
    typedef rootSolver<valueT,derivativeT> root_solver_type;

    /* implementations shall further provide:
       typedef T numerical_type;
       static const int value_dimension = dim;
    */

  protected:

    //state variables for external use:
    solver_state state;
    value_type z;
    value_type f;
    derivativeT J;
    real_type absF;

    //success/abort criteria
    real_type precisionGoal;
    real_type minStepSize;
    real_type minImprovePerStep;
    real_type minRelativeImprovePerStep;

    //user provides actual function:
    functions<value_type,derivativeT> * calc;

  public:

    rootSolver(functions_type * f_):state(NOINIT),precisionGoal(0.0),minStepSize(0.0),minImprovePerStep(0.0),minRelativeImprovePerStep(0.0),calc(f_){
      if(getPrecisionGoal()==0.0){
        setPrecisionGoal(numeric_limits<real_type>::epsilon());
      }
    }

    ~rootSolver(){}


    value_type getLastPoint(){return z;}
    value_type getLastValue(){return f;}
    solver_state getState(){return state;}
    real_type getAbsF(){return absF;}

    /// Setting a new starting point will reset the solver
    virtual void setStartPoint(value_type z_)=0;

    real_type getPrecisionGoal(){return precisionGoal;}
    real_type getMinStepSize(){return minStepSize;}
    real_type getMinImprovePerStep(){return minImprovePerStep;}
    real_type getMinRelativeImprovePerStep(){return minRelativeImprovePerStep;}

    /// \f$ |f| < precisionGoal \f$ should be the SUCCESS criterion
    virtual void setPrecisionGoal(real_type goal);

    /// if the absolute value of a step is smaller than this value,
    /// the solver is considered as STUCK
    virtual void setMinStepSize(real_type size){minStepSize=size;}

    /// if the improve of the absolute value \f$ |f| \f$ in a step is
    /// smaller than this, the solver is considered as STUCK
    virtual void setMinImprovePerStep(real_type impro){minImprovePerStep=impro;}

    /// if \f$ \frac{|f_{old}|-|f_{new}|}{|f_{old}|} \f$ in a step is
    /// smaller than this, the solver is considered as STUCK
    virtual void setMinRelativeImprovePerStep(real_type relImpro){minRelativeImprovePerStep=relImpro;}

    /// This is the main part of the solver, actually computing a step towards the solution.
    /// Should be run while (step() >= CONTINUE)
    virtual solver_state step()=0;


    /// step may optionally reset the precision goal ...
    virtual solver_state step(real_type epsilonF){
      setPrecisionGoal(epsilonF);
      return step();
    }

    /// ... and stepSize
    virtual solver_state step(real_type epsilonF, real_type epsilonZ){
      setMinStepSize(epsilonZ);
      return step(epsilonF);
    }


    /// Return SUCCESS, CONTINUE or STUCK depending on the
    /// recent progress made and the corresponding thresh holds.
    ///
    /// The default is to consider the solver as STUCK if any of the
    /// stuck criteria is fulfilled.
    virtual solver_state checkStatus(real_type newAbsF, real_type oldAbsF, real_type absDeltaZ);

  };



  //------------------------------------------------------------------------------------------------



  template <typename valueT, typename derivativeT>
  solver_state rootSolver<valueT, derivativeT>::checkStatus(real_type newAbsF, real_type oldAbsF, real_type absDeltaZ){
    if(newAbsF<=getPrecisionGoal()){
#if DEBUG>=SPAM
      cout << __FILE__ << " : SUCCESS"<<endl;
#endif
      return SUCCESS;
    }
    if(absDeltaZ < minStepSize){
#if DEBUG>=SPAM
      cout << __FILE__ << " : STUCK due to stepSize " << absDeltaZ << " < " << minStepSize << endl;
#endif
      return STUCK;
    }
    if(newAbsF > oldAbsF){
#if DEBUG>=SPAM
      cout << __FILE__ << " : REJECT"<<endl;
#endif
      return REJECT;
    }
    real_type deltaAbsF = oldAbsF - newAbsF;
    if(deltaAbsF < minImprovePerStep || deltaAbsF / oldAbsF < minRelativeImprovePerStep){
#if DEBUG>=SPAM
      cout << __FILE__ << " : STUCK due to |f| improve " << deltaAbsF << " < " << minImprovePerStep << " or relative improve " << deltaAbsF / oldAbsF  << " < " << minRelativeImprovePerStep << endl;
#endif
      return STUCK;
    }

#if DEBUG>=SPAM
    cout << __FILE__ << " : CONTINUE"<<endl;
#endif
    return CONTINUE;
  }


  //------------------------------------------------------------------------------------------------


  template <typename valueT, typename derivativeT>
  void rootSolver<valueT, derivativeT>::setPrecisionGoal(real_type goal){
    precisionGoal=goal;
    if(getMinStepSize()==0 && getMinImprovePerStep()==0&& getMinRelativeImprovePerStep()==0){//no stuck criterion set
      setMinRelativeImprovePerStep(1e-2);//default relative improve criterion is 1%
      setMinImprovePerStep(numeric_limits<real_type>::epsilon()*1e5);//much larger than current machine precision
    }
  }




}//end namespace root_solver

#endif
