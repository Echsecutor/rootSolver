/**
 * @file rootSolver.hpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0-2014-03-31
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

#include <iostream>
#include <string>
using namespace std;


//todo: do we need to align this calss if some types are aligned eigen classes...?

/**
 * This is the general interface for the set of equations \f$ f \mathbb{C}^d \to \mathbb{C}^d$ to be solved.
 * Single and multi root solver specify it further.
 *
 */
template <typename valueT, typename derivativeT>
class functions{
public:

  typedef valueT value_type;
  typedef derivativeT  derivative_type;

  virtual valueT calcF()=0;
  virtual derivativeT calcJ()=0;
  virtual void changePoint(const valueT x)=0;
  virtual valueT guessStartPoint()=0;
};


/// States of the root solvers
enum solver_state
  {
    ERR=-3, ///< ohoh...
    NOINIT, ///< no starting point given
    STUCK, ///< solver not making progress, aborting this run
    SUCCESS, ///< root was found within the given accuracy
    INIT, ///< fully initialised, no step made
    CONTINUE ///< Solving in progress, required accuracy not yet reached
  };


/**
 * This interface will be specified in the decendend multiRootSolver and singleRootSolver class.
 * See there for a detailed desribtion of how to use them.
 */
template <typename valueT, typename derivativeT>
class rootSolver{
protected:

  //state variables for external use:
  solver_state state;
  valueT z;
  valueT f;
  derivativeT J;
  double absF;

  //user provides actual function:
  functions<valueT,derivativeT> * calc;

  /// for implementing the "virtual friend <<"
  virtual void writeDatToStream(ostream &out)=0;

public:
  typedef valueT value_type;
  typedef derivativeT  derivative_type;
  typedef functions<valueT,derivativeT> functions_type;
  typedef rootSolver<valueT,derivativeT> root_solver_type;

  /* implementations shall further provide:
  typedef T numerical_type;
  static const int value_dimension = dim;
  */

  rootSolver(functions_type * f_):state(NOINIT),calc(f_){}
  ~rootSolver(){}


  value_type getLastPoint(){return z;}
  value_type getLastValue(){return f;}
  solver_state getState(){return state;}
  double getAbsF(){return absF;}

  virtual string getLastPointInDatFormat()=0;///< shall return "real(z[0])\timag(z[0])\treal(z[1])\t...\timag(z[n])"

  /// same as getLastPointInDatFormat but respects stream formatting etc.
  friend ostream& operator<< (ostream &out, rootSolver &rs){
    rs.writeDatToStream(out);
    return out;
  }

  /// Setting a new starting point will reset the solver
  virtual void setStartPoint(valueT z_)=0;

  /// This is the main part of the solver, actually computing a step towards the solution.
  /// Should be run while (step(epsilonF,epsilonZ) == MRS_CONTINUE)
  /// The first parameter sets a scale for f, defining the minimal improvement per step and maximal allowed final value for \f$ ||f||^2\f$.
  /// The first parameter similarly sets a scale for z, steps smaller than epsilonZ will be counted as stuck.
  virtual solver_state step(double epsilonF,double epsilonZ)=0;
  virtual solver_state step(double epsilonF)=0;

};



#endif
