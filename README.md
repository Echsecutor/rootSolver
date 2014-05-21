rootSolver
==========

This is a solver (frame work) for solving (complex) one or multi dimensional (non-linear) equations. It is written generically, in particular it can easily be used for complex functions and variables, other than the GSL variant http://www.gnu.org/software/gsl/manual/html_node/Multidimensional-Root_002dFinding.html  from which this project is completely independend.

Internally, [eigen](http://eigen.tuxfamily.org) is used to solve the linear equations.

## Version

1.0.2014-05-15

The Version number is formatted as "M.S.D" where M is the major release branch (backward compatibility to all non-alpha releases of the same branch is guaranteed), S is the state of this release (0 for alpha, 1 for beta, 2 for stable), and D is the date formatted as yyyy-mm-dd.


## How-To

To make use of this library you will need to include the appropriate *Solver.hpp file into your projetc. These files contain a description of how to run the respective solver. Generically you will need to derive the apropriate *functions class (interface) to implement your function, instanciate the solver passing your function, set a starting poin, and then run it like `while(Solver.step(epsilon) == CONTINUE)`.

## Solver traits

file|description
---|---
singleRootSolver.hpp | This is a simple implementation of Newtons method for one complex variable/equation using back traking.
multiRootSolver.hpp | Complex multidimensional root finding via Newtons method.
extrapolationSolver.hpp | Extends one of the elementary solvers to track the root of a function that changes continuously with an external parameter.
batchSolver.hpp | Wrapper to solve a set of equations, i.e. computing implicit functions.

Batch and extrapolation solver can be combined to give the most powerful solver implemented so far.

See the test cases for examples.


## ToDo

Before I will call a release a "beta", I want to do (at least) the following

* Actually test wether the results of test runs make sense. Currently test run through unless there are compile/runtime errors.
* A steepest gradient descent method +dogleg step is partially implemented for the multiRootSolver. These should be completed and tested.
* Whenever you see a shared_ptr I am usually working around [alignment issues with eigen](http://eigen.tuxfamily.org/dox/group__TopicUnalignedArrayAssert.html). In most cases, alignment will not give a large speed advantage. Hence get rid of unneccesary alignments by just usign the dynamic eigen classes, wherever appropriate.
