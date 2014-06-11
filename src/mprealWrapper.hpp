/**
 * @file mprealWrapper.hpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0.2014-06-11
 *
 * @section DESCRIPTION
 *
 * Specialisations for mpfr::mpreal in the root_solver environment.
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



#ifndef MPREALWRAPPER_HPP
#define MPREALWRAPPER_HPP

#include "util.hpp"
#include <mpreal.h>
#include <eigen3/Eigen/Dense>

#ifndef MPREALPRECISION
#define MPREALPRECISION 256
#endif

namespace root_solver{

  /// needed to set mpreal precision in each thread at run time (addressing https://github.com/Echsecutor/rootSolver/issues/9 )
  template<>
  class initScalarType<mpfr::mpreal>{
  public:
    static void ini(){
      mpfr::mpreal::set_default_prec(MPREALPRECISION);
    }
  };

}


namespace Eigen{
  namespace internal{

/****************************************************************************
* Implementation of cast to int                                             *
****************************************************************************/
    template<>
    struct cast_impl<mpfr::mpreal,int>
    {
      static inline int run(const mpfr::mpreal& x)
      {
	return static_cast<int>(x.toDouble());
      }
    };

  }
}

#endif
