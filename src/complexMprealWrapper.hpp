/**
 * @file complexMprealWrapper.hpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0.2014-06-11
 *
 *
 * @section DESCRIPTION
 *
 * Specialisations for root_solver::complex<mpfr::mpreal> in the root_solver environment.
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


#ifndef COMPLEXMPREALWRAPPER_HPP
#define COMPLEXMPREALWRAPPER_HPP

#include "util.hpp"
#include "mprealWrapper.hpp"
#include "complex.hpp"

namespace root_solver{

  template<>
  class initScalarType<typename root_solver::complex<mpfr::mpreal> >{
  public:
    static void ini(){
      root_solver::initScalarType<mpfr::mpreal>::ini();
      root_solver::complex<mpfr::mpreal>::PI = mpfr::const_pi(MPREALPRECISION, mpfr::mpreal::get_default_rnd());
    }
  };

}

#endif
