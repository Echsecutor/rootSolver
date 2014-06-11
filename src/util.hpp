/**
 * @file util.hpp
 * @author Sebastian Schmittner <sebastian@schmittner.pw>
 * @version 1.0.2014-06-04
 *
 *
 * @section DESCRIPTION
 *
 * Utilities, i.e. templates to be specialised appropriately.
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


#ifndef UTIL_HPP
#define UTIL_HPP

namespace root_solver{

  /// if the scalar type needs some run time initialisation in can be wrapped in a specialisation
  template<typename scalarT>
  class initScalarType{
  public:
    static void ini(){}
  };

}


#endif
