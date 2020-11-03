/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
\*---------------------------------------------------------------------------*/

#ifndef __TYPEINFO_H__
#define __TYPEINFO_H__

#include <stdexcept>
#include <typeinfo>

namespace base {

// //! @brief Reference type cast template function
// //         wraps dynamic_cast to handle bad_cast exception and generate a FatalError.
// template<typename ToType, typename FromType>
// inline ToType& dynamicCast(FromType& from) {
//   try {
//     return dynamic_cast<ToType&>(from);
//   } catch(const std::bad_cast& ex) {
//     FatalErrorInFunction << "Attempt to cast type " << typeid(FromType).name()
//       << " to type " << typeid(ToType).name() << abort(FatalError);
//     return dynamic_cast<ToType&>(from);
//   }
// }

} // namespace base

#endif // __TYPEINFO_H__
