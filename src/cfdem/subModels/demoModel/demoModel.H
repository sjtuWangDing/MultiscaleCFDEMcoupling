/*---------------------------------------------------------------------------*\
  CFDEMcoupling - Open Source CFD-DEM coupling

  CFDEMcoupling is part of the CFDEMproject
  www.cfdem.com
                              Christoph Goniva, christoph.goniva@cfdem.com
                              Copyright 2009-2012 JKU Linz
                              Copyright 2012-     DCS Computing GmbH, Linz
------------------------------------------------------------------------------
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

Description
  This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
  and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).

Class
  demoModel
\*---------------------------------------------------------------------------*/

#ifndef __DEMO_MODEL_H__
#define __DEMO_MODEL_H__

#include "cloud/cfdemCloud.H"
#include "base/runTimeSelectionTables.H"

namespace Foam {

//! @brief demoModel
class demoModel {

public:
  /* --------------------------------- base class and functions -------------------------- */

  //! @brief Runtime type information
  cfdemTypeName("demoModel")

  //! @brief Declare runtime constructor selection
  cfdemDeclareRunTimeSelection(autoPtr, demoModel, (cfdemCloud& cloud), (cloud))

  //! @brief Selector
  static autoPtr<demoModel> New(cfdemCloud& cloud, const dictionary& dict);

  //! @brief Constructor
  demoModel(cfdemCloud& cloud);

  //! @brief Destructor
  virtual ~demoModel();

  /* --------------------------------- virtual functions --------------------------------- */

  /* --------------------------------- member functions ---------------------------------- */

  /* --------------------------------- public interfaces --------------------------------- */

protected:

  cfdemCloud& cloud_;
};

} // namespace Foam

#endif // __DEMO_MODEL_H__
