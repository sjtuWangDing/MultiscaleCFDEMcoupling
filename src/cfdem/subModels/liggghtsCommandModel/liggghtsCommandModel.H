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
  liggghtsCommandModel
\*---------------------------------------------------------------------------*/

#ifndef __LIGGGHTS_COMMAND_MODEL_H__
#define __LIGGGHTS_COMMAND_MODEL_H__

#include "cloud/cfdemCloud.H"
#include "tools/runTimeSelectionTables.H"

namespace Foam {

class liggghtsCommandModel {

public:

  //! @brief Runtime type information
  cfdemTypeName("liggghtsCommandModel");

  //! @brief Declare runtime constructor selection
  cfdemDeclareRunTimeSelection(std::unique_ptr, liggghtsCommandModel, (cfdemCloud& cloud), (cloud));

  //! @brief Selector
  static std::unique_ptr<liggghtsCommandModel> New(cfdemCloud& cloud,
                                                   const dictionary& dict,
                                                   const std::string& modelName);

  //! @brief Constructor
  liggghtsCommandModel(cfdemCloud& cloud);

  //! @brief Destructor
  virtual ~liggghtsCommandModel();

public:

  virtual std::string getCommand(int index) = 0;

public:

  inline int commandLines() const { return commandLines_; }

protected:

  cfdemCloud& cloud_;

  int commandLines_;
};

} // namespace Foam

#endif // __LIGGGHTS_COMMAND_MODEL_H__
