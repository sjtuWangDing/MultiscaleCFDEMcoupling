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
  cfdemCloud class managing DEM data for CFD-DEM coupling

Class
  Foam::cfdemCloud
\*---------------------------------------------------------------------------*/

#include "cloud/cfdemCloud.H"
#include "cloud/couplingProperties.H"
#include "cloud/particleCloud.H"
#include "subModels/liggghtsCommandModel/liggghtsCommandModel.H"

namespace Foam {

cfdemCloud::cfdemCloud(const fvMesh& mesh):
  mesh_(mesh),
  couplingPropertiesDict_(
    IOobject(
      "couplingProperties", // coupling properties file name
      mesh.time().constant(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  ),
  liggghtsCommandsDict_(
    IOobject(
      "liggghtsCommands", // liggghts commands file name
      mesh.time().constant(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
    )
  ),
  cProps_(mesh, couplingPropertiesDict_, liggghtsCommandsDict_),
  pCloud_(0) {

  Info << "\nEnding of Constructing cfdemCloud Base Class Object......\n" << endl;
  Info << "\nEntry of cfdemCloud::cfdemCloud(const fvMesh&)......\n" << endl;

  for (int i = 0; i < liggghtsCommandModelList().size(); ++i) {
    liggghtsCommandModel_.emplace_back(
      liggghtsCommandModel::New(*this, liggghtsCommandsDict_, liggghtsCommandModelList()[i]));
  }



  // liggghtsCommandModel_ = new autoPtr<liggghtsCommandModel>[liggghtsCommandModelList().size()];
  // for (int i = 0; i < liggghtsCommandModelList().size(); ++i) {
  //   liggghtsCommandModel_[i] =
  //     liggghtsCommandModel::New(*this, liggghtsCommandsDict_, liggghtsCommandModelList()[i]);
  // }
}

} // namespace Foam
