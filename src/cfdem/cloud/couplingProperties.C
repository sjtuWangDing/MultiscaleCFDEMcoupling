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
  Foam::CouplingProperties
\*---------------------------------------------------------------------------*/

#include "couplingProperties.H"

const std::string CFDEM_VERSION = "cfdem-3.8.1";
const std::string LIGGGHTS_VERSION = "3.8.0";

namespace Foam {

CouplingProperties::CouplingProperties(const fvMesh& mesh,
                                       const IOdictionary& couplingPropertiesDict,
                                       const IOdictionary& liggghtsCommandsDict):
  mesh_(mesh),
  couplingPropertiesDict_(couplingPropertiesDict),
  liggghtsCommandsDict_(liggghtsCommandsDict),
  modelType_(couplingPropertiesDict.lookup("modelType")),
  forceModels_(couplingPropertiesDict.lookup("forceModels")),
  momCoupleModels_(couplingPropertiesDict.lookup("momCoupleModels")),
  liggghtsCommandModels_(couplingPropertiesDict.lookup("liggghtsCommandModels")),
  turbulenceModelType_(couplingPropertiesDict.lookup("turbulenceModelType")),
  debug_(couplingPropertiesDict.lookupOrDefault<Switch>("debug", false)),
  ignore_(couplingPropertiesDict.lookupOrDefault<Switch>("ignore", false)),
  solveFlow_(couplingPropertiesDict.lookupOrDefault<Switch>("solveFlow", true)),
  allowAdjustTimeStep_(couplingPropertiesDict.lookupOrDefault<Switch>("allowAdjustTimeStep", false)),
  verbose_(couplingPropertiesDict.lookupOrDefault<Switch>("verbose", false)) {

  Info << "CFDEM coupling version: " << CFDEM_VERSION <<  endl;
  Info << "LIGGGHTS version: " << LIGGGHTS_VERSION <<  endl;

#if __MIXCLOUD__

  fineParticleRatio_ = couplingPropertiesDict_.lookupOrDefault<double>("fineParticleRatio", 3.0);
  coarseParticleRatio_ = couplingPropertiesDict_.lookupOrDefault<double>("coarseParticleRatio", 0.33);
  usedForSolverIB_ = couplingPropertiesDict_.lookupOrDefault<Switch>("usedForSolverIB", false);
  usedForSolverPiso_ = couplingPropertiesDict_.lookupOrDefault<Switch>("usedForSolverPiso", false);
  useDynamicRefineMesh_ = couplingPropertiesDict_.lookupOrDefault<Switch>("useDynamicRefineMesh", false);

  // 读取 fixed_particle_ 与来流速度
  fixedParticle_ = couplingPropertiesDict_.lookupOrDefault<bool>("fixedParticle", false);
  flowVelocity_ = vector(couplingPropertiesDict_.lookupOrDefault<double>("U0x", 0.0),
                         couplingPropertiesDict_.lookupOrDefault<double>("U0y", 0.0),
                         couplingPropertiesDict_.lookupOrDefault<double>("U0z", 0.0));

  if (fixedParticle_) {
    Info << "Using fixed particle, reading flow_velocity_: " << flowVelocity_ << endl;
  }

  if (mag(flowVelocity_) < SMALL) {
    FatalError << "mag(flowVelocity_) is zero" << abort(FatalError);
  }

#endif // __MIXCLOUD__
}

} // namespace Foam
