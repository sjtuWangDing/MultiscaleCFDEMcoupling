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
  forceModel
\*---------------------------------------------------------------------------*/

#ifndef __FORCE_MODEL_H__
#define __FORCE_MODEL_H__

#include "cloud/cfdemCloud.H"
#include "base/runTimeSelectionTables.H"
#include "subModels/forceModel/forceSubModel.H"

namespace Foam {

//! \brief forceModel
class forceModel {

public:
  /* --------------------------------- base class and functions -------------------------- */

  //! \brief Runtime type information
  cfdemTypeName("forceModel")

  //! \brief Declare runtime constructor selection
  cfdemDeclareRunTimeSelection(std::unique_ptr, forceModel, (cfdemCloud& cloud), (cloud))

  //! \brief Selector
  static std::unique_ptr<forceModel> New(cfdemCloud& cloud,
                                         const dictionary& dict,
                                         const std::string& modelName);

  //! \brief Constructor
  forceModel(cfdemCloud& cloud);

  //! \brief Destructor
  virtual ~forceModel();

  /* --------------------------------- virtual functions --------------------------------- */

  virtual void setForce() {
    FatalError << "forceModel:setForce(): using base class function, please use derived class function\n"
      << abort(FatalError);
  }

  // virtual void setMixForce(const std::vector<double>& dimensionRatios) {
  //   FatalError << "forceModel:setMixForce(): using base class function, please use derived class function\n"
  //     << abort(FatalError);
  // }

  /* --------------------------------- member functions ---------------------------------- */

  /* --------------------------------- public interfaces --------------------------------- */

  inline base::CDTensor2& impForces() { return cloud_.pCloud().impForces(); }
  inline base::CDTensor2& expForces() { return cloud_.pCloud().expForces(); }
  inline base::CDTensor2& DEMForces() { return cloud_.pCloud().DEMForces(); }
  inline base::CDTensor2& fluidVel() { return cloud_.pCloud().fluidVel(); }
  inline base::CDTensor1& cds() { return cloud_.pCloud().cds(); }

protected:

  cfdemCloud& cloud_;

  //! \brief 当前 force model 的 forceSubModel
  forceSubModel forceSubModel_;

  //! \brief 是否激活探针
  bool useProbe_;

  //! \brief 颗粒隐式力的总和 [N]
  volVectorField impParticleForces_;

  //! \brief 颗粒显式力的总和 [N]
  volVectorField expParticleForces_;
};

} // namespace Foam

#endif // __FORCE_MODEL_H__
