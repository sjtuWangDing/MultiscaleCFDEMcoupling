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
  forceSubModel
\*---------------------------------------------------------------------------*/

#include "forceSubModel.H"
#include "forceModel.H"

namespace Foam {

//! \brief Constructor
forceSubModel::forceSubModel(cfdemCloud& cloud, forceModel& forceModel):
  cloud_(cloud),
  forceModel_(forceModel),
  switches_(),
  nu_(
    IOobject(
      "scalarViscosity",
      cloud.mesh().time().timeName(),
      cloud.mesh(),
      IOobject::NO_READ,
      IOobject::NO_WRITE
    ),
    cloud.mesh(),
    dimensionedScalar("nu0", dimensionSet(0, 2, -1, 0, 0), 1.)
  ) {}

//! \brief Destructor
forceSubModel::~forceSubModel() {}

/*!
 * \param index                  <[in] 颗粒索引
 * \param dragTot                <[in] 索引为 index 的颗粒受到的总阻力
 * \param dragEx                 <[in] 索引为 index 的颗粒受到的显式阻力
 * \param Ufluid = vector::zero  <[in] 索引为 index 的颗粒中心处流体速度(可以指定是否使用插值模型计算)
 * \param scalar Cd = 0          <[in] 颗粒阻力系数
 */
void forceSubModel::partToArray(const label& index,
                                const Foam::vector& dragTot,
                                const Foam::vector& dragEx,
                                const Foam::vector& Ufluid,
                                scalar Cd) const {
  if (false == switches_.isTrue(treatForceDEM)) {
    // CFD 与 DEM 求解器都考虑耦合力
    if (switches_.isTrue(treatForceExplicit)) {
      // 耦合力视为显式力
      #pragma unroll
      for (int j = 0; j < 3; ++j) {
        // 将耦合力累加到 expFoces_
        forceModel_.expForces()[index][j] += dragTot[j];
      }
    } else {
      // 耦合力视为隐式力
      #pragma unroll
      for (int j = 0; j < 3; ++j) {
        // 将耦合力累加到 impForces_ 和 expFoces_
        forceModel_.impForces()[index][j] += dragTot[j] - dragEx[j];  // 隐式力 = dragTot[j] - dragEx[j]
        forceModel_.expForces()[index][j] += dragEx[j];
      }
    }
  }
  // implForceDEM - true:
  // 颗粒中心处的流体速度和阻力系数都被传递到 DEM 中，从而在每个 DEM 时间步中，使用阻力系数和流体速度，与当前颗粒速度一起计算颗粒受到的阻力
  if (switches_.isTrue(implForceDEM)) {
    #pragma unroll
    for (int j = 0; j < 3; j++) {
      forceModel_.fluidVel()[index][j] = Ufluid[j];
    }
    forceModel_.cds()[index] = Cd;
  } else {
    // 直接将总阻力传递给 DEMForces_[index]
    #pragma unroll
    for (int j = 0; j < 3; j++) {
      forceModel_.DEMForces()[index][j] = dragTot[j];
    }
  }
}

void forceSubModel::readSwitches() {
  for (int i = 0; i < Switches::Num; ++i) {
    Info << "Foam::forceSubModel::readSwitches(): looking for " << Switches::NameList[i] << "...";
    
  }
}

} // namespace Foam
