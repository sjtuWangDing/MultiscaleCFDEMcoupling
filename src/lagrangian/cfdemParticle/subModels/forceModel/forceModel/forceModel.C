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

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"
#include "forceModel.H"
#include "mathExtra.H"

namespace Foam {

defineTypeNameAndDebug(forceModel, 0);

defineRunTimeSelectionTable(forceModel, dictionary);

//! @brief Constructors
forceModel::forceModel(const dictionary& dict, cfdemCloud& sm):
  dict_(dict),
  particleCloud_(sm),
  impParticleForces_
  (
    IOobject
    (
      "impParticleForces",
      sm.mesh().time().timeName(),
      sm.mesh(),
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
    ),
    sm.mesh(),
    dimensionedVector("zero", dimensionSet(1, 1, -2, 0, 0), vector(0, 0, 0))  // [N] == [kg * m / s^2]
  ),
  expParticleForces_
  (
    IOobject
    (
      "expParticleForces",
      sm.mesh().time().timeName(),
      sm.mesh(),
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
    ),
    sm.mesh(),
    dimensionedVector("zero", dimensionSet(1, 1, -2, 0, 0), vector(0, 0, 0))  // [N] == [kg * m / s^2]
  ),
  modelType_(sm.modelType()),
  probeIt_(sm.probeM().active()),
  requiresEx_(false),
  requiresShape_(false),
  requiresQuaternion_(false),
  requiresSuperquadric_(false),
  pullPushRotation_(false),
  implicitDrag_(false),
  implicitAnisotropicDrag_(false),
  implicitRotation_(false),
  forceSubModels_(0),
  forceSubModel_(new autoPtr<forceSubModel>[nrForceSubModels()]),
  voidfractionInterpolator_(NULL),
  UInterpolator_(NULL),
  vorticityInterpolator_(NULL),
  gradPInterpolator_(NULL),
  gradUInterpolator_(NULL),
  gradVoidfractionInterpolator_(NULL),
  Up1Interpolator_(NULL),
  Up2Interpolator_(NULL),
  dSauterInterpolator_(NULL),
  phiP1Interpolator_(NULL),
  phiP2Interpolator_(NULL),
  alphaInterpolator_(NULL),
  gradAlphaInterpolator_(NULL),
  TInterpolator_(NULL),
  UsInterpolator_(NULL),
  fluidScalarFieldInterpolator_(NULL),
  gradPsolidInterpolator_(NULL),
  shearRateInterpolator_(NULL),
  DDtUInterpolator_(NULL),
  divTauInterpolator_(NULL) {}

//! @brief Destructor
forceModel::~forceModel() {}

void forceModel::applyDebugSettings(bool debug) const {
  if(!debug) {
    impParticleForces_.writeOpt() = IOobject::NO_WRITE;
    expParticleForces_.writeOpt() = IOobject::NO_WRITE;
  }
}

//! @brief 重新划分 Implicit / Explicit force
void forceModel::repartitionImExForces() const {
  if (particleCloud_.imExSplitFactor() < 1.0) {
    Info << "Will re-partition split of implicit and explicit forces: alpha = " << particleCloud_.imExSplitFactor() << endl;
    // Update implicit particle
    expParticleForces_ += (1.0 - particleCloud_.imExSplitFactor()) * impParticleForces_;
    impParticleForces_ *= particleCloud_.imExSplitFactor();
  }
}

//! @brief 对于不包含任何颗粒的网格使用显式力耦合
void forceModel::treatVoidCells() const {
  if (particleCloud_.treatVoidCellsAsExplicitForce()) {
    int counter(0);
    volVectorField& Us = particleCloud_.averagingM().UsNext();
    forAll(Us, cellI) {
      if (mag(Us[cellI]) == 0.0) {  // 如果网格的局部平均速度场为 0.0，则说明网格中没有小颗粒
        expParticleForces_[cellI] += impParticleForces_[cellI];
        impParticleForces_[cellI] *= 0.0;
        counter +=1;
      }
    }
    Info << "Re-partitioned "<< counter << " cells void of particles" << endl;
  }
}

//! @brief 设置 force sub model
//! @brief 比如对于 DiFeliceDrag force model, setForceSubModels 会从 DiFeliceDragProps 中寻找并设置所有的 forceSubModels 字段
//        如果 DiFeliceDragProps 中没有指定 forceSubModels 字段, 则会默认设置为 ImEx
// @note 每一类型的力都可以指定多个种类的 forceSubModel, 但是目前只有一种 forceSubModel 的实现, 即 ImEx
void forceModel::setForceSubModels(dictionary& dict) {
  if (dict.found("forceSubModels")) {
    // 读入字典文件中指定的所有 force Sub Models 的名称
    forceSubModels_ = wordList(dict.lookup("forceSubModels"));
  } else if (dict.found("forceSubModel")) {
    FatalError << "Did you mean the forceSubModels keyword? " << abort(FatalError);
  } else {
    forceSubModels_.setSize(1, "ImEx");
  }

  delete[] forceSubModel_;
  forceSubModel_ = new autoPtr<forceSubModel>[nrForceSubModels()];
  for (int i = 0; i < nrForceSubModels(); i++) {
    forceSubModel_[i] = forceSubModel::New(dict, particleCloud_, *this, forceSubModels_[i]);
  }
}

}  // End of namespace Foam
