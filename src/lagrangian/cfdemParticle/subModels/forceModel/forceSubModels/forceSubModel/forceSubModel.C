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
#include "forceSubModel.H"
#include "forceModel.H"
#include "mathExtra.H"

namespace Foam {

defineTypeNameAndDebug(forceSubModel, 0);

defineRunTimeSelectionTable(forceSubModel, dictionary);

// @brief Constructors
forceSubModel::forceSubModel(const dictionary& dict, cfdemCloud& sm, forceModel& fm):
  dict_(dict),
  particleCloud_(sm),
  forceModel_(fm),
  // Note: 在 forceSubModel 中默认的 switch 共 11 个
  nrDefaultSwitches_(11),
  switchesNameList_(wordList(nrDefaultSwitches_)),
  switchesList_(List<Switch>(nrDefaultSwitches_)),
  switches_(List<Switch>(nrDefaultSwitches_)),
  nu_
  (
    IOobject
    (
      "scalarViscosity",
      sm.mesh().time().timeName(),
      sm.mesh(),
      IOobject::NO_READ,
      IOobject::NO_WRITE
    ),
    sm.mesh(),
    dimensionedScalar("nu0", dimensionSet(0, 2, -1, 0, 0), 1.)
  ),
  divTau_
  (
    IOobject
    (
      "divTau",
      sm.mesh().time().timeName(),
      sm.mesh(),
      IOobject::NO_READ,
      IOobject::NO_WRITE
    ),
    sm.mesh(),
    dimensionedVector("divTau", dimensionSet(1, -2, -2, 0, 0), vector::zero)
  ),
  IBDragPerV_
  (
    IOobject
    (
      "IBDragPerV",
      sm.mesh().time().timeName(),
      sm.mesh(),
      IOobject::NO_READ,
      IOobject::NO_WRITE
    ),
    sm.mesh(),
    dimensionedVector("IBDragPerV", dimensionSet(1, -2, -2, 0, 0), vector::zero)
  ),
  densityFieldName_(dict_.lookupOrDefault<word>("densityFieldName", "rho")),
  rho_(sm.mesh().lookupObject<volScalarField>(densityFieldName_)),
  verboseDiskIntervall_(1),
  verboseDiskCounter_(0),
  scaleDia_(dict_.lookupOrDefault<scalar>("scale", 1.)),
  scaleDrag_(dict_.lookupOrDefault<scalar>("scaleDrag", 1.)),
  scaleDH_(dict_.lookupOrDefault<scalar>("scaleDH", 1.)) {

  // init standard switch list
  switchesNameList_[0] = "treatForceExplicit";
  switchesNameList_[1] = "treatForceDEM";
  switchesNameList_[2] = "implForceDEM";
  switchesNameList_[3] = "verbose";
  switchesNameList_[4] = "interpolation";
  switchesNameList_[5] = "useFilteredDragModel";
  switchesNameList_[6] = "useParcelSizeDependentFilteredDrag";
  switchesNameList_[7] = "implForceDEMaccumulated";
  switchesNameList_[8] = "scalarViscosity";
  switchesNameList_[9] = "verboseToDisk";
  switchesNameList_[10] = "useCorrectedVoidage";

  // info about scaleDia being used
  if (scaleDia_ != 1) {
    Info << "forceSubModel::forceSubModel(): using scale = " << scaleDia_ << endl;
  } else if (particleCloud_.cg() != 1) {
    scaleDia_ = particleCloud_.cg();
    Info << "forceSubModel::forceSubModel(): using scale from liggghts cg = " << scaleDia_ << endl;
  }

  // info about scaleDrag being used
  if (scaleDrag_ != 1.) {
    Info << "forceSubModel::forceSubModel(): using scaleDrag = " << scaleDrag_ << endl;
  }
  if (scaleDrag_ < SMALL) {
    FatalError<< "forceSubModel::forceSubModel(): scaleDrag > 0 required" << abort(FatalError);
  }
  particleCloud_.registryM().addProperty("scaleDrag", scaleDrag_);

  // info about scaleDH being used
  if (scaleDH_ != 1) {
    Info << "forceSubModel::forceSubModel(): using scaleDH = " << scaleDH_ << endl;
  }
}

// @brief Destructor
forceSubModel::~forceSubModel() {}

// @param index                  <[in] 颗粒索引
// @param dragTot                <[in] 索引为 index 的颗粒受到的总阻力
// @param dragEx                 <[in] 索引为 index 的颗粒受到的显式阻力
// @param Ufluid = vector::zero  <[in] 索引为 index 的颗粒中心处流体速度(可以指定是否使用插值模型计算)
// @param scalar Cd = 0          <[in] 颗粒阻力系数
void forceSubModel::partToArray(const label& index,
                                const vector& dragTot,
                                const vector& dragEx,
                                const vector& Ufluid,
                                scalar Cd) const {

  if (switches_[1] == false) {
    // CFD 与 DEM 求解器都考虑耦合力
    if (switches_[0] == true) {
      // 耦合力视为显式力
      for (int j = 0; j < 3; ++j) {
        // 将耦合力累加到 particleCloud_.expFoces_
        myForceM().expForces()[index][j] += dragTot[j];
      }
    } else {
      // 耦合力视为隐式力
      for (int j = 0; j < 3; ++j) {
        // 将耦合力累加到 particleCloud_.impForces_ 和 particleCloud_.expFoces_
        myForceM().impForces()[index][j] += dragTot[j] - dragEx[j];  // 隐式力 = dragTot[j] - dragEx[j]
        myForceM().expForces()[index][j] += dragEx[j];
      }
    }
  }

  if (switches_[2] == true) {
    // implForceDEM - true: 颗粒中心处的流体速度和阻力系数都被传递到 DEM 中，从而在每个 DEM 时间步中，使用阻力系数和流体速度，与当前颗粒速度一起计算颗粒受到的阻力
    for (int j = 0; j < 3; j++) {
      myForceM().fluidVel()[index][j] = Ufluid[j];
    }
    myForceM().Cds()[index][0] = Cd;
  } else {
    // 直接将总阻力传递给 particleCloud_.DEMForces_[index]
    for (int j = 0; j < 3; j++) {
      myForceM().DEMForces()[index][j] += dragTot[j];
    }
  }
}

// @brief 计算索引为 index 的颗粒的 scale 直径
// @param d     <[in, out] 颗粒直径
// @param index <[in] 颗粒索引
void forceSubModel::scaleDia(scalar& d, int index) const {
  if (particleCloud_.cgTypeSpecificDifferent) {
    d /= particleCloud_.cg(particleCloud_.particleType(index));
  } else {
    d /= scaleDia_ / scaleDH_;
  }
}

void forceSubModel::scaleForce(vector& force, scalar& d, int index) const {
  if (particleCloud_.cgTypeSpecificDifferent) {
    double cgCurr = particleCloud_.cg(particleCloud_.particleType(index));
    force *= cgCurr * cgCurr * cgCurr;
  } else {
    force *= scaleDia_ * scaleDia_ * scaleDia_;
  }
  force *= scaleDrag_;
}

void forceSubModel::scaleCoeff(scalar& coeff, scalar& d, int index) const {
  if (particleCloud_.cgTypeSpecificDifferent) {
    double cgCurr = particleCloud_.cg(particleCloud_.particleType(index));
    coeff *= cgCurr * cgCurr * cgCurr;
  } else {
    coeff *= scaleDia_ * scaleDia_ * scaleDia_;
  }
  coeff *= scaleDrag_;
}

void forceSubModel::readSwitches() const {

  // 遍历每一个 switch
  for (int i = 0; i < switchesNameList_.size(); ++i) {
    if (switchesList_[i]) {  // 该力模型中指定了第 i 个 switch

      Info << "Foam::forceSubModel::readSwitches(): looking for " << switchesNameList_[i] << "...";

      if (dict_.found(switchesNameList_[i])) {
        // 搜索到了 switchesNameList_[i]
        Info << " found in dict.\n" << endl;
        // 将 switches_[i] 设置为 dict 所指定的
        switches_[i] = Switch(dict_.lookup(switchesNameList_[i]));
      } else {
        Info << " not found in dict, using default.\n" << endl;
      }

      if (i == 0 && treatForceExplicit() && particleCloud_.registryM().getProperty("explicitCouple_index") < 0 &&
          !treatForceDEM() && particleCloud_.modelType() != "none") {

        FatalError << "modelType: " << particleCloud_.modelType() << endl
          << "getProperty(explicitCouple_index): " 
          << particleCloud_.registryM().getProperty("explicitCouple_index") << endl
          << "switch treatForceExplicit: " << switches_[0] << endl
          << "switch treatForceDEM: " << switches_[1] << endl
          << "You are using treatForceExplicit = true here, this requres having an explicit momentum couple model!" 
          << abort(FatalError);
      } else {
        Info << switchesNameList_[i] << " = " << switches_[i] << endl;
      }
    }  // End of switchesList_[i]
  }  // End of forAll(switchesNameList_)

  if (implForceDEM()) {
    // 将 particleCloud 中的 impDEMdrag_ 设置为 true
    particleCloud_.impDEMdrag_ = true;
  }

  if (implForceDEMaccumulated()) {
    if (!implForceDEM()) {
      Warning << "please check your settings, implForceDEMaccumulated without implForceDEM does not work! (using implForceDEMaccumulated = false)" << endl;
      // 将 verbose_ 设置为 true, 强制打印调试信息
      switches_[3] = false;
    } else {
      particleCloud_.impDEMdragAcc_=true;
    }
  }

  if (scalarViscosity()) {
    Info << "Using a constant viscosity for this force model." << endl;

    dimensionedScalar nu0_("nu", dimensionSet(0, 2, -1, 0, 0), dict_.lookup("nu"));
    nu_ = volScalarField(IOobject("scalarViscosity",
                                  particleCloud_.mesh().time().timeName(),
                                  particleCloud_.mesh(),
                                  IOobject::NO_READ,
                                  IOobject::NO_WRITE),
                         particleCloud_.mesh(),
                         nu0_);
  }

  dict_.readIfPresent("verboseDiskIntervall", verboseDiskIntervall_);

  // 如果使用旧式语法，则提示查看 forceSubModel doc.
  if (dict_.found("treatExplicit") || dict_.found("treatDEM") || dict_.found("implDEM")) {
    FatalError << "You are using an old nomenclature for force model settings, please have a look at the forceSubModel doc." << abort(FatalError);
  }

  if (dict_.found("verbose")) {
    Warning << "Please make sure you use the new nomenclature for verbose force model settings, please have a look at the forceSubModel doc." << endl;
  }
}

const volScalarField& forceSubModel::nuField() const {
#ifdef compre
  nu_ = particleCloud_.turbulence().mu() / rho_;
  return nu_;
#else
  if (scalarViscosity()) {
    // 使用用户自定义动力粘度
    return nu_;
  } else {
    return particleCloud_.turbulence().nu();
  }
#endif
}

const volScalarField& forceSubModel::muField() const {
#ifdef compre
  return particleCloud_.turbulence().mu();
#else
  if (scalarViscosity()) {
    // 使用用户自定义动力粘度
    // usage of constant mu_ is still commented,
    // as not tested particleCloud_.turbulence().nu() * rho_ does not work properly
    FatalError << "Foam::forceSubModel::muField(): implementation not complete!" << abort(FatalError);
  }
  return particleCloud_.turbulence().nu() * rho_;
#endif
}

const volScalarField& forceSubModel::rhoField() const {
  return rho_;
}

// @brief 计算粘性力
// @note 用于计算 viscForceModel
const volVectorField& forceSubModel::divTauField(const volVectorField& U) const {
#ifdef compre
  const volScalarField& mu_ = muField();
  divTau_ = -fvc::laplacian(mu_, U) - fvc::div(mu_ * dev(fvc::grad(U)().T()));
#else
  const volScalarField& nu_ = nuField();
  const volScalarField& rho_ = rhoField();
  divTau_ = -fvc::laplacian(nu_ * rho_, U) - fvc::div(nu_ * rho_ * dev(fvc::grad(U)().T()));
#endif
  return divTau_;
}

// @brief 计算 IB drag
// @note 用于计算 ShirgaonkarIBModel
const volVectorField& forceSubModel::IBDragPerV(const volVectorField& U, const volScalarField& p) const {
#ifdef compre
  IBDragPerV_ = muField() * fvc::laplacian(U) - fvc::grad(p);
#else
  IBDragPerV_ = rhoField() * (nuField() * fvc::laplacian(U) - fvc::grad(p));
#endif
  return IBDragPerV_;
}

}  // End of namespace Foam
