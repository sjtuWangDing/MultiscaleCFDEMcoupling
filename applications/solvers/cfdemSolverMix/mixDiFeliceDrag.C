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
#include "mixDiFeliceDrag.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam {

defineTypeNameAndDebug(mixDiFeliceDrag, 0);

addToRunTimeSelectionTable(forceModel, mixDiFeliceDrag, dictionary);

// @brief Constructors
mixDiFeliceDrag::mixDiFeliceDrag(const dictionary& dict, cfdemCloud& sm):
  forceModel(dict, sm),
  propsDict_(dict.subDict(typeName + "Props")),
  velFieldName_(propsDict_.lookup("velFieldName")),
  U_(sm.mesh().lookupObject<volVectorField>(velFieldName_)),
  voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
  voidfraction_(sm.mesh().lookupObject<volScalarField>(voidfractionFieldName_)),
  UsFieldName_(propsDict_.lookup("granVelFieldName")),
  UsField_(sm.mesh().lookupObject<volVectorField>(UsFieldName_)) {

  // 设置激活 or 抑制颗粒探针
  if (probeIt_ && propsDict_.found("suppressProbe")) {
    probeIt_ = !Switch(propsDict_.lookup("suppressProbe"));
  }

  if (probeIt_) {
    particleCloud_.probeM().initialize(typeName, typeName + ".logDat");
    particleCloud_.probeM().vectorFields_.append("dragForce");     // first entry must the be the force
    particleCloud_.probeM().vectorFields_.append("Urel");          // other are debug
    particleCloud_.probeM().scalarFields_.append("Rep");           // other are debug
    particleCloud_.probeM().scalarFields_.append("Cd");            // other are debug
    particleCloud_.probeM().scalarFields_.append("voidfraction");  // other are debug
    particleCloud_.probeM().writeHeader();
  }

  particleCloud_.checkCG(true);

  // init force sub model
  setForceSubModels(propsDict_);

  // define switches which can be read from dict
  forceSubM(0).setSwitchesList(0, true);  // activate treatExplicit switch
  forceSubM(0).setSwitchesList(2, true);  // activate implDEM switch
  forceSubM(0).setSwitchesList(3, true);  // activate search for verbose switch
  forceSubM(0).setSwitchesList(4, true);  // activate search for interpolate switch
  forceSubM(0).setSwitchesList(8, true);  // activate scalarViscosity switch

  // read those switches defined above, if provided in dict
  forceSubM(0).readSwitches();
  // for (int iFSub = 0; iFSub < nrForceSubModels(); iFSub++) {
  //   forceSubM(iFSub).readSwitches();
  // }
}

// @brief Destructor
mixDiFeliceDrag::~mixDiFeliceDrag() {}

void mixDiFeliceDrag::calForceKernel(const int& index,
                                     const int& cellI,
                                     const volScalarField& nufField,
                                     const volScalarField& rhoField,
                                     vector& drag, vector& Ufluid, scalar& dragCoefficient,
                                     scalar& voidfraction, vector& Ur,
                                     scalar& Rep, scalar& Cd) const {
  vector position(0, 0, 0);      // 颗粒位置
  vector Us(0, 0, 0);            // 局部平均颗粒速度
  scalar magUr(0);               // 相对速度值
  scalar dReal(0);               // 颗粒的实际直径
  scalar dScale(0);              // 颗粒的 scale 直径
  scalar nuf(0);                 // 流体动力粘度
  scalar rho(0);                 // 流体密度
  scalar Xi(0);                  // 模型阻力系数 Xi = 3.7 - 0.65 * exp(-sqr(1.5 - log10(Rep)) / 2)

  drag = vector::zero;           // 总阻力 = dragCoefficient * Ur
  Ufluid = vector::zero;         // 颗粒中心处流体速度(可以指定是否使用插值模型计算)
  voidfraction = 1;              // 颗粒中心所在网格的空隙率(可以指定是否使用插值模型计算)
  Ur = vector::zero;             // 相对速度
  Rep = 0;                       // 颗粒雷诺数
  Cd = 0;                        // 流体阻力系数 Cd = sqr(0.63 + 4.8 / sqrt(Rep))
  dragCoefficient = 0;           // 颗粒阻力系数

  if (forceSubM(0).interpolation()) {  // 使用插值模型将欧拉场插值到拉格朗日场
    // 获取颗粒中心的坐标, 将颗粒中心所在网格的空隙率和流体速度插值到颗粒中心处
    position = particleCloud_.position(index);
    voidfraction = voidfractionInterpolator_().interpolate(position, cellI);
    Ufluid = UInterpolator_().interpolate(position, cellI);

    // 确保插值后颗粒中心的空隙率有意义
    if (voidfraction > 1.00) { voidfraction = 1.00; }
    if (voidfraction < SMALL) { voidfraction = SMALL; }
  } else {
    // 不使用插值模型，直接设置颗粒中心处的值为颗粒中心所在网格的值
    voidfraction = voidfraction_[cellI];
    Ufluid = U_[cellI];
  }

  Us = particleCloud_.velocity(index);
  dReal = 2 * particleCloud_.radius(index);
  dScale = 2 * particleCloud_.radius(index);
  forceSubM(0).scaleDia(dScale, index);  // 计算颗粒的 scale 直径

  // 初始化
  Ur = Ufluid - Us;
  magUr = mag(Ur);
  nuf = nufField[cellI];
  rho = rhoField[cellI];
  Rep = 0;
  Cd = 0;
  dragCoefficient = 0;

  if (magUr > 0) {
    // 计算颗粒雷诺数
    Rep = dScale * voidfraction * magUr / (nuf + SMALL);
    // 计算流体阻力系数
    Cd = sqr(0.63 + 4.8 / sqrt(Rep));
    // 计算模型阻力系数
    Xi = 3.7 - 0.65 * exp(-sqr(1.5 - log10(Rep)) / 2);
    // 计算颗粒阻力系数
    dragCoefficient = 0.125 * Cd * rho * M_PI * dScale * dScale * pow(voidfraction, (2 - Xi)) * magUr;
    if (modelType_ == "B") {
      dragCoefficient /= voidfraction;
    }
    forceSubM(0).scaleCoeff(dragCoefficient, dReal, index);
    // 计算总阻力
    drag = dragCoefficient * Ur;
  }
}

// Member Functions
void mixDiFeliceDrag::setForce() const {
  Info << "Setting mixDiFeliceDrag force..." << endl;

  const volScalarField& nufField = forceSubM(0).nuField();
  const volScalarField& rhoField = forceSubM(0).rhoField();
  label cellI(0);                // 颗粒中心所在网格的索引
  vector drag(0, 0, 0);          // 总阻力 = dragCoefficient * Ur
  vector Ufluid(0, 0, 0);        // 颗粒中心处流体速度(可以指定是否使用插值模型计算)
  scalar dragCoefficient(0);     // 颗粒阻力系数
  scalar voidfraction(1);        // 颗粒中心所在网格的空隙率(可以指定是否使用插值模型计算)
  vector Ur(0, 0, 0);            // 相对速度
  scalar Rep(0);                 // 颗粒雷诺数
  scalar Cd(0);                  // 流体阻力系数 Cd = sqr(0.63 + 4.8 / sqrt(Rep))

  // #include "resetVoidfractionInterpolator.H"
  // #include "resetUInterpolator.H"
  UInterpolator_.reset(
    interpolation<vector>::New(propsDict_.lookupOrDefault("UInterpolationType", word("cellPointFace")), U_).ptr()
  );
  voidfractionInterpolator_.reset(
    interpolation<scalar>::New(propsDict_.lookupOrDefault("voidfractionInterpolationType", word("cellPoint")), voidfraction_).ptr()
  );

  #include "setupProbeModel.H"

  for (int index = 0; index <  particleCloud_.numberOfParticles(); index++) {
    cellI = particleCloud_.cellIDs()[index][0];
    if (cellI >= 0) {  // particle Found
      calForceKernel(index, cellI, nufField, rhoField, drag, Ufluid, dragCoefficient, voidfraction, Ur, Rep, Cd);

      Pout << "mixDiFeliceDrag force: " << drag[0] << ", " << drag[1] << ", " << drag[2] << endl;

      // Set value fields and write the probe
      if (probeIt_) {
        #include "setupProbeModelfields.H"
        // Note: for other than ext one could use vValues.append(x)
        // instead of setSize
        vValues.setSize(vValues.size() + 1, drag);  // first entry must the be the force
        vValues.setSize(vValues.size() + 1, Ur);
        sValues.setSize(sValues.size() + 1, Rep);
        sValues.setSize(sValues.size() + 1, Cd);
        sValues.setSize(sValues.size() + 1, voidfraction);
        particleCloud_.probeM().writeProbe(index, sValues, vValues);
      }
    }  // End of cellI >= 0
    // 将每个颗粒的总阻力以及总显式阻力传递到 cfdemCloud 中的 impForces()[index] 和 expForces()[index] 中
    forceSubM(0).partToArray(index, drag, vector::zero, Ufluid, dragCoefficient);
  }  // End of index

  Info << "Setting mixDiFeliceDrag force - done\n" << endl;
}

void mixDiFeliceDrag::setMixForce(const std::vector<double>& dimensionRatios) const {
  if (dimensionRatios.size() == 0) { return setForce(); }
  Info << "Setting mixDiFeliceDrag force..." << endl;

  const volScalarField& nufField = forceSubM(0).nuField();
  const volScalarField& rhoField = forceSubM(0).rhoField();
  label cellI(0);                // 颗粒中心所在网格的索引
  vector drag(0, 0, 0);          // 总阻力 = dragCoefficient * Ur
  vector Ufluid(0, 0, 0);        // 颗粒中心处流体速度(可以指定是否使用插值模型计算)
  scalar dragCoefficient(0);     // 颗粒阻力系数
  scalar voidfraction(1);        // 颗粒中心所在网格的空隙率(可以指定是否使用插值模型计算)
  vector Ur(0, 0, 0);            // 相对速度
  scalar Rep(0);                 // 颗粒雷诺数
  scalar Cd(0);                  // 流体阻力系数 Cd = sqr(0.63 + 4.8 / sqrt(Rep))

  // #include "resetVoidfractionInterpolator.H"
  // #include "resetUInterpolator.H"
  UInterpolator_.reset(
    interpolation<vector>::New(propsDict_.lookupOrDefault("UInterpolationType", word("cellPointFace")), U_).ptr()
  );
  voidfractionInterpolator_.reset(
    interpolation<scalar>::New(propsDict_.lookupOrDefault("voidfractionInterpolationType", word("cellPoint")), voidfraction_).ptr()
  );

  #include "setupProbeModel.H"

  for (int index = 0; index <  particleCloud_.numberOfParticles(); index++) {
    if (particleCloud_.checkFAndMParticle(dimensionRatios[index])) {
      cellI = particleCloud_.cellIDs()[index][0];
      if (cellI >= 0) {  // particle Found
        calForceKernel(index, cellI, nufField, rhoField, drag, Ufluid, dragCoefficient, voidfraction, Ur, Rep, Cd);

        Pout << "mixDiFeliceDrag force: " << drag[0] << ", " << drag[1] << ", " << drag[2] << endl;

        // Set value fields and write the probe
        if (probeIt_) {
          #include "setupProbeModelfields.H"
          // Note: for other than ext one could use vValues.append(x)
          // instead of setSize
          vValues.setSize(vValues.size() + 1, drag);  // first entry must the be the force
          vValues.setSize(vValues.size() + 1, Ur);
          sValues.setSize(sValues.size() + 1, Rep);
          sValues.setSize(sValues.size() + 1, Cd);
          sValues.setSize(sValues.size() + 1, voidfraction);
          particleCloud_.probeM().writeProbe(index, sValues, vValues);
        }
      }  // End of cellI >= 0
      // 将每个颗粒的总阻力以及总显式阻力传递到 cfdemCloud 中的 impForces()[index] 和 expForces()[index] 中
      forceSubM(0).partToArray(index, drag, vector::zero, Ufluid, dragCoefficient);
    }  // End of fine and middle paritcle
  }  // End of index

  Info << "Setting mixDiFeliceDrag force - done\n" << endl;
}

}  // End of namespace Foam
