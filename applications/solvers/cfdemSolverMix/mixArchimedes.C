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
#include "mixArchimedes.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam {

defineTypeNameAndDebug(mixArchimedes, 0);

addToRunTimeSelectionTable(forceModel, mixArchimedes, dictionary);

// \brief Constructors
mixArchimedes::mixArchimedes(const dictionary& dict, cfdemCloud& sm):
  forceModel(dict, sm),
  propsDict_(dict.subDict(typeName + "Props")),
  twoDimensional_(false),
  gravityFieldName_(propsDict_.lookup("gravityFieldName")),
#if defined(version21) || defined(version16ext)
  g_(sm.mesh().lookupObject<uniformDimensionedVectorField>(gravityFieldName_))
#elif defined(version15)
  g_(dimensionedVector(sm.mesh().lookupObject<IOdictionary>("environmentalProperties").lookup(gravityFieldName_)).value())
#endif
{
  // Append the field names to be probed
  // suppress particle probe
  if (probeIt_ && propsDict_.found("suppressProbe")) {
    probeIt_ = !Switch(propsDict_.lookup("suppressProbe"));
  }
  if (probeIt_) {
    particleCloud_.probeM().initialize(typeName, typeName + ".logDat");
    particleCloud_.probeM().vectorFields_.append("archimedesForce");
    particleCloud_.probeM().scalarFields_.append("Vp");
    particleCloud_.probeM().writeHeader(); 
  }

  if (propsDict_.found("twoDimensional")) {
    twoDimensional_ = true;
    Info << "2-dimensional simulation - make sure DEM side is 2D" << endl;
  }

  // init force sub model
  setForceSubModels(propsDict_);

  // define switches which can be read from dict (default = false)
  forceSubM(0).setSwitchesList(3, true);  // activate search for verbose switch

    // formally, according to Zhou et al. B would have treatForceDEM=false
    // , but we do not have g-body force term in NSE --> so treatForceDEM=true always?
    //if (modelType_=="A" || modelType_=="Bfull")
    //    forceSubM(0).setSwitches(1,true);       // Archimedes only on DEM side (treatForceDEM=true)
    //else if (modelType_=="B")
    //    forceSubM(0).setSwitches(1,false);      // Archimedes on CFD and DEM side (treatForceDEM=false)

  // 强制设置在 mixArchimedes 模型中, treatForceDEM = true
  // 即仅在 DEM 中考虑耦合力，即仅考虑流体对颗粒的作用，不考虑颗粒对流体的反作用
  forceSubM(0).setSwitches(1, true); // Archimedes only on DEM side (treatForceDEM = true)

  // for (int iFSub = 0; iFSub < nrForceSubModels(); iFSub++) {
  //   forceSubM(iFSub).readSwitches();
  // }
  forceSubM(0).readSwitches();
}

// \brief Destructor
mixArchimedes::~mixArchimedes() {}

void mixArchimedes::calForceKernel(const int& index,
                                   const int& cellI,
                                   vector& force) const {
  if (twoDimensional_) {
    scalar r = particleCloud_.radius(index);
    force = -g_.value() * forceSubM(0).rhoField()[cellI] * r * r * M_PI;
    Warning << "mixArchimedes::calForceKernel() : this functionality is not tested!" << endl;
  } else {
    // 真实颗粒直径
    scalar dReal = particleCloud_.d(index);
    // 定义修正直径
    scalar dScale = particleCloud_.d(index);
    // scale diameter
    forceSubM(0).scaleDia(dScale, index);
    // 颗粒体积
    scalar Vs = dScale * dScale * dScale * M_PI / 6.0;
    // 计算 Archimedes force
    force = -g_.value() * forceSubM(0).rhoField()[cellI] * Vs;
    forceSubM(0).scaleForce(force, dReal);
    // Pout << "-g_.value(): " << -g_.value() << endl;
    // Pout << "dScale: " << dScale << endl;
    // Pout << "Vs: " << Vs << endl;
    // Pout << "forceSubM(0).rhoField()[" << cellI << "]: " << forceSubM(0).rhoField()[cellI] << endl;
  }
}

void mixArchimedes::setForce() const {
  Info << "Setting mixArchimedes force..." << endl;

  vector force(0, 0, 0);
  #include "setupProbeModel.H"

  for (int index = 0; index < particleCloud_.numberOfParticles(); ++index) {
    label cellI = particleCloud_.cellIDs()[index][0];
    force = vector::zero;
    if (cellI >= 0) {
      calForceKernel(index, cellI, force);

      Pout << "mixArchimedes force: " << force[0] << ", " << force[1] << ", " << force[2] << endl;

      // Set value fields and write the probe
      if (probeIt_) {
        #include "setupProbeModelfields.H"
        // Note: for other than ext one could use vValues.append(x) instead of setSize
        vValues.setSize(vValues.size() + 1, force);  // first entry must the be the force
        sValues.setSize(sValues.size() + 1, particleCloud_.particleVolume(index)); 
        particleCloud_.probeM().writeProbe(index, sValues, vValues);
      }
    }  // End of cellI >= 0
    // 将每个颗粒的总阻力以及总显式阻力传递到 cfdemCloud 中的 impForces()[index] 和 expForces()[index] 中
    forceSubM(0).partToArray(index, force, vector::zero);
  }  // End of index

  Info << "Setting mixArchimedes force - done" << endl;
}

void mixArchimedes::setMixForce(const std::vector<double>& dimensionRatios) const {
  if (dimensionRatios.size() == 0) { return setForce(); }
  Info << "Setting mixArchimedes force..." << endl;

  vector force(0, 0, 0);
  #include "setupProbeModel.H"

  for (int index = 0; index < particleCloud_.numberOfParticles(); ++index) {
    force = vector::zero;
    if (particleCloud_.checkFAndMParticle(dimensionRatios[index])) {
      label cellI = particleCloud_.cellIDs()[index][0];
      if (cellI >= 0) {
        calForceKernel(index, cellI, force);

        Pout << "mixArchimedes force: " << force[0] << ", " << force[1] << ", " << force[2] << endl;

        // Set value fields and write the probe
        if (probeIt_) {
          #include "setupProbeModelfields.H"
          // Note: for other than ext one could use vValues.append(x) instead of setSize
          vValues.setSize(vValues.size() + 1, force);  // first entry must the be the force
          sValues.setSize(sValues.size() + 1, particleCloud_.particleVolume(index)); 
          particleCloud_.probeM().writeProbe(index, sValues, vValues);
        }
      }  // End of cellI >= 0
      // 将每个颗粒的总阻力以及总显式阻力传递到 cfdemCloud 中的 impForces()[index] 和 expForces()[index] 中
      forceSubM(0).partToArray(index, force, vector::zero);
    }  // End of fine and middle particle
  }  // End of index

  Info << "Setting mixArchimedes force - done" << endl;
}

}  // End of namespace Foam
