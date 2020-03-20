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
#include "mixArchimedesIB.H"
#include "addToRunTimeSelectionTable.H"
#include "voidFractionModel.H"

namespace Foam {

defineTypeNameAndDebug(mixArchimedesIB, 0);

addToRunTimeSelectionTable(forceModel, mixArchimedesIB, dictionary);

// @brief Constructors
mixArchimedesIB::mixArchimedesIB(const dictionary& dict,
                                 cfdemCloud& sm):
  forceModel(dict, sm),
  propsDict_(dict.subDict(typeName + "Props")),
  twoDimensional_(false),
  // 从 dict 读取体积分数场的名称
  volumefractionFieldName_(propsDict_.lookup("volmefracionFieldName")),
  // 获取体积分数场的引用
  volumefractions_(sm.mesh().lookupObject<volScalarField> (volumefractionFieldName_)),
  gravityFieldName_(propsDict_.lookup("gravityFieldName")),
#if defined(version21) || defined(version16ext)
  g_(sm.mesh().lookupObject<uniformDimensionedVectorField>(gravityFieldName_))
#elif  defined(version15)
  g_(dimensionedVector(sm.mesh().lookupObject<IOdictionary>("environmentalProperties").lookup(gravityFieldName_)).value())
#endif
{
  // Append the field names to be probed
  particleCloud_.probeM().initialize(typeName, typeName + ".logDat");
  // First entry must the be the force
  particleCloud_.probeM().vectorFields_.append("mixArchimedesIBForce");
  particleCloud_.probeM().writeHeader();

  if (propsDict_.found("twoDimensional")) {
    twoDimensional_ = true;
    Info << "2-dimensional simulation - make sure DEM side is 2D" << endl;
  }

  // Init force sub model
  setForceSubModels(propsDict_);

  // 激活搜索 switch treatForceExplicit
  forceSubM(0).setSwitchesList(0, true);

  // 强制设置 treatForceExplicit = true
  forceSubM(0).setSwitches(0, true);

  // 从 dict 读入需要搜索的 switches
  forceSubM(0).readSwitches();

  // 强制设置 switch treatForceDEM = true
  // 仅在 DEM 中考虑 Archimedes, 即仅考虑流体对颗粒的作用, 不考虑颗粒对流体的反作用
  forceSubM(0).setSwitches(1, true);
  Info << "accounting for Archimedes only on DEM side!" << endl;

  particleCloud_.checkCG(true);
}

// @brief Destructor
mixArchimedesIB::~mixArchimedesIB() {}

void mixArchimedesIB::setForce() const {
  FatalError << "mixArchimedesIB.H: setForce() only used to cover pure virtual function in base class"
    << abort(FatalError);
}

void mixArchimedesIB::setForce(const std::vector<double>& dimensionRatios) const {

  Info << "Setting mixArchimedesIB force......" << endl;

  #include "setupProbeModel.H"

  for (int index = 0; index <  particleCloud_.numberOfParticles(); ++index) {

    if (particleCloud_.checkCoarseParticle(dimensionRatios[index])) {

      vector force = vector::zero;

      for(int subCell = 0; subCell < particleCloud_.cellsPerParticle()[index][0]; subCell++) {

        label cellI = particleCloud_.cellIDs()[index][subCell];
        if (cellI >= 0) {  // cell found
          force += -g_.value() * forceSubM(0).rhoField()[cellI] * particleCloud_.mesh().V()[cellI] * (1 - volumefractions_[cellI]);
        }  // End of cell Found
      }  // End of loop all subCell

      if (probeIt_) {
        #include "setupProbeModelfields.H"
        // Note: for other than ext one could use vValues.append(x)
        // instead of setSize
        vValues.setSize(vValues.size() + 1, force);
        particleCloud_.probeM().writeProbe(index, sValues, vValues);
      }

      // set force on particle
      if (twoDimensional_) {
        Warning << "mixArchimedesIB model doesn't work for 2D right now!!\n" << endl;
      }

      // write particle based data to global array
      forceSubM(0).partToArray(index, force, vector::zero);

    }  // End of check coarse particles
  }  // End of loop all particles
}

}  // End of namespcae Foam
