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
#include "./mixShirgaonkarIB.H"
#include "addToRunTimeSelectionTable.H"
#include "voidFractionModel.H"

namespace Foam {

defineTypeNameAndDebug(mixShirgaonkarIB, 0);

addToRunTimeSelectionTable(forceModel, mixShirgaonkarIB, dictionary);

// @brief 构造函数
mixShirgaonkarIB::mixShirgaonkarIB(const dictionary& dict,
                                   cfdemCloud& sm):
  forceModel(dict,sm),
  propsDict_(dict.subDict(typeName + "Props")),
  twoDimensional_(false),
  depth_(1),
  velFieldName_(propsDict_.lookup("velFieldName")),
  U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
  pressureFieldName_(propsDict_.lookup("pressureFieldName")),
  p_(sm.mesh().lookupObject<volScalarField> (pressureFieldName_)),
  // 从 dict 读取体积分数场的名称
  volumefractionFieldName_(propsDict_.lookup("volumefractionFieldName")),
  // 获取体积分数场的引用
  volumefractions_(sm.mesh().lookupObject<volScalarField>(volumefractionFieldName_)),
  // 从 dict 读取 lmpf 的名称
  lmpfFieldName_(propsDict_.lookup("lmpfFieldName")),
  // 获取 lmpf 的引用
  lmpf_(sm.mesh().lookupObject<volVectorField>(lmpfFieldName_)),
  useTorque_(false) {

  // Append the field names to be probed
  particleCloud_.probeM().initialize(typeName, typeName + ".logDat");
  // First entry must the be the force
  particleCloud_.probeM().vectorFields_.append("dragForce");
  particleCloud_.probeM().writeHeader();

  if (propsDict_.found("twoDimensional")) {
    twoDimensional_ = true;
    depth_ = readScalar(propsDict_.lookup("depth"));
    Info << "2-dimensional simulation - make sure DEM side is 2D" << endl;
    Info << "depth of domain is assumed to be :" << depth_ << endl;
  }

  if (propsDict_.found("useTorque")) {
    useTorque_ = true;
  }

  // Init force sub model
  setForceSubModels(propsDict_);

  // Define switches which can be read from dict
  forceSubM(0).setSwitchesList(0, true);  // activate treatExplicit switch
  forceSubM(0).setSwitchesList(3, true);  // activate search for verbose switch

  // Read those switches defined above, if provided in dict
  for (int iFSub = 0; iFSub < nrForceSubModels(); iFSub++) {
    forceSubM(iFSub).readSwitches();
  }

  particleCloud_.checkCG(false);
}

// @brief Destructor
mixShirgaonkarIB::~mixShirgaonkarIB() {}

void mixShirgaonkarIB::calForceKernel(const int& index,
                                      const volVectorField& IBDrag,
                                      vector& drag,
                                      vector& torque) const {
  drag = vector::zero;
  torque = vector::zero;
  vector sumLmpf_1 = vector::zero;
  vector sumLmpf_2 = vector::zero;
  vector sumLmpf_3 = vector::zero;
  vector sumLmpf_4 = vector::zero;
  vector sumLmpf_5 = vector::zero;
  vector positionCenter = particleCloud_.position(index);

  double vfCell = 0.0, rhoCell = 0.0, vCell = 0.0;

  // 获取颗粒速度
  vector Us = particleCloud_.velocity(index);
  // 获取颗粒半径
  double radius = particleCloud_.radius(index);

  for (int subCell = 0; subCell < particleCloud_.cellsPerParticle()[index][0]; subCell++) {
    label cellI = particleCloud_.cellIDs()[index][subCell];
    if (cellI >= 0) {
      // 计算雷诺数
      double Re = mag(Us) * (2.0 * radius) / forceSubM(0).nuField()[cellI];
      vector cellCenter = particleCloud_.mesh().C()[cellI];
      vfCell = volumefractions_[cellI];
      rhoCell = forceSubM(0).rhoField()[cellI];
      vCell = particleCloud_.mesh().V()[cellI];

      sumLmpf_1 += lmpf_[cellI] * rhoCell * vCell;
      sumLmpf_2 += lmpf_[cellI] * rhoCell * vCell * (1 - vfCell);
      if (vfCell > 0.0 && vfCell < 1.0) {
        sumLmpf_3 += lmpf_[cellI] * rhoCell * vCell;
        sumLmpf_4 += lmpf_[cellI] * rhoCell * vCell * (1 - vfCell);
      } else {
        sumLmpf_5 += lmpf_[cellI] * rhoCell * vCell;
      }
      // 计算阻力与力矩
      drag += IBDrag[cellI] * vCell * (1 - vfCell);
      torque += (cellCenter - positionCenter) ^ IBDrag[cellI] * vCell * (1 - vfCell);
      if (Re < 50) {
        drag -= lmpf_[cellI] * rhoCell * vCell;
        torque -= (cellCenter - positionCenter) ^ lmpf_[cellI] * rhoCell * vCell;
      } else {
        drag -= lmpf_[cellI] * rhoCell * vCell * (1 - vfCell);
        torque -= (cellCenter - positionCenter) ^ lmpf_[cellI] * rhoCell * vCell * (1 - vfCell);
      }
      // drag += (IBDrag[cellI] - lmpf_[cellI] * forceSubM(0).rhoField()[cellI]) * particleCloud_.mesh().V()[cellI] * (1 - volumefractions_[cellI]);
      // torque -= (cellCenter - positionCenter) ^ (IBDrag[cellI] - lmpf_[cellI] * forceSubM(0).rhoField()[cellI]) * particleCloud_.mesh().V()[cellI] * (1 - volumefractions_[cellI]);
    }
  }  // End of loop subCell
  if (twoDimensional_) { drag /= depth_; }

  if (mag(sumLmpf_1) > 0.0) {
    Pout << "sumLmpf_1: " << sumLmpf_1 << endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if (mag(sumLmpf_2) > 0.0) {
    Pout << "sumLmpf_2: " << sumLmpf_2 << endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if (mag(sumLmpf_3) > 0.0) {
    Pout << "sumLmpf_3: " << sumLmpf_3 << endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if (mag(sumLmpf_4) > 0.0) {
    Pout << "sumLmpf_4: " << sumLmpf_4 << endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if (mag(sumLmpf_5) > 0.0) {
    Pout << "sumLmpf_5: " << sumLmpf_5 << endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void mixShirgaonkarIB::setForce() const {
  Info << "Setting mixShirgaonkarIB force..." << endl;

  // IBDrag dimensionSet = (1, -2, -2, 0, 0)
  volVectorField IBDrag = forceSubM(0).IBDragPerV(U_, p_);
  vector drag;
  vector torque;
  #include "setupProbeModel.H"

  for (int index = 0; index < particleCloud_.numberOfParticles(); index++) {
    drag = vector::zero;
    torque = vector::zero;
    calForceKernel(index, IBDrag, drag, torque);

    if (probeIt_) {
      #include "setupProbeModelfields.H"
      // Note: for other than ext one could use vValues.append(x)
      // instead of setSize
      vValues.setSize(vValues.size() + 1, drag);
      particleCloud_.probeM().writeProbe(index, sValues, vValues);
    }
    Info << "mixShirgaonkarIB" << index << ": " << drag[0] << ", " << drag[1] << ", " << drag[2] << endl;
    // write particle based data to global array
    forceSubM(0).partToArray(index, drag, vector::zero);

    if (useTorque_) {
      for (int j = 0; j < 3; j++) {
        particleCloud_.DEMTorques()[index][j] = torque[j];
      }
    }
    if (forceSubM(0).verbose()) {
      Info << "impForces = (" << impForces()[index][0] << ", " << impForces()[index][1] << ", " <<impForces()[index][2] << ")" << endl;
    }
  }  // End of index

  Info << "Setting mixShirgaonkarIB force - done" << endl;
}

void mixShirgaonkarIB::setMixForce(const std::vector<double>& dimensionRatios) const {
  if (dimensionRatios.size() == 0) { return setForce(); }
  Info << "Setting mixShirgaonkarIB force..." << endl;

  volVectorField IBDrag = forceSubM(0).IBDragPerV(U_, p_);
  vector drag;
  vector torque;
  #include "setupProbeModel.H"

  for (int index = 0; index < particleCloud_.numberOfParticles(); index++) {
    if (particleCloud_.needSetFieldForCoarseParticle(index, false, dimensionRatios)) {
      drag = vector::zero;
      torque = vector::zero;
      calForceKernel(index, IBDrag, drag, torque);

      if (probeIt_) {
        #include "setupProbeModelfields.H"
        // Note: for other than ext one could use vValues.append(x)
        // instead of setSize
        vValues.setSize(vValues.size() + 1, drag);
        particleCloud_.probeM().writeProbe(index, sValues, vValues);
      }
      if (mag(drag) > 0.0) {
        Pout << "mixShirgaonkarIB_" << index << ": " << drag[0] << ", " << drag[1] << ", " << drag[2] << endl;
      }
      // write particle based data to global array
      forceSubM(0).partToArray(index, drag, vector::zero);
      if (useTorque_) {
        for (int j = 0; j < 3; j++) {
          particleCloud_.DEMTorques()[index][j] = torque[j];
        }
      }
    }  // End of coarse particle
  }  // End of index

  Info << "Setting mixShirgaonkarIB force - done" << endl;
}

}  // End of namespace Foam
