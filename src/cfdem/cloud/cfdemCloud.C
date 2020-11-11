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
#include "subModels/dataExchangeModel/dataExchangeModel.H"
#include "subModels/liggghtsCommandModel/liggghtsCommandModel.H"
#include "subModels/averagingModel/averagingModel.H"

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
  pCloud_(0),
  dataExchangeModel_(dataExchangeModel::New(*this, couplingPropertiesDict_)),
  averagingModel_(averagingModel::New(*this, couplingPropertiesDict_)) {

  Info << "\nEnding of Constructing cfdemCloud Base Class Object......\n" << endl;
  Info << "\nEntry of cfdemCloud::cfdemCloud(const fvMesh&)......\n" << endl;

  for (int i = 0; i < liggghtsCommandModelList().size(); ++i) {
    // liggghtsCommandModel::New() 函数返回的是 std::unique_ptr
    liggghtsCommandModels_.emplace_back(
      liggghtsCommandModel::New(*this, liggghtsCommandsDict_, liggghtsCommandModelList()[i]));
  }
}

/*!
 * \brief 更新函数
 * \note used for cfdemSolverPiso
 * \param VoidF  <[in, out] 小颗粒空隙率场
 * \param Us     <[in, out] 局部平均小颗粒速度场
 * \param U      <[in] 流体速度场
 */
bool cfdemCloud::evolve(volScalarField& VoidF,
                        volVectorField& Us,
                        volVectorField& U) {
  Info << "\nFoam::cfdemCloud::evolve(), used for cfdemSolverPiso......\n" << endl;
  if (!ignore()) {
    if (!writeTimePassed_ && mesh_.time().outputTime()) {
      writeTimePassed_ = true;
    }
    if (dataExchangeM().doCoupleNow()) {
      Info << "evolve coupling..." << endl;
      // couple() 函数执行 liggghts 脚本，并获取新的颗粒数量
      pCloud_.setNumberOfParticles(dataExchangeM().couple());
      Info << "get number of particles: " << pCloud_.numberOfParticles()<< " at coupling step: "
        << dataExchangeM().couplingStep() << endl;

      // 重置局部平均颗粒速度
      averagingM().resetVectorAverage(averagingM().UsPrev(),
                                      averagingM().UsNext(), false);
    }
  }
  Info << "Foam::cfdemCloud::evolve() - done\n" << endl;
}

const std::vector<std::shared_ptr<liggghtsCommandModel>>& cfdemCloud::liggghtsCommandModels() const {
  return liggghtsCommandModels_;
}

const std::vector<std::shared_ptr<forceModel>>& cfdemCloud::forceModels() const {
  return forceModels_;
}

const std::vector<std::shared_ptr<momCoupleModel>>& cfdemCloud::momCoupleModels() const {
  return momCoupleModels_;
}

const dataExchangeModel& cfdemCloud::dataExchangeM() const {
  // Foam::autoPtr<T> 中定义了 inline operator const T&() const;
  return dataExchangeModel_;
}

dataExchangeModel& cfdemCloud::dataExchangeM() {
  return const_cast<dataExchangeModel&>(static_cast<const dataExchangeModel&>(dataExchangeModel_));
}

const averagingModel& cfdemCloud::averagingM() const {
  return averagingModel_;
}

averagingModel& cfdemCloud::averagingM() {
  return const_cast<averagingModel&>(static_cast<const averagingModel&>(averagingModel_));
}

} // namespace Foam
