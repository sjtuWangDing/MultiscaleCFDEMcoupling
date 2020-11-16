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
  Foam::twoWayMPI
\*---------------------------------------------------------------------------*/

#include "twoWayMPI.H"
#include "subModels/liggghtsCommandModel/liggghtsCommandModel.H"

namespace Foam {

cfdemDefineTypeName(twoWayMPI)

cfdemAddToNewFunctionMap(dataExchangeModel, twoWayMPI)

twoWayMPI::twoWayMPI(cfdemCloud& cloud):
  dataExchangeModel(cloud),
  lmp_(nullptr),
  subPropsDict_(cloud.couplingPropertiesDict().subDict(typeName_ + "Props")) {

  Info << "Starting up LIGGGHTS for first time execution."<<endl;
  MPI_Comm_dup(MPI_COMM_WORLD, &lgs_comm_);

  // open LIGGGHTS input script
  const fileName lgsPath(subPropsDict_.lookup("liggghtsPath"));
  lmp_ = new LAMMPS_NS::LAMMPS(0, nullptr, lgs_comm_);
  lmp_->input->file(lgsPath.c_str());

  // get DEM time step size
  DEMts_ = lmp_->update->dt;
  checkTimeStepSize();
}

//! \return 当前耦合时间步中颗粒的数量
int twoWayMPI::couple() {
  Info << "dataExchangeModel " << typeName() << ": Starting up CFD-DEM couping" << endl;
  Info << "Starting up LIGGGHTS" << endl;
  couplingStep_ += 1;
  for (const std::shared_ptr<liggghtsCommandModel>& model: cloud_.liggghtsCommandModels()) {
    liggghtsCommandModel* lgs = model.get();
    // 在 runCommand 中会判断当前的 couplingStep_ 是否满足时间步的要求
    if (lgs->runCommand(couplingStep_)) {
      int commandLines = lgs->commandLines();
      for (int j = 0; j < commandLines; ++j) {
        std::string command = lgs->getCommand(j);
        Info << "Executing command: " << command << endl;
        lmp_->input->one(command.c_str());
      }
    }
  }
  Info << "LIGGGHTS finished" << endl;
  // 返回颗粒数量
  return liggghts_get_maxtag(lmp_);
}

//! \brief Allocate and destroy for 2-D double array
void twoWayMPI::destroy(double**& array) {
}
void twoWayMPI::destroy(double**& array, int) {

}
void twoWayMPI::allocateArray(double**& array, double initVal, int width, int length) {

}
void twoWayMPI::allocateArray(double**& array, double initVal, int width, const char* length/* = "nparticles" */) {

}

//! \brief Allocate and destroy for 2-D int array
void twoWayMPI::destroy(int**& array) {

}
void twoWayMPI::destroy(int**& array, int) {

}
void twoWayMPI::allocateArray(int**& array, int initVal, int width, int length) {

}
void twoWayMPI::allocateArray(int**& array, int initVal, int width, const char* length/* = "nparticles" */) {

}

//! \brief Allocate and destroy for 1-D double array
void twoWayMPI::destroy(double*& array) {

}
void twoWayMPI::allocateArray(double*& array, double initVal, int length) {

}

//! \brief Allocate and destroy for 1-D int array
void twoWayMPI::destroy(int*& array) {

}
void twoWayMPI::allocateArray(int*& array, int initVal, int length) {

}

void twoWayMPI::getData(const std::string& dataName,
              const std::string& dataType,
              double** const& field,
              label step) {

}
void twoWayMPI::getData(const std::string& dataName,
              const std::string& dataType,
              int** const& field,
              label step) {

}
void twoWayMPI::giveData(const std::string& dataName,
              const std::string& dataType,
              double** const& field,
              const char* fieldType/* = "double" */) {

}

} // namespace Foam
