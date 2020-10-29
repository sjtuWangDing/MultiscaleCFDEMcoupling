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

#include "dataExchangeModel.H"

namespace Foam {

defineTypeNameAndDebug(dataExchangeModel, 0);

defineRunTimeSelectionTable(dataExchangeModel, dictionary);

//! @brief Constructor
dataExchangeModel::dataExchangeModel(const dictionary& dict,
                                     cfdemCloud& sm):
  dict_(dict),
  particleCloud_(sm),
  maxNumberOfParticles_(0),
  couplingStep_(0),
  DEMts_(-1.),
  couplingInterval_(readScalar(dict_.lookup("couplingInterval"))),
  // 在初始化 dataExchangeModel 的时候，记录下当前的流体时间步为 timeIndexOffset
  timeIndexOffset_(particleCloud_.mesh().time().timeIndex()) {}

//! @brief Destructor
dataExchangeModel::~dataExchangeModel() {}

//! @brief Allocate and destroy for 2-D array
//! @note DType = double
void dataExchangeModel::destroy(double**& array) const {
  if (array == NULL) { return; }
  delete [] array;
  array = NULL;
}

void dataExchangeModel::destroy(double**& array, int) const {
  if (array == NULL) { return; }
  delete [] array[0];
  delete [] array;
  array = NULL;
}

void dataExchangeModel::allocateArray(double**& array, double initVal,
                                      int width, int length) const {
  // allocate and init array
  destroy(array, -1);
  array = new double*[length];
  double *data = new double[width * length];
  std::fill_n(data, width * length, initVal);

  for (int i = 0; i < length; i++) {
    array[i] = data + i * width;
  }
}

void dataExchangeModel::allocateArray(double**& array, double initVal,
                                      int width, const char* length) const {
  int len = 0;
  if (strcmp(length, "nparticles") == 0) {
    len = particleCloud_.numberOfParticles();
  } else if (strcmp(length,"nbodies") == 0) {
    len = particleCloud_.numberOfClumps();
  } else {
    FatalError << "call allocateArray with length, nparticles or nbodies!\n" << abort(FatalError);
  }
  allocateArray(array, initVal, width, len);
}

//! @brief Allocate and destroy for 2-D array
//! @note DType = int
void dataExchangeModel::destroy(int**& array) const {
  if (array == NULL) { return; }
  delete [] array;
  array = NULL;
}

void dataExchangeModel::destroy(int**& array, int) const {
  if (array == NULL) { return; }
  delete [] array[0];
  delete [] array;
  array = NULL;
}

void dataExchangeModel::allocateArray(int**& array, int initVal,
                                      int width, int length) const {
  // allocate and init array
  destroy(array, -1);
  array = new int*[length];
  int *data = new int[width * length];
  std::fill_n(data, width * length, initVal);

  for (int i = 0; i < length; i++) {
    array[i] = data + i * width;
  }
}

void dataExchangeModel::allocateArray(int**& array, int initVal,
                                      int width, const char* length) const {
  int len = 0;
  if (strcmp(length, "nparticles") == 0) {
    len = particleCloud_.numberOfParticles();
  } else if (strcmp(length,"nbodies") == 0) {
    len = particleCloud_.numberOfClumps();
  } else {
    FatalError << "call allocateArray with length, nparticles or nbodies!\n" << abort(FatalError);
  }
  allocateArray(array, initVal, width, len);
}

//! @brief Allocate and destroy for 1-D array
//! @note DType = double
void dataExchangeModel::destroy(double*& array) const {
  delete [] array;
  array = NULL;
}

void dataExchangeModel::allocateArray(double*& array, double initVal, int length) const {
  destroy(array);
  // allocate and init array
  array = new double[length];
  std::fill_n(array, length, initVal);
}

//! @brief Allocate and destroy for 1-D array
//! @note DType = int
void dataExchangeModel::destroy(int*& array) const {
  delete [] array;
  array = NULL;
}

void dataExchangeModel::allocateArray(int*& array, int initVal, int length) const {
  destroy(array);
  // allocate and init array
  array = new int[length];
  std::fill_n(array, length, initVal);
}

bool dataExchangeModel::couple(int i) const {
  bool coupleNow = false;
  if (doCoupleNow()) {
    couplingStep_ += 1;
    coupleNow = true;
  }
  return coupleNow;
}

}  // End of namespace Foam
