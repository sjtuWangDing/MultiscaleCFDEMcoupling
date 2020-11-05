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
#include "mixTwoWayMPI.H"
#include "addToRunTimeSelectionTable.H"
#include "clockModel.H"
#include "pair.h"
#include "force.h"
#include "forceModel.H"

namespace Foam {

defineTypeNameAndDebug(mixTwoWayMPI, 0);

addToRunTimeSelectionTable(dataExchangeModel, mixTwoWayMPI, dictionary);

// @brief Constructors
mixTwoWayMPI::mixTwoWayMPI(const dictionary& dict,
                           cfdemCloud& sm):
  dataExchangeModel(dict, sm),
  propsDict_(dict.subDict(typeName + "Props")),
  lmp(NULL) {

  Info << "mixTwoWayMPI::mixTwoWayMPI(): Starting up LIGGGHTS for first time execution\n" << endl;

  MPI_Comm_dup(MPI_COMM_WORLD, &comm_liggghts);

  // read path from dictionary
  // fileName 为 liggghts 的执行脚本
  const fileName liggghtsPath(propsDict_.lookup("liggghtsPath"));

  // 开始执行 liggghts 脚本
  Info << "Executing liggghts input script "<< liggghtsPath.c_str() << "" << endl;
  lmp = new LAMMPS_NS::LAMMPS(0, NULL, comm_liggghts);
  lmp->input->file(liggghtsPath.c_str());

  // get DEM time step size
  DEMts_ = lmp->update->dt;
  checkTSsize();
}

// @brief Destructor
mixTwoWayMPI::~mixTwoWayMPI() {
  delete lmp;
}

void mixTwoWayMPI::getData(word name,
                           word type,
                           double** const& field,
                           label step) const {
  char* charName = wordToChar(name);
  char* charType = wordToChar(type);
  data_liggghts_to_of(charName, charType, lmp, (void*&)field, (char*)"double");
}

void mixTwoWayMPI::getData(word name,
                           word type,
                           int** const& field,
                           label step) const {
  char* charName = wordToChar(name);
  char* charType = wordToChar(type);
  data_liggghts_to_of(charName, charType, lmp, (void*&)field, (char*)"int");
}

void mixTwoWayMPI::giveData(word name,
                            word type,
                            double ** const& field,
                            const char* datatype) const {
  char* charName = wordToChar(name);
  char* charType = wordToChar(type);
  char* charDatatype = const_cast<char*>(datatype);
  data_of_to_liggghts(charName, charType, lmp, (void*)field, charDatatype);
}

// @brief Allocate and destroy for 2-D array
// @note DType = double
void mixTwoWayMPI::destroy(double**& array) const {
  if (array == NULL) { return; }
  free(array);
  array = NULL;
}

void mixTwoWayMPI::destroy(double**& array, int) const {
  if (array == NULL) { return; }
  free(array[0]);
  free(array);
  array = NULL;
}

// @brief Allocate and destroy for 1-D array
// @note DType = double
void mixTwoWayMPI::destroy(double*& array) const {
  if (array == NULL) { return; }
  free(array);
  array = NULL;
}

void mixTwoWayMPI::allocateArray(double**& array,
                                 double initVal,
                                 int width,
                                 int length) const {
  allocate_external_double(array, width, length, initVal, lmp);
}

void mixTwoWayMPI::allocateArray(double**& array,
                                 double initVal,
                                 int width,
                                 const char* length) const {
  char* charLength = const_cast<char*>(length);
  allocate_external_double(array, width, charLength, initVal, lmp);
}

// @brief Allocate and destroy for 2-D array
// @note DType = int
void mixTwoWayMPI::destroy(int**& array) const {
  if (array == NULL) { return; }
  free(array);
  array = NULL;
}

void mixTwoWayMPI::destroy(int**& array, int) const {
  if (array == NULL) { return; }
  free(array[0]);
  free(array);
  array = NULL;
}

// @brief Allocate and destroy for 1-D array
// @note DType = int
void mixTwoWayMPI::destroy(int*& array) const {
  if (array == NULL) { return; }
  free(array);
  array = NULL;
}

// @brief 释放离散内存
void mixTwoWayMPI::destroyDiscreteMemory(double** const& array, int len) const {
  if (array == NULL || len == 0) { return; }
  double**& ptr = const_cast<double**&>(array);
  for (int i = 0; i < len; ++i) {
    if (array[i] != NULL) {
      free(ptr[i]);
      ptr[i] = NULL;
    }
  }
  free(ptr);
  ptr = NULL;
}

// @brief 释放离散内存
void mixTwoWayMPI::destroyDiscreteMemory(int** const& array, int len) const {
  if (array == NULL || len == 0) { return; }
  int**& ptr = const_cast<int**&>(array);
  for (int i = 0; i < len; ++i) {
    if (array[i] != NULL) {
      free(ptr[i]);
      ptr[i] = NULL;
    }
  }
  free(ptr);
  ptr = NULL;
}

void mixTwoWayMPI::allocateArray(int**& array,
                                 int initVal,
                                 int width,
                                 int length) const {
  allocate_external_int(array, width, length, initVal, lmp);
}

void mixTwoWayMPI::allocateArray(int**& array,
                                 int initVal,
                                 int width,
                                 const char* length) const {
  char* charLength = const_cast<char*>(length);
  allocate_external_int(array, width, charLength, initVal, lmp);
}

bool mixTwoWayMPI::couple(int i) const {
  bool coupleNow = false;

  if (i == 0) {  // 如果 i == 0, 则表示当前时间步是耦合时间步
    couplingStep_ += 1;
    coupleNow = true;

    // 开始执行 liggghts 脚本
    Info << "mixTwoWayMPI::couple(): Starting up LIGGGHTS..." << endl;

    // check if liggghtsCommandModels with exaxt timing are being run
    bool exactTiming(false);
    int runComNr = -10;
    DynamicList<scalar> interruptTimes(0);
    DynamicList<int> DEMstepsToInterrupt(0);
    DynamicList<int> lcModel(0);

    forAll (particleCloud_.liggghtsCommandModelList(), i) {
      // Check if exact timing is needed
      // get time for execution
      // store time for execution in list
      if (particleCloud_.liggghtsCommand()[i]().exactTiming()) {
        exactTiming = true;
        DynamicList<scalar> h = particleCloud_.liggghtsCommand()[i]().executionsWithinPeriod(TSstart(),TSend());

        forAll (h, j) {
          // save interrupt times (is this necessary)
          interruptTimes.append(h[j]);

          // calc stepsToInterrupt
          DEMstepsToInterrupt.append(DEMstepsTillT(h[j]));

          // remember which liggghtsCommandModel to run
          lcModel.append(i);
        }

        // make cumulative
        label len = DEMstepsToInterrupt.size();
        label ind(0);
        forAll (DEMstepsToInterrupt, i) {
          ind = len - i - 1;
          if (ind > 0) {
            DEMstepsToInterrupt[ind] -= DEMstepsToInterrupt[ind - 1];
          }
        }

        Info << "Foam::twoWayMPI::couple(i): interruptTimes = " << interruptTimes << endl;
        Info << "Foam::twoWayMPI::couple(i): DEMstepsToInterrupt = " << DEMstepsToInterrupt << endl;
        Info << "Foam::twoWayMPI::couple(i): lcModel = " << lcModel << endl;
      }

      if (particleCloud_.liggghtsCommand()[i]().type() == "runLiggghts") {
        runComNr=i;
      }
    }  // End of loop liggghtsCommandModelList

    // models with exact timing exists
    label commandLines(0);

    if (exactTiming) {
      // extension for more liggghtsCommands active the same time:
      //    sort interrupt list within this run period
      //    keep track of corresponding liggghtsCommand
      int DEMstepsRun = 0;

      forAll (interruptTimes, j) {
        // set run command till interrupt
        DEMstepsRun += DEMstepsToInterrupt[j];
        particleCloud_.liggghtsCommand()[runComNr]().set(DEMstepsToInterrupt[j]);
        const char* command = particleCloud_.liggghtsCommand()[runComNr]().command(0);
        Info << "Executing run command: " << command << endl;
        lmp->input->one(command);

        // run liggghts command with exact timing
        command = particleCloud_.liggghtsCommand()[lcModel[j]]().command(0);
        Info << "Executing command: " << command << endl;
        lmp->input->one(command);
      }

      // do the run
      if (particleCloud_.liggghtsCommand()[runComNr]().runCommand(couplingStep())) {
        particleCloud_.liggghtsCommand()[runComNr]().set(couplingInterval() - DEMstepsRun);
        const char* command = particleCloud_.liggghtsCommand()[runComNr]().command(0);
        Info << "Executing run command: " << command << endl;
        lmp->input->one(command);
      }

      // do the other non exact timing models
      forAll (particleCloud_.liggghtsCommandModelList(), i) {
        if(!particleCloud_.liggghtsCommand()[i]().exactTiming() &&
            particleCloud_.liggghtsCommand()[i]().runCommand(couplingStep())) {
          commandLines = particleCloud_.liggghtsCommand()[i]().commandLines();
          for(int j = 0; j < commandLines; j++) {
            const char* command = particleCloud_.liggghtsCommand()[i]().command(j);
            Info << "Executing command: " << command << endl;
            lmp->input->one(command);
          }
        }
      }
    } else {
      forAll (particleCloud_.liggghtsCommandModelList(), i) {
        if (particleCloud_.liggghtsCommand()[i]().runCommand(couplingStep())) {
          commandLines = particleCloud_.liggghtsCommand()[i]().commandLines();
          for (int j = 0; j < commandLines; j++) {
            const char* command = particleCloud_.liggghtsCommand()[i]().command(j);
            Info << "Executing command: " << command << endl;
            lmp->input->one(command);
          }
        }
      }
    }  // End of exactTiming
    Info << "mixTwoWayMPI::couple(): LIGGGHTS - done\n" << endl;

    // give number of particles to cloud
    double newNpart = liggghts_get_maxtag(lmp);
    setNumberOfParticles(newNpart);
    setNumberOfClumps(-2);

    // re-allocate arrays of cloud
    Info << "\nmixTwoWayMPI::couple(0): invoke particleCloud_.reAllocArrays()..." << endl;
    particleCloud_.reAllocArrays();
    Info << "mixTwoWayMPI::couple(0): particleCloud_.reAllocArrays() - done\n" << endl;
  }  // End of i == 0
  return coupleNow;
}

int mixTwoWayMPI::getNumberOfParticles() const {
  return liggghts_get_maxtag(lmp);
}

int mixTwoWayMPI::getNumberOfClumps() const {
#ifdef multisphere
  return liggghts_get_maxtag_ms(lmp);
#endif
  Warning << "liggghts_get_maxtag_ms(lmp) is not available here!" << endl;
  return -1;
}

int mixTwoWayMPI::getNumberOfTypes() const {
#ifdef multisphere
  return liggghts_get_ntypes_ms(lmp);
#endif
  Warning << "liggghts_get_maxtag_ms(lmp) is not available here!" << endl;
return -1;
}

double* mixTwoWayMPI::getTypeVol() const {
#ifdef multisphere
  return liggghts_get_vclump_ms(lmp);
#endif
  Warning << "liggghts_get_vclump_ms(lmp) is not available here!" << endl;
  return NULL;
}

}  // End of namespace Foam
