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

    two way DEM-CFD coupling via MPI

Class
    mixTwoWayMPI

SourceFiles
    mixTwoWayMPI.C

\*---------------------------------------------------------------------------*/

#ifndef mixTwoWayMPI_H
#define mixTwoWayMPI_H

#include "dataExchangeModel.H"
#include "liggghtsCommandModel.H"
#include "OFstream.H"
#include "sys/stat.h"
#include "pair.h"
#include "force.h"
#include "forceModel.H"

// LAMMPS/LIGGGHTS
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lammps.h>  // LAMMPS include files
#include <input.h>
#include <atom.h>
#include <library.h>
#include <error.h>
#include <library_cfd_coupling.h>
#include <update.h>

namespace Foam {

class mixTwoWayMPI: public dataExchangeModel {

private:

  dictionary propsDict_;

  MPI_Comm comm_liggghts;

  inline char* wordToChar(word& inWord) const {
    return const_cast<char*>(inWord.c_str());
  }

protected:

  LAMMPS_NS::LAMMPS* lmp;

public:

  //! @brief Runtime type information
  TypeName("mixTwoWayMPI");

  //! @brief Constructors
  mixTwoWayMPI(const dictionary& dict, cfdemCloud& sm);

  //! @brief Destructor
  ~mixTwoWayMPI();

  void getData(word name,
               word type,
               double ** const& field,
               label step) const;

  void getData(word name,
               word type,
               int ** const& field,
               label step) const;

  void giveData(word name,
                word type,
                double ** const& field,
                const char* datatype) const;

  //! @brief Allocate and destroy for 2-D array
  //! @note DType = double
  virtual void destroy(double*& array) const;
  virtual void destroy(double**& array) const;
  virtual void destroy(double**& array, int) const;
  virtual void allocateArray(double**& array, double initVal, int width, int length) const;
  virtual void allocateArray(double**& array, double initVal, int width, const char* length = "nparticles") const;

  //! @brief Allocate and destroy for 2-D array
  //! @note DType = int
  virtual void destroy(int*& array) const;
  virtual void destroy(int**& array) const;
  virtual void destroy(int**& array, int) const;
  virtual void allocateArray(int**& array, int initVal, int width, int length) const;
  virtual void allocateArray(int**& array, int initVal, int width, const char* length = "nparticles") const;

  //! @brief 释放离散内存
  virtual void destroyDiscreteMemory(double** const& array, int len) const;

  //! @brief 释放离散内存
  virtual void destroyDiscreteMemory(int** const& array, int len) const;

  bool couple(int) const;

  int getNumberOfParticles() const;

  int getNumberOfClumps() const;

  int getNumberOfTypes() const;

  double* getTypeVol() const;

  inline word myType() const { return typeName; }

  inline void setCG() const { setCGs(lmp); }

};  // class mixTwoWayMPI

}  // namespace Foam

#endif  // mixTwoWayMPI_H
