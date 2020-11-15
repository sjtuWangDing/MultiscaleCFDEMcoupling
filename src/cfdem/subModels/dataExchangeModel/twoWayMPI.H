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

#ifndef __TWO_WAY_MPI_H__
#define __TWO_WAY_MPI_H__

#include "dataExchangeModel.H"

#include "mpi.h"
#include <lammps.h> // LAMMPS include files
#include <input.h>
#include <atom.h>
#include <library.h>
#include <error.h>
#include <library_cfd_coupling.h>
#include <update.h>

namespace Foam {

class twoWayMPI: public dataExchangeModel {

public:

  //! \brief Runtime type information
  cfdemTypeName("twoWayMPI")

  //! \brief Constructor
  twoWayMPI(cfdemCloud& cloud);

  //! \brief Destructor
  ~twoWayMPI();

public:

  //! \return 当前耦合时间步中颗粒的数量
  int couple();

  //! \brief Allocate and destroy for 2-D double array
  void destroy(double**& array);
  void destroy(double**& array, int);
  void allocateArray(double**& array, double initVal, int width, int length);
  void allocateArray(double**& array, double initVal, int width, const char* length = "nparticles");

  //! \brief Allocate and destroy for 2-D int array
  void destroy(int**& array);
  void destroy(int**& array, int);
  void allocateArray(int**& array, int initVal, int width, int length);
  void allocateArray(int**& array, int initVal, int width, const char* length  = "nparticles");

  //! \brief Allocate and destroy for 1-D double array
  void destroy(double*& array);
  void allocateArray(double*& array, double initVal, int length);

  //! \brief Allocate and destroy for 1-D int array
  void destroy(int*& array);
  void allocateArray(int*& array, int initVal, int length);

  void getData(const std::string& dataName,
               const std::string& dataType,
               double** const& field,
               label step);
  void getData(const std::string& dataName,
               const std::string& dataType,
               int** const& field,
               label step);
  void giveData(const std::string& dataName,
                const std::string& dataType,
                double** const& field,
                const char* fieldType = "double");

private:

  LAMMPS_NS::LAMMPS* lmp_;

  MPI_Comm lgs_comm_;

  dictionary subPropsDict_;
};

} // namespace Foam

#endif // __TWO_WAY_MPI_H__
