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

Class
    dataExchangeModel

SourceFiles
    dataExchangeModel.C

\*---------------------------------------------------------------------------*/

#ifndef dataExchangeModel_H
#define dataExchangeModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"

#include "mpi.h"
#include "math.h"
#include <lammps.h>
#include <pair.h>
#include <force.h>
#include <atom.h>

namespace Foam {

class dataExchangeModel {

protected:

  const dictionary& dict_;

  cfdemCloud& particleCloud_;

  int maxNumberOfParticles_;

  //! @brief 耦合计数
  mutable int couplingStep_;

  //! @brief DEM 时间步长
  scalar DEMts_;

  //! @brief 耦合间隔(单位：DEM 时间步)
  int couplingInterval_;

  const int timeIndexOffset_;

public:

  //! @brief Runtime type information
  TypeName("dataExchangeModel");

  //! @brief Declare runtime constructor selection table
  declareRunTimeSelectionTable(
    autoPtr,
    dataExchangeModel,
    dictionary,
    (
      const dictionary& dict,
      cfdemCloud& sm
    ),
    (dict, sm)
  );

  //! @brief Selector
  static autoPtr<dataExchangeModel> New(const dictionary& dict,
                                        cfdemCloud& sm);

  //! @brief Constructor
  dataExchangeModel(const dictionary& dict,
                    cfdemCloud& sm);

  //! @brief Destructor
  virtual ~dataExchangeModel();

  void setNumberOfParticles(int numberOfParticles) const {
    particleCloud_.setNumberOfParticles(numberOfParticles);
  }

  void setNumberOfClumps(int numberOfClumps) const {
    particleCloud_.setNumberOfClumps(numberOfClumps);
  }

  const int& maxNumberOfParticles() const {
    return maxNumberOfParticles_;
  }

  virtual void getData(word name, word type, double** const& field, label step) const = 0;

  virtual void getData(word name, word type, int** const& field, label step) const = 0;

  template <typename DType>
  void getData(word name, word type, DType** const& field) const {
    getData(name, type, field, couplingStep_ - 1);
  }

  virtual void giveData(word name, word type, double** const& field,
                        const char* datatype = "double") const = 0;

  //! @brief Allocate and destroy for 2-D array
  // @note DType = double
  virtual void destroy(double**& array) const;
  virtual void destroy(double**& array, int) const;
  virtual void allocateArray(double**& array, double initVal, int width, int length) const;
  virtual void allocateArray(double**& array, double initVal, int width, const char* length = "nparticles") const;

  //! @brief Allocate and destroy for 2-D array
  // @note DType = int
  virtual void destroy(int**& array) const;
  virtual void destroy(int**& array, int) const;
  virtual void allocateArray(int**& array, int initVal, int width, int length) const;
  virtual void allocateArray(int**& array, int initVal, int width, const char* length  = "nparticles") const;

  //! @brief Allocate and destroy for 1-D array
  // @note DType = double
  virtual void destroy(double*& array) const;
  virtual void allocateArray(double*& array, double initVal, int length) const;

  //! @brief Allocate and destroy for 1-D array
  // @note DType = int
  virtual void destroy(int*& array) const;
  virtual void allocateArray(int*& array, int initVal, int length) const;

  //! @brief 释放离散内存
  virtual void destroyDiscreteMemory(double** const& array, int len) const {
    FatalError << "dataExchangeModel::destroyDiscreteMemory(): using base class function, please use derived class function\n"
      << abort(FatalError);
  }

  //! @brief 释放离散内存
  virtual void destroyDiscreteMemory(int** const& array, int len) const {
    FatalError << "dataExchangeModel::destroyDiscreteMemory(): using base class function, please use derived class function\n"
      << abort(FatalError);
  }

  virtual bool couple(int) const;

  //! @brief 因为耦合时间步长 = 流体时间步长的整数倍，所以 timeStepFraction() 用于计算每个流体时间步在耦合时间步中的所占比例，
  // Eg: 如果 couplingTime() == 3 * CFDts，那么每一个耦合时间步由 3 个流体时间步构成，那么这三个流体时间步的 timeStepFraction() 分别返回 0, 0.333333, 0.666666
  scalar timeStepFraction() const {
    auto CFDts = particleCloud_.mesh().time().deltaT().value();
    auto timeIndex = particleCloud_.mesh().time().timeIndex();
    // (timeIndex - timeIndexOffset_) * CFDts: 流体时间步截止到 timeIndex 的计算时间
    // (couplingStep_ - 1) * couplingTime(): 耦合时间步截止到 (couplingStep_ - 1) 的计算时间
    return ((timeIndex - timeIndexOffset_) * CFDts - (couplingStep_ - 1) * couplingTime()) / couplingTime();
  }

  int couplingStep() const { return couplingStep_; }

  const scalar& DEMts() const { return DEMts_; }

  int couplingInterval() const { return couplingInterval_; }

  int timeIndexOffset() const { return timeIndexOffset_; }

  //! @brief 计算耦合间隔时长
  scalar couplingTime() const { return couplingInterval_ * DEMts_; }

  scalar TSstart() const {
    return particleCloud_.mesh().time().startTime().value() + (couplingStep_ - 1) * couplingTime();
  }

  scalar TSend() const {
    return particleCloud_.mesh().time().startTime().value() + couplingStep_ * couplingTime();
  }

  int DEMstepsTillT(scalar t) const {
    return (t - (particleCloud_.mesh().time().value() - couplingTime()) + SMALL) / DEMts_;
  }

  void checkTSsize() const {
    auto CFDts = particleCloud_.mesh().time().deltaT().value();
    // 耦合时间步长应该 >= CFD 时间步长
    if (CFDts > couplingTime() + SMALL) {
      Info << "particleCloud_.mesh().time().deltaT().value() = " << CFDts << endl;
      Info << "couplingInterval_ = " << couplingInterval_ << endl;
      Info << "DEMts_ = " << DEMts_ << endl;
      FatalError << "\nError - CFD time-step bigger than coupling time (= DEM time step * coupling interval)!\n" << abort(FatalError);
    }
    // 耦合时间步长应该为CFD 时间步长的整数倍
    if (std::fabs((round(couplingTime() / CFDts) * CFDts) - couplingTime()) > SMALL) {
      Info << "particleCloud_.mesh().time().deltaT().value() = " << CFDts << endl;
      Info << "couplingInterval_ = " << couplingInterval_ << endl;
      Info << "DEMts_ = " << DEMts_ << endl;
      FatalError << "\nWarning - Coupling time (= DEM time step * coupling interval) is not a multiple of  CFD time-step!\n"<< abort(FatalError);
    }
    if (!particleCloud_.allowCFDsubTimestep()) {
      if (CFDts < couplingTime() + SMALL) {
        FatalError << "\nYour models require: CFD time-step = coupling interval (= DEM time step * coupling interval)! \n" << abort(FatalError);
      }
    }
    if (CFDts < couplingTime() - SMALL) {
      Warning << "You are using sub-time-steps (i.e. CFD TS < coupling time)! Check your settings properly." << endl;
    }
  }

  void readDEMtsfromDict(dictionary& propsDict) {
    DEMts_ = readScalar(propsDict.lookup("DEMts"));
    checkTSsize();
  }

  //! @brief 因为耦合时间步长 = 流体时间步长的整数倍，所以判断当前流体时间步是否同时也是耦合时间步
  bool doCoupleNow() const {
    // couple at first sub TS(耦合发生在耦合时间步中的第一个流体时间步中)
    auto CFDts = particleCloud_.mesh().time().deltaT().value();
    auto timeIndex = particleCloud_.mesh().time().timeIndex();
    if ((timeIndex - timeIndexOffset_) * CFDts - SMALL > couplingStep_ * couplingTime()) {
      return true;
    } else {
      return false;
    }
  }

  void setPositions(label n, double* pos) const {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < 3; ++j) {
        particleCloud_.positions_[i][j] = pos[i * 3 + j];
      }
    }
  }

  void setCellIDs(label n, int* ID) const {
    for (int i = 0; i < n; i++) {
      particleCloud_.cellIDs_[i][0] = ID[i];
    }
  }

  void setPositionsCM(label n, double* pos) const {
    particleCloud_.setPositionsCM(n, pos);
  }

  void setCellIDsCM(label n, int* ID) const {
    particleCloud_.setCellIDsCM(n, ID);
  }

  virtual word myType() const = 0;

  virtual void setCG() const { FatalError << "\nsetCG() not executed correctly!\n" << abort(FatalError); }

  void setCGs(LAMMPS_NS::LAMMPS *lmp) const {
    particleCloud_.setCG(lmp->force->cg(int(1e14)));
    Info << "LIGGGHTS hosts " << lmp->atom->ntypes << " atom types." << endl;
#if defined(typeSecificCG)
    Warning << "\nYou are using untested code (type specific coarse graining)!!!\n"
      << "\nThis will most probably give wrong results for most models!!!\n" << endl;
    for (int iTyp = 1; iTyp <= lmp->atom->ntypes; iTyp++) {
      particleCloud_.setCGTypeSpecific(iTyp, lmp->force->cg(iTyp));
    }
#endif
  };

  int getNumberOfParticles() const {
    FatalError << "ask for nr of particles - which is not supported for this dataExchange model" << abort(FatalError);
    return -1;
  }

  int getNumberOfClumps() const {
    FatalError << "ask for nr of clumps - which is not supported for this dataExchange model" << abort(FatalError);
    return -1;
  }

  int getNumberOfTypes() const {
    FatalError << "ask for nr of types - which is not supported for this dataExchange model" << abort(FatalError);
    return -1;
  }

  double* getTypeVol() const {
    FatalError << "ask for type volume - which is not supported for this dataExchange model" << abort(FatalError);
    return NULL;
  }

};  // End of class dataExchangeModel

}  // End of namespace Foam

#endif  // dataExchangeModel_H
