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
  dataExchangeModel
\*---------------------------------------------------------------------------*/

#ifndef __DATA_EXCHANGE_MODEL_H__
#define __DATA_EXCHANGE_MODEL_H__

#include "fvCFD.H"
#include "cloud/cfdemCloud.H"

#include "tools/runTimeSelectionTables.H"

namespace Foam {

class dataExchangeModel {

public:

  //! @brief Runtime type information
  cfdemTypeName("dataExchangeModel");
  //! @brief Declare runtime constructor selection
  cfdemDeclareRunTimeSelection(autoPtr, dataExchangeModel, (cfdemCloud& cloud), (cloud));
  //! @brief Selector
  static autoPtr<dataExchangeModel> New(cfdemCloud& cloud);

  //! @brief Constructor
  dataExchangeModel(cfdemCloud& cloud);
  //! @brief Destructor
  virtual ~dataExchangeModel();

public:

  virtual bool couple(int);

  virtual void getData(const std::string& dataName,
                       const std::string& dataType,
                       double** const& field,
                       label step) = 0;
  virtual void getData(const std::string& dataName,
                       const std::string& dataType,
                       int** const& field,
                       label step) = 0;
  virtual void giveData(const std::string& dataName,
                        const std::string& dataType,
                        double** const& field,
                        const char* fieldType = "double") = 0;

  //! @brief Allocate and destroy for 2-D double array
  virtual void destroy(double**& array);
  virtual void destroy(double**& array, int);
  virtual void allocateArray(double**& array, double initVal, int width, int length);
  virtual void allocateArray(double**& array, double initVal, int width, const char* length = "nparticles");

  //! @brief Allocate and destroy for 2-D int array
  virtual void destroy(int**& array);
  virtual void destroy(int**& array, int);
  virtual void allocateArray(int**& array, int initVal, int width, int length);
  virtual void allocateArray(int**& array, int initVal, int width, const char* length  = "nparticles");

  //! @brief Allocate and destroy for 1-D double array
  virtual void destroy(double*& array);
  virtual void allocateArray(double*& array, double initVal, int length);

  //! @brief Allocate and destroy for 1-D int array
  virtual void destroy(int*& array);
  virtual void allocateArray(int*& array, int initVal, int length);

  //! @brief 释放离散内存
  virtual void destroyDiscreteMemory(double** const& array, int length) {
    FatalError << "dataExchangeModel::destroyDiscreteMemory(): "
      << "using base class function, please use derived class function\n"
      << abort(FatalError);
  }
  //! @brief 释放离散内存
  virtual void destroyDiscreteMemory(int** const& array, int length) {
    FatalError << "dataExchangeModel::destroyDiscreteMemory(): "
      << "using base class function, please use derived class function\n"
      << abort(FatalError);
  }

public:

  template<typename DType>
  void getData(const std::string& dataName,
               const std::string& dataType,
               DType** const& field) {
    getData(dataName, dataType, field, couplingStep_ - 1);
  }
  /*!
   * @brief 因为耦合时间步长 = 流体时间步长的整数倍，所以 timeStepFraction() 用于计算每个流体时间步在耦合时间步中的所占比例，
   *        如果 couplingTime() == 3 * CFDts，那么每一个耦合时间步由 3 个流体时间步构成，
   *        那么这三个流体时间步的 timeStepFraction() 分别返回 0, 0.333333, 0.666666
   */
  double timeStepFraction() const {
    auto CFDts = cloud_.mesh().time().deltaT().value();
    auto timeIndex = cloud_.mesh().time().timeIndex();
    // (timeIndex - timeIndexOffset_) * CFDts: 流体时间步截止到 timeIndex 的计算时间
    // (couplingStep_ - 1) * couplingTime(): 耦合时间步截止到 (couplingStep_ - 1) 的计算时间
    return ((timeIndex - timeIndexOffset_) * CFDts - (couplingStep_ - 1) * couplingTime()) / couplingTime();
  }

public:

  inline int couplingStep() const { return couplingStep_; }
  inline int couplingInterval() const { return couplingInterval_; }
  inline double DEMts() const { return DEMts_; }

protected:

  cfdemCloud& cloud_;

  const dictionary& dict_;

  //! @brief 耦合时间步的计数
  int couplingStep_;

  //! @brief 耦合间隔，单位：DEM 时间步
  int couplingInterval_;

  //! @brief DEM 的时间步长
  double DEMts_;

  // int maxNumberOfParticles_;



  // const int timeIndexOffset_;

};

} // namespace Foam

#endif // __DATA_EXCHANGE_MODEL_H__