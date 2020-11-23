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

#include "cloud/cfdemCloud.H"
#include "base/runTimeSelectionTables.H"

namespace Foam {

class dataExchangeModel {

public:
  /* --------------------------------- base class and functions -------------------------- */

  //! \brief Runtime type information
  cfdemTypeName("dataExchangeModel")

  //! \brief Declare runtime constructor selection
  cfdemDeclareRunTimeSelection(autoPtr, dataExchangeModel, (cfdemCloud& cloud), (cloud))

  //! \brief Selector
  static autoPtr<dataExchangeModel> New(cfdemCloud& cloud, const dictionary& dict);

  //! \brief Constructor
  dataExchangeModel(cfdemCloud& cloud);

  //! \brief Destructor
  virtual ~dataExchangeModel();

  /* --------------------------------- virtual functions --------------------------------- */

  //! \return 当前耦合时间步中颗粒的数量
  virtual int couple() = 0;

  //! \brief Allocate and destroy for 2-D double array
  virtual void destroy(double**& array) {};
  virtual void destroy(double**& array, int) {};
  virtual void allocateArray(double**& array, double initVal, int width, int length) {};
  virtual void allocateArray(double**& array, double initVal, int width, const char* length = "nparticles") {};

  //! \brief Allocate and destroy for 2-D int array
  virtual void destroy(int**& array) {};
  virtual void destroy(int**& array, int) {};
  virtual void allocateArray(int**& array, int initVal, int width, int length) {};
  virtual void allocateArray(int**& array, int initVal, int width, const char* length  = "nparticles") {};

  //! \brief Allocate and destroy for 1-D double array
  virtual void destroy(double*& array) {};
  virtual void allocateArray(double*& array, double initVal, int length) {};

  //! \brief Allocate and destroy for 1-D int array
  virtual void destroy(int*& array) {};
  virtual void allocateArray(int*& array, int initVal, int length) {};

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

  //! \brief 释放离散内存
  virtual void destroyDiscreteMemory(double** const& array, int length) {
    FatalError << "dataExchangeModel::destroyDiscreteMemory(): "
      << "using base class function, please use derived class function\n"
      << abort(FatalError);
  }
  //! \brief 释放离散内存
  virtual void destroyDiscreteMemory(int** const& array, int length) {
    FatalError << "dataExchangeModel::destroyDiscreteMemory(): "
      << "using base class function, please use derived class function\n"
      << abort(FatalError);
  }

  /* --------------------------------- member functions ---------------------------------- */

  template<typename DType>
  void getData(const std::string& dataName,
               const std::string& dataType,
               DType** const& field) {
    getData(dataName, dataType, field, couplingStep_ - 1);
  }

  /*!
   * \brief 因为耦合时间步长 = 流体时间步长的整数倍，所以 timeStepFraction() 用于计算每个流体时间步在耦合时间步中的所占比例，
   *        如果 couplingTime() == 3 * CFDts，那么每一个耦合时间步由 3 个流体时间步构成，
   *        那么这三个流体时间步的 timeStepFraction() 分别返回 0, 0.333333, 0.666666
   */
  double timeStepFraction() const {
    scalar CFDts = cloud_.mesh().time().deltaT().value();
    label timeIndex = cloud_.mesh().time().timeIndex();
    // (timeIndex - timeIndexOffset_) * CFDts: 耦合计算后，CFD 的计算时间
    // (couplingStep_ - 1) * couplingTime(): 耦合时间步截止到 (couplingStep_ - 1) 的计算时间
    return ((timeIndex - timeIndexOffset_) * CFDts - (couplingStep_ - 1) * couplingTime()) / couplingTime();
  }

  //! \brief 检查时间步长是否满足耦合要求
  void checkTimeStepSize() const {
    scalar CFDts = cloud_.mesh().time().deltaT().value(); // 获取 CFD 时间步长

    // CFD 时间步长 > 耦合时间步长
    if (CFDts > couplingTime() + SMALL) {
      Info << "cloud_.mesh().time().deltaT().value() = " << CFDts << endl;
      Info << "couplingInterval = " << cloud_.cProps().couplingInterval() << endl;
      Info << "DEMts = " << DEMts_ << endl;
      FatalError << "\nError - CFD time-step bigger than coupling time (= DEM time step * coupling interval)!\n" << abort(FatalError);
    }

    // 耦合时间步长应该为 CFD 时间步长的整数倍
    if (std::fabs((round(couplingTime() / CFDts) * CFDts) - couplingTime()) > SMALL) {
      Info << "cloud_.mesh().time().deltaT().value() = " << CFDts << endl;
      Info << "couplingInterval = " << cloud_.cProps().couplingInterval() << endl;
      Info << "DEMts = " << DEMts_ << endl;
      FatalError << "\nWarning - Coupling time (= DEM time step * coupling interval) is not a multiple of  CFD time-step!\n"<< abort(FatalError);
    }

    // CFD 时间步长 < 耦合时间步长
    if (CFDts < couplingTime() - SMALL) {
      // 如果不允许使用 CFD sub time step 则报错
      if (!cloud_.cProps().allowCFDsubTimeStep()) {
        FatalError << "\nYour models require: CFD time-step = coupling interval (= DEM time step * coupling interval)! \n" << abort(FatalError);
      }
      Warning << "You are using sub-time-steps (i.e. CFD TS < coupling time)! Check your settings properly." << endl;
    }
  }

  //! \brief 因为耦合时间步长 = 流体时间步长的整数倍，所以判断当前流体时间步是否同时也是耦合时间步
  bool doCoupleNow() const {
    // couple at first sub TS(耦合发生在耦合时间步中的第一个流体时间步中)
    scalar CFDts = cloud_.mesh().time().deltaT().value();
    label timeIndex = cloud_.mesh().time().timeIndex();
    // (timeIndex - timeIndexOffset_) * CFDts: 耦合计算后，CFD 的计算时间
    // 如果当前流体时间步不是耦合时间步，则 CFD 计算时间一定小于 couplingStep_ * couplingTime() 计算出的时间
    if ((timeIndex - timeIndexOffset_) * CFDts - SMALL > couplingStep_ * couplingTime()) {
      return true;
    } else {
      return false;
    }
  }

  /* --------------------------------- public interfaces --------------------------------- */

  inline int couplingStep() const { return couplingStep_; }

  inline scalar DEMts() const { return DEMts_; }

  //! \brief 计算耦合间隔的秒数
  inline scalar couplingTime() const { return cloud_.cProps().couplingInterval() * DEMts_; }

  //! \brief 当前耦合时间步的起始时间（全局时间）
  inline scalar TSstart() const {
    return cloud_.mesh().time().startTime().value() + (couplingStep_ - 1) * couplingTime();
  }

  //! \brief 当前耦合时间步的结束时间（全局时间）
  inline scalar TSend() const {
    return cloud_.mesh().time().startTime().value() + couplingStep_ * couplingTime();
  }

  //! \brief 计算指定时间的 DEM 步
  inline int DEMstepsTillT(scalar t) const {
    return (t - (cloud_.mesh().time().value() - couplingTime()) + SMALL) / DEMts_;
  }

protected:

  cfdemCloud& cloud_;

  //! \brief 初始化 dataExchangeModel 的时候，即耦合开始的时候，记录当前的流体时间步，通常从 0 开始
  const int timeIndexOffset_;

  //! \brief 耦合时间步的计数，每执行一次 coupling 记录一次
  int couplingStep_;

  //! \brief DEM 的时间步长
  scalar DEMts_;
};

} // namespace Foam

#endif // __DATA_EXCHANGE_MODEL_H__