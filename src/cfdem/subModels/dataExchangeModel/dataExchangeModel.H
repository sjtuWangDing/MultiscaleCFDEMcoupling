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

  //! \brief 耦合时间步计数器
  class CouplingStepCounter {
  public:
    CouplingStepCounter(const dataExchangeModel& model)
      : model_(model) {
      Info << "Current coupling time step = " << model_.couplingStep() << endl;
    }
    ~CouplingStepCounter() {
      const_cast<dataExchangeModel&>(model_).couplingStep_ += 1;
    }
  private:
    const dataExchangeModel& model_;
  };

  //! \brief 耦合时间步计数器
  template<bool UseTwoWayMPI, typename Tensor>
  void allocate(Tensor& tensor, const typename Tensor::TData& initVal);

  template<int nDim, typename DType, typename Device CFDEM_DEFAULT_DEVICE>
  inline void reAlloc(base::Tensor<nDim, DType, Device>& tensor,
                      DType**& pData,
                      int length, int width) const {
    CHECK(2 == nDim || 1 == nDim) << "tensor realloc: error tensor dimension";
    static base::simple_alloc<DType, base::malloc_alloc> alloc1;
    static base::simple_alloc<DType*, base::malloc_alloc> alloc2;
    if (nullptr == pData && nullptr == tensor.dptr_) {
      // 直接调用 LIGGGHTS 分配内存到 pData
      allocateArray(pData, length, width);
    } else if (nullptr != pData && nullptr != tensor.dptr_) {
      CHECK_EQ(pData[0], tensor.dptr_) << "tensor realloc: tensor not match with pData";
      // 释放 pData[0] 指向的内存
      alloc1.deallocate(pData[0], tensor.mSize());
      // 释放 pData 指向的内存
      alloc2.deallocate(pData, tensor.shape_[0]);
    } else {
      CHECK(false) << "tensor realloc: bad realloc";
    }
    // 调用 LIGGGHTS 分配内存到 pData
    allocateArray(pData, length, width);
    // 将 pData[0] 转移至 tensor
    tensor.dptr_ = pData[0];
    tensor.stride_ = 1 == nDim ? length : width;
    tensor.shape_ = 1 == nDim ? base::makeShape<nDim>({ length }) : base::makeShape<nDim>({ length, width });
  }

  /* --------------------------------- virtual functions --------------------------------- */

  //! \return 当前耦合时间步中颗粒的数量
  virtual int couple() = 0;
  /*!
   * \brief 从 LAMMPS 中获取数据（数据类型为 double）
   *        Eg: getData("x", "vector-atom", array);
   * \param dataName 数据名称
   * \param dataType 数据类型
   * \param array 内存地址（LIGGGHTS code 要求 array 为二级指针）
   */
  virtual void getData(const std::string& dataName,
                       const std::string& dataType,
                       double** array) const = 0;
  /*!
   * \brief 从 LAMMPS 中获取数据（数据类型为 int）
   * \param dataName 数据名称
   * \param dataType 数据类型
   * \param array 内存地址（LIGGGHTS code 要求 array 为二级指针）
   */
  virtual void getData(const std::string& dataName,
                       const std::string& dataType,
                       int** array) const = 0;
  /*!
   * \brief 传递数据到 LAMMPS 中（数据类型为 double）
   *        Eg: giveData("dragforce", "vector-atom", array);
   * \param dataName 数据名称
   * \param dataType 数据类型
   * \param array 内存地址（LIGGGHTS code 要求 array 为二级指针）
   */
  virtual void giveData(const std::string& dataName,
                        const std::string& dataType,
                        double** field) const = 0;
  /*!
   * \brief 传递数据到 LAMMPS 中（数据类型为 int）
   * \param dataName 数据名称
   * \param dataType 数据类型
   * \param array 内存地址（LIGGGHTS code 要求 array 为二级指针）
   */
  virtual void giveData(const std::string& dataName,
                        const std::string& dataType,
                        int** field) const = 0;
  //! \brief Allocate for 2-D double array
  virtual void allocateArray(double**& array, int length, int width, double initVal = 0.0) const {}
  //! \brief Allocate 2-D int array
  virtual void allocateArray(int**& array, int length, int width, int initVal = 0) const {}

  // //! \brief Allocate and destroy for 2-D double array
  // virtual void destroy(double**& array) {};
  // virtual void destroy(double**& array, int) {};
  // virtual void allocateArray(double**& array, double initVal, int width, int length) {};
  // virtual void allocateArray(double**& array, double initVal, int width, const char* length = "nparticles") {};

  // //! \brief Allocate and destroy for 2-D int array
  // virtual void destroy(int**& array) {};
  // virtual void destroy(int**& array, int) {};
  // virtual void allocateArray(int**& array, int initVal, int width, int length) {};
  // virtual void allocateArray(int**& array, int initVal, int width, const char* length  = "nparticles") {};

  // //! \brief Allocate and destroy for 1-D double array
  // virtual void destroy(double*& array) {};
  // virtual void allocateArray(double*& array, double initVal, int length) {};

  // //! \brief Allocate and destroy for 1-D int array
  // virtual void destroy(int*& array) {};
  // virtual void allocateArray(int*& array, int initVal, int length) {};

  // virtual void getData(const std::string& dataName,
  //                      const std::string& dataType,
  //                      double** const& field,
  //                      label step) = 0;
  // virtual void getData(const std::string& dataName,
  //                      const std::string& dataType,
  //                      int** const& field,
  //                      label step) = 0;
  // virtual void giveData(const std::string& dataName,
  //                       const std::string& dataType,
  //                       double** const& field,
  //                       const char* fieldType = "double") = 0;

  // //! \brief 释放离散内存
  // virtual void destroyDiscreteMemory(double** const& array, int length) {
  //   FatalError << "dataExchangeModel::destroyDiscreteMemory(): "
  //     << "using base class function, please use derived class function\n"
  //     << abort(FatalError);
  // }
  // //! \brief 释放离散内存
  // virtual void destroyDiscreteMemory(int** const& array, int length) {
  //   FatalError << "dataExchangeModel::destroyDiscreteMemory(): "
  //     << "using base class function, please use derived class function\n"
  //     << abort(FatalError);
  // }

  /* --------------------------------- member functions ---------------------------------- */

  /*!
   * \brief 因为耦合时间步长 = 流体时间步长的整数倍，所以 timeStepFraction() 用于计算每个流体时间步在耦合时间步中的所占比例，
   *        如果 couplingTime() == 3 * CFDts，那么每一个耦合时间步由 3 个流体时间步构成，
   *        那么这三个流体时间步的 timeStepFraction() 分别返回 0, 0.333333, 0.666666
   */
  inline double timeStepFraction() const {
    // 如果当前流体时间步是耦合时间步中第一个流体时间步
    if (checkValidCouplingStep()) {
      return 0;
    }
    scalar CFDts = cloud_.mesh().time().deltaT().value();
    label timeIndex = cloud_.mesh().time().timeIndex();
    // (timeIndex - timeIndexOffset_ - 1) * CFDts: 当前流体时间步起始时间，Eg: 0.00, 0.01, 0.02, 0.03
    // (couplingStep_ - 1) * couplingTime(): 耦合时间步截止到 (couplingStep_ - 1) 的计算时间
    return ((timeIndex - timeIndexOffset_ - 1) * CFDts - (couplingStep_ - 1) * couplingTime()) / couplingTime();
  }

  /*!
   * \brief 检查时间步长是否满足耦合要求（在初始化模型的时候调用）
   * \note (1) 耦合时间步长应该大于 or 等于 CFD 时间步长
   *       (2) 耦合时间步长应该为 CFD 时间步长的整数倍
   *       (3) 如果耦合时间步长 != CFD 时间步长，则要求 allowUseSubCFDTimeStep() 为 true
   */
  inline void checkTimeStepSize() const {
    scalar CFDts = cloud_.mesh().time().deltaT().value(); // 获取 CFD 时间步长
    Info << "CFD time step size: " << CFDts << endl;
    // 耦合时间步长 >= CFD 时间步长
    if (CFDts > couplingTime() + Foam::SMALL) {
      Info << "CFD time step size = " << CFDts << endl;
      Info << "DEM time step size = " << DEMts_ << endl;
      Info << "Coupling interval = " << cloud_.cProps().couplingInterval() << endl;
      FatalError << "CFD time-step bigger than coupling time (= DEM time step * coupling interval)!\n"
        << abort(FatalError);
    }
    // 耦合时间步长应该为 CFD 时间步长的整数倍
    if (std::fabs((std::round(couplingTime() / CFDts) * CFDts) - couplingTime()) > Foam::SMALL) {
      Info << "CFD time step size = " << CFDts << endl;
      Info << "DEM time step size = " << DEMts_ << endl;
      Info << "Coupling interval = " << cloud_.cProps().couplingInterval() << endl;
      FatalError << "Coupling time (= DEM time step * coupling interval) is not a multiple of  CFD time-step!\n"<< abort(FatalError);
    }
    // 如果耦合时间步长 != CFD 时间步长
    if (couplingTime() > CFDts + Foam::SMALL) {
      // 如果不允许使用 CFD sub time step 则报错
      if (!cloud_.cProps().allowUseSubCFDTimeStep()) {
        FatalError << "Your models require: CFD time-step = coupling interval (= DEM time step * coupling interval)! \n" << abort(FatalError);
      }
      Warning << "You are using sub-time-steps (i.e. CFD time-step < coupling time)" << endl;
    }
  }

  /*!
   * \brief 因为耦合时间步长 = 流体时间步长的整数倍，而耦合发生在耦合时间步中的第一个流体时间步中，
   *        所以判断当前流体时间步是否同时也是耦合时间步
   */
  inline bool checkValidCouplingStep() const {
    scalar CFDts = cloud_.mesh().time().deltaT().value();
    // OpenFOAM 中 timeIndex 从 1 开始
    label timeIndex = cloud_.mesh().time().timeIndex();
    // (timeIndex - timeIndexOffset_ - 1) * CFDts: 当前流体时间步起始时间，Eg: 0.00, 0.01, 0.02, 0.03
    // couplingStep_ * couplingTime(): 当前耦合时间步的起始时间，Eg: 0.00, 0.03, 0.06
    if ((timeIndex - timeIndexOffset_ - 1) * CFDts + Foam::SMALL > couplingStep_ * couplingTime()) {
      return true;
    }
    return false;
  }

  /* --------------------------------- public interfaces --------------------------------- */

  inline int couplingStep() const { return couplingStep_; }

  inline double DEMts() const { return DEMts_; }

  //! \brief 计算耦合间隔的秒数
  inline double couplingTime() const { return cloud_.cProps().couplingInterval() * DEMts_; }

  //! \brief 当前耦合时间步的起始时间（全局时间）
  inline double TSstart() const {
    return cloud_.mesh().time().startTime().value() + (couplingStep_ - 1) * couplingTime();
  }

  //! \brief 当前耦合时间步的结束时间（全局时间）
  inline double TSend() const {
    return cloud_.mesh().time().startTime().value() + couplingStep_ * couplingTime();
  }

  //! \brief 计算指定时间的 DEM 步
  inline int DEMstepsTillT(double t) const {
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