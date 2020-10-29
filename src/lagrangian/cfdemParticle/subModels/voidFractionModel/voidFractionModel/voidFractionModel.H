#ifndef voidFractionModel_H
#define voidFractionModel_H

#include "fvCFD.H"
#include "cfdemCloud.H"

#ifndef __VF_MIX__
#define __VF_MIX__ 1
#endif

namespace Foam {

//! @brief 空隙率模型
class voidFractionModel {

protected:
  const dictionary& dict_;

  cfdemCloud& particleCloud_;

  //! @brief 小颗粒空隙率
  mutable volScalarField voidfractionPrev_;

  mutable volScalarField voidfractionNext_;

#if __VF_MIX__
  //! @brief 大颗粒体积分数
  mutable volScalarField volumefractionPrev_;

  mutable volScalarField volumefractionNext_;
#endif

  //! @brief 颗粒覆盖的网格数量
  //! @note cellsPerParticle_[index][0]
  mutable double** cellsPerParticle_;

  //! @brief 颗粒覆盖的网格数量的上限
  int maxCellsPerParticle_;

  //! @brief fine颗粒覆盖的网格数量的上限
  int maxCellsNumPerFineParticle_;

  //! @brief middle 颗粒覆盖的网格数量的上限
  int maxCellsNumPerMiddleParticle_;

  //! @brief coarse 颗粒覆盖的网格数量的上限
  int maxCellsNumPerCoarseParticle_;

  //! @brief CFD网格中的颗粒所占体积可以通过 weight_参数进行调整
  //        Volumn of particle = weight_ * diameter ^ 3 * pi / 6
  // @defalut 1.0
  mutable scalar weight_;

  //! @brief 颗粒的影响区域通过 porosity_ 使其增大，但是在计算 voidFraction 的时候保持体积不变
  // @defalut 1.0
  mutable scalar porosity_;

  bool requiresSuperquadric_;

public:

  //! @brief Runtime type information
  TypeName("voidFractionModel");

  //! @brief Declare runtime constructor selection table
  declareRunTimeSelectionTable
  (
    autoPtr,
    voidFractionModel,
    dictionary,
    (
      const dictionary& dict,
      cfdemCloud& sm
    ),
    (dict, sm)
  );

  //! @brief Constructors
  voidFractionModel(const dictionary& dict,
                    cfdemCloud& sm);

  //! @brief Destructor
  virtual ~voidFractionModel();

  //! @brief Selector
  static autoPtr<voidFractionModel> New(const dictionary& dict,
                                        cfdemCloud& sm);

  //! @brief 设置空隙率(纯虚函数)
  virtual void setvoidFraction(double** const&,
                               double**&,
                               double**&,
                               double**&,
                               double**&) const = 0;

  //! @brief 计算被 middle 颗粒影响到的网格编号
  // @note cellAffectIDs 用于计算 middle particle 受到的流体阻力
  virtual void getAffectedCellIDs(const std::vector<int>& middleParticleIndexs,
                                  std::vector<std::vector<int> >& affectedCellIDs,
                                  const scalar& searchScale = 4) const {
    FatalError << "voidFractionModel::getAffectedCellIDs(): using base class function, please use derived class function\n"
      << abort(FatalError);
  }

  //! @brief 计算颗粒尺寸与其周围网格平均尺寸的比值
  virtual void getDimensionRatios(std::vector<double>& dimensionRatios,
                                  std::vector<double>& globalDimensionRatios,
                                  std::vector<double>& sumCellsNumbers,
                                  std::vector<double>& sumCellsVolumes,
                                  const std::vector<double*>& centreCellIDs) {
    FatalError << "voidFractionModel::getDimensionRatios(): using base class function, please use derived class function\n"
      << abort(FatalError);
  }

  //! @brief 设置空隙率以及颗粒体积分数前必须调用该初始化函数
  virtual void voidFractionModelInit(double**& particleWeights,
                                     double**& particleVolumes,
                                     double**& particleV,
                                     const std::vector<double>& dimensionRatios) const {}

  //! @brief 设置空隙率以及流体体积分数
  // @note for cfdemSolverMix used
  virtual void setVoidFractionAndVolumeFraction(double** const& mask,
                                                double**& volumefractions,
                                                double**& voidfractions,
                                                double**& particleWeights,
                                                double**& particleVolumes,
                                                double**& particleV,
                                                const std::vector<double>& dimensionRatios) const {}

  void applyDebugSettings(bool debug) const;

  tmp<volScalarField> voidFractionInterp() const;

  //! @brief 重置空隙率
  void resetVoidFractions() const;

#if __VF_MIX__
  //! @brief 重置体积分数
  void resetVolumeFractions() const;
#endif

  void reAllocArrays() const;

  void reAllocArrays(int nP) const;

  void checkWeightNporosity(dictionary& propsDict) const;

  virtual double pointInParticle(int index, vector positionCenter, vector point, double scale) const;

  virtual double pointInParticle(int index, vector positionCenter, vector point) const;

  virtual void setParticleType(label type) const {}

  virtual bool checkParticleType(label type) const { return true; }
  
  virtual double minPeriodicDistance(int index,
                                     vector cellCentrePosition,
                                     vector positionCenter,
                                     boundBox globalBb,
                                     vector& minPeriodicPos,
                                     vector dirCheckRange = vector(1, 1, 1)) const;

  // Access Functions

  inline volScalarField& voidFractionPrev() const { return voidfractionPrev_; }

  inline volScalarField& voidFractionNext() const { return voidfractionNext_; }

  inline scalar weight() const { return weight_; }

  inline scalar porosity() const { return porosity_; }

  inline double ** const& cellsPerParticle() const { return cellsPerParticle_; }

  inline int maxCellsPerParticle() const { return maxCellsPerParticle_; }

  inline int maxCellsNumPerFineParticle() const { return maxCellsNumPerFineParticle_; }

  inline int maxCellsNumPerMiddleParticle() const { return maxCellsNumPerMiddleParticle_; }

  inline int maxCellsNumPerCoarseParticle() const { return maxCellsNumPerCoarseParticle_; }

  inline bool requiresSuperquadric() const { return requiresSuperquadric_; }

#if __VF_MIX__
  inline volScalarField& volumeFractionPrev() const { return volumefractionPrev_; }

  inline volScalarField& volumeFractionNext() const { return volumefractionNext_; }
#endif

};  // End of class voidFractionModel

}  // End of namespace Foam

#endif  // voidFractionModel_H
