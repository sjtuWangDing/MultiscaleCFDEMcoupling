#ifndef mixVoidFraction_H
#define mixVoidFraction_H

#include "voidFractionModel.H"

namespace Foam {

// @brief 混合颗粒场空隙率模型
class mixVoidFraction: public voidFractionModel {

private:
  // @brief 当颗粒的体积在CFD单元的尺寸范围内时，应使用分割(divided)空隙率模型。粒子的半径为 R ，颗粒的体积被划分为 29 个等体积的不重叠区域
  static const int numberOfMarkerPoints = 29;

  dictionary propsDict_;

  // @brief 是否打印详细信息
  bool verbose_;

  // @brief 是否允许在处理器边界修正，这需要使用 engineIB
  Switch procBoundaryCorrection_;

  // @brief 空隙率的最小值
  const scalar alphaMin_;

  // @brief true - 空隙率取 alphaMin_
  //        false - 空隙率没有被 alphaMin_修正
  mutable bool alphaLimited_;

  // @brief 由于空隙率限制而导致颗粒体积的累积损失
  mutable scalar tooMuch_;

  bool cfdemUseOnly_;

  // @brief 标志点坐标
  // @note 这里 offsets = (标志点中心到颗粒中心的距离 / 半径)，不是实际坐标
  vector offsets[numberOfMarkerPoints];

  // @brief 在调用 setvoidFraction (for mix solver) 函数前是否调用过 mixVoidFractionInit 函数初始化
  mutable bool hasInit_;

  // // @brief 体积扩大因数
  // // @note scaleUpVol_ >= 1
  // const scalar scaleUpVol_;

public:
  // @brief Runtime type information
  TypeName("mix");

  // @brief 构造函数
  mixVoidFraction(const dictionary& dict, cfdemCloud& sm);

  // @brief 析构函数
  ~mixVoidFraction() {}

  // @brief 设置空隙率
  // @note for other solvers used
  virtual void setvoidFraction(double** const& mask,
                               double**& voidfractions,
                               double**& particleWeights,
                               double**& particleVolumes,
                               double**& particleV) const;

  // @brief 计算被 middle 颗粒影响到的网格编号
  // @note cellAffectIDs 用于计算 middle particle 受到的流体阻力, 搜索范围默认是 4 倍颗粒半径
  virtual void getAffectedCellIDs(const std::vector<int>& middleParticleIndexs,
                                  std::vector<std::vector<int> >& affectedCellIDs,
                                  const scalar& searchScale = 4) const;

  // @brief 计算颗粒尺寸与其周围网格平均尺寸的比值, 并将颗粒索引按照颗粒尺寸归类
  virtual void getDimensionRatios(std::vector<double>& dimensionRatios);

  // @brief 设置空隙率以及颗粒体积分数前必须调用该初始化函数
  virtual void voidFractionModelInit(double**& voidfractions,
                                     double**& particleWeights,
                                     double**& particleVolumes,
                                     double**& particleV,
                                     const std::vector<double>& dimensionRatios) const;

  // @brief 设置空隙率以及流体体积分数
  // @note for cfdemSolverMix used
  void setVoidFractionAndVolumeFraction(double** const& mask,
                                        double**& volumefractions,
                                        double**& voidfractions,
                                        double**& particleWeights,
                                        double**& particleVolumes,
                                        double**& particleV,
                                        const std::vector<double>& dimensionRatios) const;

  // @brief 设置索引为 index 的单个颗粒的空隙率
  void setvoidFractionForSingleParticle(const int index,
                                        double**& particleWeights,
                                        double**& particleVolumes,
                                        double**& particleV) const;

  // @brief 设置索引为 index 的单个颗粒的体积分数场
  virtual void setVolumeFractionForSingleParticle(const int index) const;

  /*
   * @brief (模糊搜索!!)构建颗粒覆盖的所有网格的哈希集合
   * @note 设置为递归函数，通过哈希器将网格编号转换为哈希值，并存入 set 中以便于搜索
   * @param index     <[in] 颗粒索引
   * @param position  <[in] 颗粒中心位置
   * @param cellID    <[in] 递归循环中要检索网格编号
   * @param hashSett  <[in, out] 需要构建的哈希集
   * @param scale     <[in] 颗粒半径扩大系数
   */
  void buildLabelHashSetForDimensionRatio(const int index,
                                          const vector& position,
                                          const label& cellID,
                                          labelHashSet& hashSett,
                                          const double& scale) const;

  /*
  * @brief 构建颗粒覆盖的所有网格的哈希集合
  * @note 设置为递归函数,  通过哈希器将网格编号转换为哈希值,  并存入 set 中以便于搜索
  * @param index          <[in] 颗粒索引
  * @param position       <[in] 颗粒中心位置
  * @param cellID         <[in] 递归循环中要检索网格编号
  * @param hashSett       <[in, out] 需要构建的哈希集
  * @param initialInsert
  */
  void buildLabelHashSetForVolumeFractions(int index,
                                           const vector position,
                                           const label cellID,
                                           labelHashSet& hashSett,
                                           bool initialInsert) const;

private:

  // @brief 获取 Corona Point 坐标
  vector getCoronaPointPosition(vector const& position,
                                vector const& cellCentrePosition,
                                scalar const corona) const;

  // @brief 计算距离系数
  double segmentParticleIntersection(int index,
                                     vector positionCenter,
                                     vector pointInside,
                                     vector pointOutside) const;

};  // End of class mixVoidFraction

}  // End of namespace Foam

#endif  // mixVoidFraction_H
