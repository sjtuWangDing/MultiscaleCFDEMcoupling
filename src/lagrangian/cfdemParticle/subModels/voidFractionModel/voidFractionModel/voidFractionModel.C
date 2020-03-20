#include "error.H"
#include "voidFractionModel.H"
#include "dataExchangeModel.H"

namespace Foam {
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(voidFractionModel, 0);

defineRunTimeSelectionTable(voidFractionModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

voidFractionModel::voidFractionModel(const dictionary& dict,
                                     cfdemCloud& sm):
  dict_(dict),
  particleCloud_(sm),
  voidfractionPrev_(
    IOobject(
      "voidfractionPrev",
      sm.mesh().time().timeName(),
      sm.mesh(),
      IOobject::READ_IF_PRESENT,  // MUST_READ,
      IOobject::AUTO_WRITE
    ),
    // sm.mesh().lookupObject<volScalarField>("voidfraction")
    sm.mesh(),
    dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1)
  ),
  voidfractionNext_(
    IOobject(
      "voidfractionNext",
      sm.mesh().time().timeName(),
      sm.mesh(),
      IOobject::READ_IF_PRESENT,  // MUST_READ,
      IOobject::AUTO_WRITE
    ),
    // sm.mesh().lookupObject<volScalarField>("voidfraction")
    sm.mesh(),
    dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1)
  ),
#if __VF_MIX__
  volumefractionPrev_(
    IOobject(
      "volumefractionPrev",
      sm.mesh().time().timeName(),
      sm.mesh(),
      IOobject::READ_IF_PRESENT, // MUST_READ,
      IOobject::AUTO_WRITE
    ),
    // sm.mesh().lookupObject<volScalarField>("volumefraction")
    sm.mesh(),
    dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1)
  ),
  volumefractionNext_(
    IOobject(
      "volumefractionNext",
      sm.mesh().time().timeName(),
      sm.mesh(),
      IOobject::READ_IF_PRESENT, // MUST_READ,
      IOobject::AUTO_WRITE
    ),
    // sm.mesh().lookupObject<volScalarField> ("volumefraction")
    sm.mesh(),
    dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1)
  ),
#endif  // __VF_MIX__
  cellsPerParticle_(NULL),
  maxCellsPerParticle_(1),
  maxCellsNumPerFineParticle_(1),
  maxCellsNumPerMiddleParticle_(1),
  maxCellsNumPerCoarseParticle_(1),
  weight_(1.),
  porosity_(1.),
  requiresSuperquadric_(false) {
  // 分配 cellsPerParticle_ 内存
    particleCloud_.dataExchangeM().allocateArray(
    cellsPerParticle_,   // <[in, out] 指定分配的二维数组指针
    1,                   // <[in] 指定数组的初值
    1                    // <[in] width
  );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

voidFractionModel::~voidFractionModel() {
  particleCloud_.dataExchangeM().destroy(cellsPerParticle_, 1);
}

// * * * * * * * * * * * * * *  Member Functions  * * * * * * * * * * * ** * //

void voidFractionModel::applyDebugSettings(bool debug) const {
  if (!debug) {
    voidfractionPrev_.writeOpt() = IOobject::NO_WRITE;
    voidfractionNext_.writeOpt() = IOobject::NO_WRITE;
#if __VF_MIX__
    volumefractionPrev_.writeOpt() = IOobject::NO_WRITE;
    volumefractionNext_.writeOpt() = IOobject::NO_WRITE;
#endif
  }
}

// @brief 计算空隙率时间插值
tmp<volScalarField> voidFractionModel::voidFractionInterp() const {
  scalar tsf = particleCloud_.dataExchangeM().timeStepFraction();
  Info << "using voidfraction blend, tsf=" << tsf << endl;
  return tmp<volScalarField>(
    new volScalarField("alpha_voidFractionModel", (1 - tsf) * voidfractionPrev_ + tsf * voidfractionNext_)
  );
}

// @brief 重置空隙率
void voidFractionModel::resetVoidFractions() const {
  voidfractionPrev_ == voidfractionNext_;
  voidfractionNext_ == dimensionedScalar("one", voidfractionNext_.dimensions(), 1.);
}

#if __VF_MIX__
// @brief 重置体积分数
void voidFractionModel::resetVolumeFractions() const {
  volumefractionPrev_ == volumefractionNext_;
  volumefractionNext_ == dimensionedScalar("one", volumefractionNext_.dimensions(), 1.);
}
#endif

void voidFractionModel::reAllocArrays() const {
  if (particleCloud_.numberOfParticlesChanged()) {
    // get arrays of new length
    particleCloud_.dataExchangeM().allocateArray(cellsPerParticle_, 1, 1);
  }
}

void voidFractionModel::reAllocArrays(int nP) const {
  if (particleCloud_.numberOfParticlesChanged()) {
    // get arrays of new length
    particleCloud_.dataExchangeM().allocateArray(cellsPerParticle_, 1, 1, nP);
  }
}

// @brief 判断某一点是否在颗粒中
// @return < 0 -- 在颗粒中
//         > 0 -- 在颗粒外
double voidFractionModel::pointInParticle(int index,
                                          vector positionCenter,
                                          vector point,
                                          double scale) const {
  scalar radius =  particleCloud_.radius(index);
  if (radius > SMALL) {
    scalar pointDistSq = magSqr(point - positionCenter);
    return pointDistSq / (scale * scale * radius * radius) - 1.0;
  } else {
    return 0.;
  }
}

double voidFractionModel::pointInParticle(int index,
                                          vector positionCenter,
                                          vector point) const {
  return pointInParticle(index, positionCenter, point, 1.0);
}

void voidFractionModel::checkWeightNporosity(dictionary& propsDict) const {
  if (propsDict.found("weight")) { 
    weight_ = readScalar(propsDict.lookup("weight"));
  }
  if (propsDict.found("porosity")) {
    porosity_ = readScalar(propsDict.lookup("porosity"));
  }
}

// @brief 如果边界是周期性边界，获取流场中某一点到某个颗粒或者其镜像的最小距离
//        并返回颗粒的位置坐标(可能是镜像颗粒的坐标，也可能是原颗粒的坐标)
double voidFractionModel::minPeriodicDistance(
  int index,                  // <[in] 颗粒的编号
  vector cellCentrePosition,  // <[in] 网格中心的位置矢量
  vector positionCenter,      // <[in] 颗粒中心的位置矢量
  boundBox globalBd,          // <[in] Boundary Box
  vector& minPeriodicPos,     // <[in, out] 与网格中心点距离最小的(镜像)颗粒坐标
  vector dirCheckRange        // <[in] 检查的方向，默认是 (1, 1, 1)
) const{
  /*
    * 如果 dirCheckRange = (1, 0, 0)
    *                                  globalBb
    *                         --------------------
    *                         |                  |
    *                         |                  |
    *                         |                  |
    * (颗粒镜像位置1) *         |  *(颗粒实际位置)   |  * (颗粒镜像位置2)
    *                         |                  |
    *                         |                + |
    *                         |           (mesh) |
    *                         -------------------- ------> x 方向
    * 如图所示，如果 x 方向上为周期边界条件，则对于图中位置处的 mesh，
    * 距离其最近的颗粒位置应该是 镜像位置2，所以参数 minPeriodicPos 就
    * 被设置为镜像位置2 
    */
  double f = 999e32;
  vector positonCenterPeriodic;

  for(int xDir = - static_cast<int>(dirCheckRange[0]);
      xDir <= static_cast<int>(dirCheckRange[0]); ++ xDir){

    positonCenterPeriodic[0] = positionCenter[0] + static_cast<double>(xDir) * (globalBd.max()[0] - globalBd.min()[0]);

    for(int yDir = - static_cast<int>(dirCheckRange[1]);
        yDir <= static_cast<int>(dirCheckRange[1]); ++ yDir){

      positonCenterPeriodic[1] = positionCenter[1] + static_cast<double>(yDir) * (globalBd.max()[1] - globalBd.min()[1]);

      for(int zDir = - static_cast<int>(dirCheckRange[2]);
          zDir <= static_cast<int>(dirCheckRange[2]); ++ zDir){

        positonCenterPeriodic[2] = positionCenter[2] + static_cast<double>(zDir) * (globalBd.max()[2] - globalBd.min()[2]);

        // 此时 positonCenterPeriodic 就是颗粒的某个镜像坐标
        // 当 xDir = yDir = zDir = 0 时，positonCenterPeriodic 才是颗粒本身的坐标
        if (pointInParticle(index, positonCenterPeriodic, cellCentrePosition) < f) {
          f = pointInParticle(index, positonCenterPeriodic, cellCentrePosition);
          minPeriodicPos = positonCenterPeriodic;
        }
      }
    }
  }
  return f;
}

}  // End of namespace Foam
