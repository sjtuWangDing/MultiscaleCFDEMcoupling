#include "fileName.H"
#include "cfdemCloudMix.H"
#include "global.H"
#include "forceModel.H"
#include "locateModel.H"
#include "momCoupleModel.H"
#include "meshMotionModel.H"
#include "voidFractionModel.H"
#include "dataExchangeModel.H"
#include "IOModel.H"
#include "probeModel.H"
#include "registryModel.H"
#include "averagingModel.H"
#include "clockModel.H"
#include "smoothingModel.H"
#include "liggghtsCommandModel.H"

#include "mpi.h"
#include "IOmanip.H"
#include "OFversion.H"

namespace Foam {

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

cfdemCloudMix::cfdemCloudMix(const fvMesh& mesh):
  cfdemCloud(mesh),
  angularVelocities_(NULL),
  DEMTorques_(NULL),
  skipLagrangeToEulerMapping_(false),
  skipAfter_(false),
  timeStepsToSkip_(0),
  calculateTortuosity_(false),
  flowDir_(vector()),
  frontMeshRefineField_(NULL),
  frontMeshRefine_(false),
  pRefCell_(0),
  pRefValue_(0),
  haveEvolvedOnce_(false),
  volumefractions_(NULL) {

  Info << "\nEnding of Constructing cfdemCloudMix Class Object......\n" << endl;
  Info << "\nEntry of cfdemCloudMix::cfdemCloudMix......\n" << endl;

  if (usedForSolverPiso() == false) {
    pRefCell_ = readLabel(mesh.solutionDict().subDict("PISO").lookup("pRefCell"));
    pRefValue_ = readScalar(mesh.solutionDict().subDict("PISO").lookup("pRefValue"));

    if (this->couplingProperties().found("skipLagrangeToEulerMapping")) {
      Info << "Will skip lagrange-to-Euler mapping...\n" << endl;
      skipLagrangeToEulerMapping_=true;
    }

    if (this->couplingProperties().found("timeStepsBeforeSkipping")) {
      skipAfter_ = true;
      timeStepsToSkip_ = readScalar(this->couplingProperties().lookup("timeStepsBeforeSkipping"));
      Info << "Will skip LagrangeToEuler mapping after " << timeStepsToSkip_ << " time steps\n" <<  endl;
    }

    if (this->couplingProperties().found("tortuosity")) {
      calculateTortuosity_ = true;
      flowDir_ = this->couplingProperties().subDict("tortuosity").lookup("flowDirection");
      flowDir_ = flowDir_ / mag(flowDir_);
      Info << "Will calculate tortuosity in the mean flow direction "
        << "(" << flowDir_[0] << ", " << flowDir_[1] << ", " << flowDir_[2] << ")" << endl;
    }

    // Must check for walls in case of checkPeriodicCells
    // periodic check will mirror particles and probing points to ensure proper behavior near processor bounds
    if (checkPeriodicCells_) {
      // Enforce reading of the blocking for periodic checks
      // ??????????????????????????????????
      if (readBool(this->couplingProperties().subDict("wall_blockPeriodicityCheck").lookup("x"))) {
        wall_periodicityCheckRange_[0] = 0;
      }
      if (readBool(this->couplingProperties().subDict("wall_blockPeriodicityCheck").lookup("y"))) {
        wall_periodicityCheckRange_[1] = 0;
      }
      if (readBool(this->couplingProperties().subDict("wall_blockPeriodicityCheck").lookup("z"))) {
        wall_periodicityCheckRange_[2] = 0;
      }

      if (this->couplingProperties().found("wall_periodicityCheckTolerance")) {
        wall_periodicityCheckTolerance_ = readScalar(this->couplingProperties().lookup("wall_periodicityCheckTolerance"));
      }
    }
  }  // End of usedForSolverPiso() == false
  Info << "\ncfdemCloudMix::cfdemCloudMix - done\n" << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cfdemCloudMix::~cfdemCloudMix() {
  if (usedForSolverPiso()) {
    return ;
  } else if (usedForSolverIB()) {
    dataExchangeM().destroy(angularVelocities_, 3);
    dataExchangeM().destroy(dragPrev_, 3);
    dataExchangeM().destroy(DEMTorques_, 3);
  } else {
    dataExchangeM().destroy(angularVelocities_, 3);
    dataExchangeM().destroy(dragPrev_, 3);
    dataExchangeM().destroy(DEMTorques_, 3);
    dataExchangeM().destroyDiscreteMemory(cellIDs_, numberOfParticles_);
    dataExchangeM().destroyDiscreteMemory(volumefractions_, numberOfParticles_);
  }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void cfdemCloudMix::getDEMdata() {
  cfdemCloud::getDEMdata();
  if (usedForSolverPiso() == false) {
    dataExchangeM().getData("omega", "vector-atom", angularVelocities_);
  }
}

void cfdemCloudMix::giveDEMdata() {
  cfdemCloud::giveDEMdata();
  if (usedForSolverPiso() == false) {
    dataExchangeM().giveData("hdtorque", "vector-atom", DEMTorques_);
  }
}

void cfdemCloudMix::doCoupleDebug() {
  if (verbose_) {
    int timeIndex = mesh().time().timeIndex();
    int timeIndexOffset = dataExchangeM().timeIndexOffset();
    scalar CFDts = mesh().time().deltaT().value();         // 流体时间步长
    int cfdStep = timeIndex - timeIndexOffset;             // 发生耦合的流体时间步数
    scalar couplingTime = dataExchangeM().couplingTime();  // 耦合时间步长
    int couplingStep = dataExchangeM().couplingStep();     // 耦合时间步数
    Info << "\tcurrent cfdStep: " << cfdStep << endl;
    Info << "\tcurrent couplingStep: " << couplingStep << endl;
    Info << "\tcfdStep * CFDts: " << cfdStep * CFDts - SMALL << endl;
    Info << "\tcouplingStep * couplingTime: " << couplingStep * couplingTime << endl;
  }
}

void cfdemCloudMix::findCellDebug() {
  if (verbose_) {
    for (int index = 0; index < numberOfParticles(); ++index) {
      Info << "\tparticle index: " << index << ", position: " << position(index) << ", cellIDs: "
        << cellIDs()[index][0] << ", cell centre: " << mesh().C()[cellIDs()[index][0]] << endl;
    }
  }
}

// @brief 强制重新分配内存
void cfdemCloudMix::mixForceReAllocArrays() {
  if (numberOfParticlesChanged_ && !arraysReallocated_) {
    size_t fmWidth = std::max(voidFractionM().maxCellsNumPerFineParticle(),
                              voidFractionM().maxCellsNumPerMiddleParticle());
    size_t cWidth = voidFractionM().maxCellsNumPerCoarseParticle();

    dataExchangeM().allocateArray(voidfractions_, 1., fmWidth);
    dataExchangeM().allocateArray(particleWeights_, 0., fmWidth);
    dataExchangeM().allocateArray(particleVolumes_, 0., fmWidth);

    // 分配 cellIDs_ 和 volumefractions_ 内存(double**)
    // 首先将宽度设置为 1, 获取 double* 指针的内存, 同时释放 double 内存
    dataExchangeM().allocateArray(cellIDs_, 0., 1, numberOfParticles_);
    dataExchangeM().allocateArray(volumefractions_, 0., 1, numberOfParticles_);
    for (int i = 0; i < numberOfParticles_; ++i) {
      dataExchangeM().destroy(cellIDs_[i]);
      dataExchangeM().destroy(volumefractions_[i]);
    }
    for (int i = 0; i < numberOfParticles_; ++i) {
      // 获取当前颗粒的尺度
      double globalRatio = globalDimensionRatios_[i];
      double** temp_cellIDs = NULL;
      double** temp_volumefractions = NULL;
      if (checkFAndMParticle(globalRatio)) {
        dataExchangeM().allocateArray(temp_cellIDs, 0, fmWidth, 1);
        dataExchangeM().allocateArray(temp_volumefractions, 0, 1, 1);
      } else {
        dataExchangeM().allocateArray(temp_cellIDs, 0, cWidth, 1);
        dataExchangeM().allocateArray(temp_volumefractions, 0, cWidth, 1);
      }
      cellIDs_[i] = temp_cellIDs[0];
      volumefractions_[i] = temp_volumefractions[0];
    }
  }
  arraysReallocated_ = true;
}

bool cfdemCloudMix::reAllocArrays() const {
  if (numberOfParticlesChanged_ && !arraysReallocated_) {
    // get arrays of new length
    dataExchangeM().allocateArray(positions_, 0., 3);
    dataExchangeM().allocateArray(fluidVel_, 0., 3);
    dataExchangeM().allocateArray(fAcc_, 0., 3);
    dataExchangeM().allocateArray(impForces_, 0., 3);
    dataExchangeM().allocateArray(expForces_, 0., 3);
    dataExchangeM().allocateArray(DEMForces_, 0., 3);
    dataExchangeM().allocateArray(Cds_, 0., 1);
    dataExchangeM().allocateArray(radii_, 0., 1);
    dataExchangeM().allocateArray(particleV_, 0., 1);
    dataExchangeM().allocateArray(velocities_, 0., 3);

    if (usedForSolverPiso()) {
      size_t width = std::max(voidFractionM().maxCellsNumPerFineParticle(),
                              voidFractionM().maxCellsNumPerMiddleParticle());
      dataExchangeM().allocateArray(cellIDs_, -1., width);
      dataExchangeM().allocateArray(voidfractions_, 1., width);
      dataExchangeM().allocateArray(particleWeights_, 0., width);
      dataExchangeM().allocateArray(particleVolumes_, 0., width);
      arraysReallocated_ = true;
    } else if (!usedForSolverPiso() && usedForSolverIB()) {
      dataExchangeM().allocateArray(cellIDs_, -1., voidFractionM().maxCellsPerParticle());
      dataExchangeM().allocateArray(voidfractions_, 1., voidFractionM().maxCellsPerParticle());
      dataExchangeM().allocateArray(particleWeights_, 0., voidFractionM().maxCellsPerParticle());
      dataExchangeM().allocateArray(particleVolumes_, 0., voidFractionM().maxCellsPerParticle());
      dataExchangeM().allocateArray(dragPrev_, 0, 3);
      dataExchangeM().allocateArray(DEMTorques_, 0, 3);
      dataExchangeM().allocateArray(angularVelocities_, 0, 3);
      arraysReallocated_ = true;
    } else {
      dataExchangeM().allocateArray(dragPrev_, 0, 3);
      dataExchangeM().allocateArray(DEMTorques_, 0, 3);
      dataExchangeM().allocateArray(angularVelocities_, 0, 3);
      // 此处不能设置 arraysReallocated_ 为 true
      // 在 mixForceReAllocArrays 中还需要分配内存
    }

    if (namesFieldsUserCFDEMToExt.size() != particleDatFieldsUserCFDEMToExt.size()){
      allocateParticleDatFieldsUserCFDEMToExt();
    } else {
      reAllocateParticleDatFieldsUserCFDEMToExt();
    }
    return true;
  }
  return false;
}

bool cfdemCloudMix::reAllocArrays(int nP, bool forceRealloc) const {
  if ((numberOfParticlesChanged_ && !arraysReallocated_) || forceRealloc) {
    // get arrays of new length
    dataExchangeM().allocateArray(positions_, 0., 3, nP);
    dataExchangeM().allocateArray(velocities_, 0., 3, nP);
    dataExchangeM().allocateArray(fluidVel_, 0., 3, nP);
    dataExchangeM().allocateArray(impForces_, 0., 3, nP);
    dataExchangeM().allocateArray(expForces_, 0., 3, nP);
    dataExchangeM().allocateArray(DEMForces_, 0., 3, nP);
    dataExchangeM().allocateArray(Cds_, 0., 1, nP);
    dataExchangeM().allocateArray(radii_, 0., 1, nP);
    dataExchangeM().allocateArray(particleV_, 0., 1, nP);
    // dataExchangeM().allocateArray(cellIDs_, 0., voidFractionM().maxCellsPerParticle(), nP);

    if (usedForSolverPiso()) {
      size_t width = std::max(voidFractionM().maxCellsNumPerFineParticle(), voidFractionM().maxCellsNumPerMiddleParticle());
      dataExchangeM().allocateArray(cellIDs_, -1., width, nP);
      dataExchangeM().allocateArray(voidfractions_, 1., width, nP);
      dataExchangeM().allocateArray(particleWeights_, 0., width, nP);
      dataExchangeM().allocateArray(particleVolumes_, 0., width, nP);
      arraysReallocated_ = true;
    } else if (!usedForSolverPiso() && usedForSolverIB()) {
      dataExchangeM().allocateArray(cellIDs_, -1., voidFractionM().maxCellsPerParticle(), nP);
      dataExchangeM().allocateArray(voidfractions_, 1., voidFractionM().maxCellsPerParticle(), nP);
      dataExchangeM().allocateArray(particleWeights_, 0., voidFractionM().maxCellsPerParticle(), nP);
      dataExchangeM().allocateArray(particleVolumes_, 0., voidFractionM().maxCellsPerParticle(), nP);
      dataExchangeM().allocateArray(angularVelocities_, 0, 3, nP);
      dataExchangeM().allocateArray(dragPrev_, 0, 3, nP);
      dataExchangeM().allocateArray(DEMTorques_, 0, 3, nP);
      arraysReallocated_ = true;
    } else {
      FatalError << "cfdemCloudMix::reAllocArrays(): not implement when usedForSolverPiso_ == false and usedForSolverIB == false." << abort(FatalError);
    }

    if (namesFieldsUserCFDEMToExt.size() != particleDatFieldsUserCFDEMToExt.size()) {
      allocateParticleDatFieldsUserCFDEMToExt();
    } else {
      reAllocateParticleDatFieldsUserCFDEMToExt();
    }
    return true;
  }
  return false;
}

void cfdemCloudMix::mixResetFieldKernel() {
  // 重置局部平均颗粒速度
  averagingM().UsPrev() == averagingM().UsNext();
  averagingM().UsNext() == dimensionedVector("zero",  averagingM().UsNext().dimensions(), vector::zero);
  Info << "Reset Us fields - done" << endl;

  // 重置空隙率场
  voidFractionM().resetVoidFractions();
  Info << "Reset voidfraction fields - done" << endl;

  // 重置体积分数场
  voidFractionM().resetVolumeFractions();
  Info << "Reset volumefraction fields - done" << endl;

  // 重置隐式力场
  averagingM().resetVectorAverage(forceM(0).impParticleForces(),
                                  forceM(0).impParticleForces(), true);
  Info << "Reset implicit force fields - done" << endl;

  // 重置显式力场
  averagingM().resetVectorAverage(forceM(0).expParticleForces(),
                                  forceM(0).expParticleForces(), true);
  Info << "Reset Explicit force fields - done" << endl;

  // 重置颗粒速度影响因数场
  averagingM().UsWeightField() == dimensionedScalar("zero", averagingM().UsWeightField().dimensions(), 0.0);
  Info << "Reset Us weight fields - done" << endl;

  // 重置单位体积动量交换场
  for (int i = 0; i < momCoupleModels_.size(); i++) {
    // momCoupleM(i).KslPrev_ == momCoupleM(i).KslNext_;
    // momCoupleM(i).KslNext_ ==  dimensionedScalar("zero", momCoupleM(i).KslNext_.dimensions(), 0.0);
    momCoupleM(i).resetMomSourceField();
  }
  Info << "Reset Ksl fields - done" << endl;
}

/*! \brief 处理流场与颗粒场之间的耦合 */
void cfdemCloudMix::mixCouplingKernel() {
  // 获取 DEM 数据(颗粒半径, 颗粒位置矢量, 颗粒速度...)
  if (verbose_) { Info << "cfdemCloudMix::getDEMdata()..." << endl; }
  getDEMdata();
  if (verbose_) { Info << "cfdemCloudMix::getDEMdata() - done\n" << endl; }

  // 计算颗粒尺寸与其周围网格平均尺寸的比值, 并将颗粒索引按照颗粒尺寸归类
  if (verbose_) { Info << "voidFractionModel::getDimensionRatios()..." << endl; }
  std::vector<double> temp_vec(numberOfParticles_, -1);
  std::vector<double*> cellIDs_vec(numberOfParticles_, NULL);
  for (int i = 0; i < numberOfParticles_; ++i) { cellIDs_vec[i] = temp_vec.data() + i; }
  double** cellIDs_ptr = cellIDs_vec.data();
  locateM().findCell(NULL, positions_, cellIDs_ptr, numberOfParticles());  // 搜索颗粒中心所在网格编号
  const_cast<voidFractionModel&>(voidFractionM()).getDimensionRatios(dimensionRatios_,
                                                                     globalDimensionRatios_,
                                                                     sumCellsNumbers_,
                                                                     sumCellsVolumes_,
                                                                     cellIDs_vec);
  if (verbose_) { Info << "voidFractionModel::getDimensionRatios() - done\n" << endl; }

  // 根据颗粒尺度分配内存
  if (verbose_) { Info << "cfdemCloudMix::mixForceReAllocArrays()..." << endl; }
  mixForceReAllocArrays();
  if (verbose_) { Info << "cfdemCloudMix::mixForceReAllocArrays() - done\n" << endl; }

  // 定位颗粒中心所在网格的编号(创建临时空间存放颗粒中心所在网格的编号)
  if (verbose_) { Info << "cfdemCloudMix::findCell()..." << endl; }
  findCells();
  if (verbose_) { Info << "cfdemCloudMix::findCell() - done\n" << endl; }

  // 设置小颗粒空隙率场
  if (verbose_) Info << "cfdemCloudMix::setvoidFraction()..." << endl;
  voidFractionM().voidFractionModelInit(particleWeights_,
                                        particleVolumes_,
                                        particleV_,
                                        dimensionRatios_);
  voidFractionM().setVoidFractionAndVolumeFraction(NULL,
                                                   volumefractions_,
                                                   voidfractions_,
                                                   particleWeights_,
                                                   particleVolumes_,
                                                   particleV_,
                                                   dimensionRatios_);
  if (verbose_) Info << "cfdemCloudMix::setvoidFraction() - done\n" << endl;

  // 设置局部平均颗粒速度场，将小颗粒的速度 velocities_ 映射到局部平局小颗粒速度，本质上是将拉格朗日场映射到欧拉场
  // setVectorAverages();
  if (verbose_) Info << "cfdemCloudMix::setMixVectorAverage()..." << endl;
  averagingM().setMixVectorAverage(velocities_,
                                   particleWeights_,
                                   dimensionRatios_,
                                   averagingM().UsNext(),
                                   averagingM().UsWeightField(),
                                   NULL);
  if (verbose_) Info << "cfdemCloudMix::setMixVectorAverage() - done\n" << endl;

  // smoothen "next" fields
  if (verbose_) Info << "averageModel::smoothen()..." << endl;
  smoothingM().dSmoothing();
  smoothingM().smoothen(voidFractionM().voidFractionNext());
  // only smoothen if we use implicit force coupling in cells void of particles
  // because we need unsmoothened Us field to detect cells for explicit force coupling
  if (!treatVoidCellsAsExplicitForce()) {
    smoothingM().smoothenReferenceField(averagingM().UsNext());
  }
  if (verbose_) Info << "averageModel::smoothen() - done\n" << endl;
  MPI_Barrier(MPI_COMM_WORLD);
}

/*!
 * \brief 更新函数
 * \param alphaVoidfraction    <[in, out] 小颗粒空隙率场
 * \param alphaVolumefraction  <[in, out] 大颗粒体积分数场
 * \param Us                   <[in, out] 局部平均小颗粒速度场
 * \param U                    <[in] 流体速度场
 * \param interFace            <[in, out] 界面场
 */
bool cfdemCloudMix::evolve(volScalarField& alphaVoidfraction,
                           volScalarField& alphaVolumefraction,
                           volVectorField& Us,
                           volVectorField& U,
                           volScalarField& interFace) {
  Info << "\nFoam::cfdemCloudMix::evolve()...\n" << endl;
  numberOfParticlesChanged_ = false;
  arraysReallocated_ = false;
  bool doCouple = false;

  if (!ignore()) {

    // 判断是否跳过 mapping 步骤
    if (skipAfter_ && timeStepsToSkip_ < 1) {
      skipLagrangeToEulerMapping_ = true;
    }

    if (!writeTimePassed_ && mesh_.time().outputTime()) {
      writeTimePassed_ = true;
    }

    doCoupleDebug();
    // 因为耦合时间步长 = 流体时间步长的整数倍，所以判断当前流体时间步是否同时也是耦合时间步
    if (dataExchangeM().doCoupleNow()) {
      Info << "do coupling now...\n" << endl;

      // 当 dataExchangeModel 为 twoWayMPI or mixTwoWayMPI 时, couple 函数从 DEM 求解器读取颗粒数量, 同时调用 reAllocArray 函数分配内存
      // 注意: 在 reAllocArrays 函数中针对所有尺度的颗粒分配内存
      doCouple = dataExchangeM().couple(0);

      if (verbose_) {
        Info << "number of particles: " << numberOfParticles() << endl;                      // 颗粒数量
        Info << "current couplingStep: " << dataExchangeM().couplingStep() << endl << endl;  // 当前耦合时间步
      }

      // 重新设置所有的 field
      Info << "ResetVolFields now..." << endl;
      mixResetFieldKernel();
      Info << "ResetVolFields - done\n" << endl;

      // do coupling
      Info << "mixCouplingKernel now..." << endl;
      mixCouplingKernel();
      Info << "mixCouplingKernel - done\n" << endl;

      if (verbose_) { Info << "cfdemCloudMix::setInterFace()..." << endl; }
      setInterFace(interFace);
      if (verbose_) { Info << "cfdemCloudMix::setInterFace() - done\n" << endl; }

      Info << "Do coupling - done\n" << endl;
    }  // End of doCoupleNow

    // 获取 timeStepFraction，以防止使用自适应流体时间步
    if (dataExchangeM().timeStepFraction() > 1.001 ||
        dataExchangeM().timeStepFraction() < -0.001) {
      FatalError << "cfdemCloud::dataExchangeM().timeStepFraction() = "<< dataExchangeM().timeStepFraction()
        << " must be >1 or <0 : This might be due to the fact that you used a adjustable CFD time step. Please use a fixed CFD time step."
        << abort(FatalError);
    }

    // 更新 voidFractionField, volumeFractionField
    alphaVoidfraction = voidFractionM().voidFractionNext();
    alphaVolumefraction = voidFractionM().volumeFractionNext();
    if (dataExchangeM().couplingStep() < 2) {
      alphaVoidfraction.oldTime() = alphaVoidfraction; // supress volume src
      alphaVoidfraction.oldTime().correctBoundaryConditions();
    }
    alphaVoidfraction.correctBoundaryConditions();

    // calc ddt(voidfraction)
    calcDdtVoidfraction(alphaVoidfraction, Us);

    // 更新局部平均小颗粒速度场
    Us = averagingM().UsInterp();
    Us.correctBoundaryConditions();

    if (doCouple) {  // 如果当前时间步是耦合时间步
      if (verbose_) { Info << "cfdemCloudMix::setForce()..." << endl; }
      // setForces();
      resetArray(impForces_, numberOfParticles_, 3, 0.0);  // 重置颗粒对流体的隐式力
      resetArray(expForces_, numberOfParticles_, 3, 0.0);  // 重置颗粒对流体的显式力
      resetArray(DEMForces_, numberOfParticles_, 3, 0.0);  // 重置流体对颗粒的总作用力
      resetArray(fluidVel_, numberOfParticles_, 3, 0.0);   // 重置小颗粒所在网格的流体速度
      resetArray(Cds_, numberOfParticles_, 1, 0.0);        // 重置阻力系数
      //reset all USER-defined particle fields
      zeroizeParticleDatFieldsUserCFDEMToExt();
      for (int i = 0; i < cfdemCloud::nrForceModels(); i++) {
        cfdemCloud::forceM(i).setMixForce(dimensionRatios_);
      }
      if (verbose_) { Info << "cfdemCloudMix::setForce() - done\n" << endl; }

      if (verbose_) { Info << "cfdemCloudMix::calcMultiphaseTurbulence()..." << endl; }
      calcMultiphaseTurbulence();
      if (verbose_) { Info << "cfdemCloudMix::calcMultiphaseTurbulence() - done\n" << endl; }

      // get next force field
      if (verbose_) { Info << "cfdemCloudMix::setParticleForceField()..." << endl; }
      // setParticleForceField();
      averagingM().setVectorSum(forceM(0).impParticleForces(), impForces_,
                                particleWeights_, dimensionRatios_);
      averagingM().setVectorSum(forceM(0).expParticleForces(), expForces_,
                                particleWeights_, dimensionRatios_);
      if (verbose_) { Info << "cfdemCloudMix::setParticleForceField() - done\n" << endl; }

      // write DEM data
      if (verbose_) { Info << "cfdemCloudMix::giveDEMdata()..." << endl; }
      giveDEMdata();
      if (verbose_) { Info << "cfdemCloudMix::giveDEMdata() - done\n" << endl; }

      dataExchangeM().couple(1);
      haveEvolvedOnce_ = true;
    }

    // do particle IO
    IOM().dumpDEMdata();
    if (skipAfter_) {
      timeStepsToSkip_ -= 1;
      Info << "Will skip LagrangeToEuler mapping after " << timeStepsToSkip_ << " time steps" <<  endl;
    }
  }  // End of ignore
  Info << "Foam::cfdemCloudMix::evolve() - done\n" << endl;
  return doCouple;
}

// @brief 更新函数
// @param alpha  <[in, out] 小颗粒空隙率场
// @param Us     <[in, out] 局部平均小颗粒速度场
// @param U      <[in] 流体速度场
// @note used for cfdemSolverPiso
bool cfdemCloudMix::evolve(volScalarField& alpha,
                           volVectorField& Us,
                           volVectorField& U) {
  Info << "\n// * * * * * * * * * * * * * * * " << "Foam::cfdemCloudMix::evolve()" << " * * * * * * * * * * * * * * * //\n" << endl;
  numberOfParticlesChanged_ = false;
  arraysReallocated_ = false;
  bool doCouple = false;

  if (!ignore()) {
    if (!writeTimePassed_ && mesh_.time().outputTime()) {
      writeTimePassed_ = true;
    }
    // 因为耦合时间步长 = 流体时间步长的整数倍，所以判断当前流体时间步是否同时也是耦合时间步
    doCoupleDebug();
    if (dataExchangeM().doCoupleNow()) {
      Info << "\n// * * * * * * * * * * " << "Coupling" << " * * * * * * * * * * //\n" << endl;
      // 当 dataExchangeModel 为 twoWayMPI 时, couple(0) 函数从 DEM 求解器读取颗粒数量, 同时调用 reAllocArray 函数分配内存
      dataExchangeM().couple(0);
      Info << "number of particles = " << numberOfParticles() << "\n" << endl;
      doCouple = true;

      if (verbose_) {
        Info << "\n// * * * * * * * * * * " << "ResetVolFields" << " * * * * * * * * * * //\n" << endl;
        Info << "number of particles = " << numberOfParticles() << endl;
        // dataExchangeM().couplingStep() 返回当前耦合时间步
        Info << "current couplingStep: " << dataExchangeM().couplingStep() << endl;
      }

      // 重置局部平均颗粒速度
      // averagingM().UsPrev() == averagingM().UsNext();
      // averagingM().UsNext() == dimensionedVector("zero",  averagingM().UsNext().dimensions(), vector::zero);
      averagingM().resetVectorAverage(averagingM().UsPrev(),
                                      averagingM().UsNext(), false);

      // 重置小颗粒空隙率场
      // voidfractionPrev_ == voidfractionNext_;
      // voidfractionNext_ == dimensionedScalar("one", voidfractionNext_.dimensions(), 1.);
      resetVoidFraction();

      // 重置隐式力场
      // forceM(0).impParticleForces() ==
      //   dimensionedVector("zero", forceM(0).impParticleForces().dimensions(), vector::zero);
      averagingM().resetVectorAverage(forceM(0).impParticleForces(),
                                      forceM(0).impParticleForces(), true);

      // 重置显式力场
      // forceM(0).expParticleForces() ==
      //   dimensionedVector("zero", forceM(0).expParticleForces().dimensions(), vector::zero);
      averagingM().resetVectorAverage(forceM(0).expParticleForces(),
                                      forceM(0).expParticleForces(), true);

      // 重置颗粒速度影响因数场
      // averagingM().resetWeightFields();
      averagingM().UsWeightField() == dimensionedScalar("zero", averagingM().UsWeightField().dimensions(), 0.0);

      // 重置动量交换场
      for (int i = 0; i < momCoupleModels_.size(); i++) {
        // momCoupleM(i).KslPrev_ == momCoupleM(i).KslNext_;
        // momCoupleM(i).KslNext_ ==  dimensionedScalar("zero", momCoupleM(i).KslNext_.dimensions(),  0.0);
        momCoupleM(i).resetMomSourceField();
      }
      if (verbose_) {
        Info << "\n// * * * * * * * * * * " << "ResetVolFields - done" << " * * * * * * * * * * //\n" << endl;
      }

      // End of reset fields, start of evolve fields!

      // 获取 DEM 数据
      if (verbose_) { Info << "cfdemCloudMix::getDEMdata()..." << endl; }
      getDEMdata();
      if (verbose_) { Info << "cfdemCloudMix::getDEMdata() - done\n" << endl; }

      // search cellIDs of particles
      if (verbose_) { Info << "cfdemCloudMix::findCell()..." << endl; }
      // locateM().findCell(NULL, positions_, cellIDs_, numberOfParticles());
      findCells();
      if (verbose_) { Info << "cfdemCloudMix::findCell() - done\n" << endl; }

      // // 计算颗粒尺寸与其周围网格平均尺寸的比值, 并将颗粒索引按照颗粒尺寸归类
      // if (verbose_) { Info << "\nvoidFractionModel::getDimensionRatios()..." << endl; }
      // const_cast<voidFractionModel&>(voidFractionM()).getDimensionRatios(dimensionRatios_,
      //                                                                     globalDimensionRatios_,
      //                                                                     sumCellsNumbers_,
      //                                                                     sumCellsVolumes_);
      // if (verbose_) { Info << "voidFractionModel::getDimensionRatios() - done" << endl; }

      // 设置小颗粒空隙率场
      if (verbose_) Info << "cfdemCloudMix::setvoidFraction()..." << endl;
      voidFractionM().setvoidFraction(NULL,
                                      voidfractions_,
                                      particleWeights_,
                                      particleVolumes_,
                                      particleV_);
      if (verbose_) Info << "cfdemCloudMix::setvoidFraction() - done\n" << endl;

      // 设置局部平均颗粒速度场，将小颗粒的速度 velocities_ 映射到局部平局小颗粒速度，本质上是将拉格朗日场映射到欧拉场
      // setVectorAverages();
      if (verbose_) Info << "cfdemCloudMix::setVectorAverage()..." << endl;
      averagingM().setVectorAverage(averagingM().UsNext(),
                                    velocities_,
                                    particleWeights_,
                                    averagingM().UsWeightField(),
                                    NULL,
                                    NULL,
                                    false);
      if (verbose_) Info << "cfdemCloudMix::setVectorAverage() - done\n" << endl;

      // smoothen "next" fields
      smoothingM().dSmoothing();
      smoothingM().smoothen(voidFractionM().voidFractionNext());
      // only smoothen if we use implicit force coupling in cells void of particles
      // because we need unsmoothened Us field to detect cells for explicit force coupling
      if (!treatVoidCellsAsExplicitForce()) {
        smoothingM().smoothenReferenceField(averagingM().UsNext());
      }
      Info << "\n// * * * * * * * * * * " << "Coupling - done" << " * * * * * * * * * * //\n" << endl;
    }  // End of doCoupleNow

    // 获取 timeStepFraction，以防止使用自适应流体时间步
    if (verbose_) { Info << "timeStepFraction() = " << dataExchangeM().timeStepFraction() << endl; }
    if (dataExchangeM().timeStepFraction() > 1.001 ||
       dataExchangeM().timeStepFraction() < -0.001) {
      FatalError << "cfdemCloud::dataExchangeM().timeStepFraction() = "<< dataExchangeM().timeStepFraction()
        << " must be >1 or <0 : This might be due to the fact that you used a adjustable CFD time step. Please use a fixed CFD time step."
        << abort(FatalError);
    }

    // update voidFractionField
    setAlpha(alpha);
    if (dataExchangeM().couplingStep() < 2) {
      alpha.oldTime() = alpha; // supress volume src
      alpha.oldTime().correctBoundaryConditions();
    }
    alpha.correctBoundaryConditions();

    // calc ddt(voidfraction)
    calcDdtVoidfraction(alpha, Us);

    // update mean particle velocity Field
    Us = averagingM().UsInterp();
    Us.correctBoundaryConditions();

    if (doCouple) {
      // set particles forces
      if (verbose_) { Info << "cfdemCloudMix::setForce()..." << endl; }
      setForces();
      if (verbose_) { Info << "cfdemCloudMix::setForce() - done\n" << endl; }

      if (verbose_) { Info << "cfdemCloudMix::ncalcMultiphaseTurbulence()..." << endl; }
      calcMultiphaseTurbulence();
      if (verbose_) { Info << "cfdemCloudMix::calcMultiphaseTurbulence() - done\n" << endl; }

      // get next force field
      if (verbose_) { Info << "cfdemCloudMix::nsetParticleForceField()..." << endl; }
      setParticleForceField();
      if (verbose_) { Info << "cfdemCloudMix::setParticleForceField() - done\n" << endl; }

      // write DEM data
      if (verbose_) { Info << "cfdemCloudMix::ngiveDEMdata()..." << endl; }
      giveDEMdata();
      if (verbose_) { Info << "cfdemCloudMix::giveDEMdata() - done\n" << endl; }

      dataExchangeM().couple(1);
    }
    if (verbose_) {
      // #include "debugInfo.H"
    }
    // do particle IO
    IOM().dumpDEMdata();
  }  // End of ignore
  Info << "\n// * * * * * * * * * * * * * * * " << "Foam::cfdemCloudMix::evolve() - done" << " * * * * * * * * * * * * * * * //\n" << endl;
  return doCouple;
}

void cfdemCloudMix::calcVelocityCorrection(volScalarField& p,
                                           volVectorField& U,
                                           volScalarField& phiIB,
                                           volScalarField& voidfraction) {
  setParticleVelocity(U);

  // make field divergence free - set reference value in case it is needed
  fvScalarMatrix phiIBEqn(
    fvm::laplacian(phiIB) == fvc::div(U) + fvc::ddt(voidfraction)
  );

  if(phiIB.needReference()) {
    phiIBEqn.setReference(pRefCell_, pRefValue_);
  }
  
  phiIBEqn.solve();

  U = U - fvc::grad(phiIB);
  U.correctBoundaryConditions();

  // correct the pressure as well
  // do we have to account for rho here?
  p = p + phiIB / U.mesh().time().deltaT();
  p.correctBoundaryConditions();

  if (couplingProperties_.found("checkinterface")) {
    Info << "checking no-slip on interface..." << endl;
    // #include "checkInterfaceVelocity.H"
  }
}

// @brief 设置 coarse particle 中的流体速度
void cfdemCloudMix::setParticleVelocity(volVectorField& U) {
  label cellI = 0;
  vector particleVel(0,0,0);
  vector relativeVec(0, 0, 0);
  vector rotationVel(0, 0, 0);
  vector angularVel(0, 0, 0);

  for (int index = 0; index < numberOfParticles(); index++) {

    bool particleNeedSet = false;
    if (usedForSolverIB()) {
      particleNeedSet = true;
    } else if (dimensionRatios_.empty() != false) {
      particleNeedSet = checkCoarseParticle(dimensionRatios_[index]);
    }

    if (particleNeedSet) {
      for (int subCell = 0; subCell < cellsPerParticle()[index][0]; ++subCell) {
        cellI = cellIDs()[index][subCell];

        if (cellI >= 0) {
          for (int i = 0;i < 3; i++) {
            relativeVec[i] = U.mesh().C()[cellI][i] - position(index)[i];
          }
          for (int i = 0; i < 3; i++) {
            angularVel[i] = angularVelocities()[index][i];
          }
          // 计算转动速度
          rotationVel = angularVel ^ relativeVec;
          // 计算颗粒速度
          for (int i = 0; i < 3; ++i) {
            particleVel[i] = velocities()[index][i] + rotationVel[i];
          }
          // 设置网格速度
          if (usedForSolverIB()) {
            U[cellI] = (1 - voidfractions_[index][subCell]) * particleVel + voidfractions_[index][subCell] * U[cellI];
          } else {
            // ??????????????????????????????????????????
          }
        }  // End of cellI >= 0
      }  // End of loop subCells
    }  // End of particleNeedSet
  }  // End of loop all particles
  U.correctBoundaryConditions();
}

// @brief 确定颗粒周围细化网格的区域(每个方向的尺寸都是颗粒尺寸的两倍)
void cfdemCloudMix::setInterFace(volScalarField& interFace) {

  interFace == dimensionedScalar("zero", interFace.dimensions(), 0.);

  const boundBox& globalBb = mesh().bounds();

  for (int index = 0; index < numberOfParticles(); index++) {

    bool particleNeedSet = false;

    if (usedForSolverIB()) {
      particleNeedSet = true;
    } else if (dimensionRatios_.empty() != false) {
      particleNeedSet = checkCoarseParticle(dimensionRatios_[index]);
    }

    if (particleNeedSet) {
      vector particlePos = position(index);
      double skin = 2.0;
      // 遍历当前处理器上的所有网格
      forAll (mesh_.C(), cellI) {
        vector cellPos = mesh_.C()[cellI];
        if (checkPeriodicCells_) {

          // Some cells may be located on the other side of a periodic boundary.
          // In this case, the particle center has to be mirrored in order to correctly
          // evaluate the interpolation points.
          vector minPeriodicParticlePos = particlePos;
          voidFractionM().minPeriodicDistance(index,
                                              cellPos,
                                              particlePos,
                                              globalBb,
                                              minPeriodicParticlePos,
                                              wall_periodicityCheckRange());
          particlePos = minPeriodicParticlePos;
        }
        double value = voidFractionM().pointInParticle(index, particlePos, cellPos, skin);
        if (value <= 0.0) {
          interFace[cellI] = value + 1.0;
        }
      }  // End of loop all cells
    }  // End of particleNeedSet
  }  // End of loop all particles
}

// @brief 更新函数
// @param alpha      <[in, out] 大颗粒体积分数场
// @param interFace  <[in, out]
bool cfdemCloudMix::evolve(volScalarField& alpha,
                           volScalarField& interFace) {

  Info << "\n// * * * * * * * * * * * * * * * " << "Foam::cfdemCloudMix::evolve()" << " * * * * * * * * * * * * * * * //\n" << endl;
  numberOfParticlesChanged_ = false;
  arraysReallocated_ = false;
  bool doCouple = false;

  if (!ignore()) {
    if (skipAfter_ && timeStepsToSkip_ < 1) {
      skipLagrangeToEulerMapping_ = true;
    }

    if (!writeTimePassed_ && mesh_.time().outputTime()) {
      writeTimePassed_ = true;
    }

    doCoupleDebug();
    if (dataExchangeM().doCoupleNow()) {
      Info << "\n// * * * * * * * * * * " << "Coupling" << " * * * * * * * * * * //\n" << endl;

      dataExchangeM().couple(0);
      doCouple = true;

      if (!skipLagrangeToEulerMapping_ || !haveEvolvedOnce_) {

        // 获取 DEM 数据
        if (verbose_) { Info << "cfdemCloudMix::getDEMdata()..." << endl; }
        getDEMdata();
        if (verbose_) { Info << "cfdemCloudMix::getDEMdata() - done\n" << endl; }

        // search cellIDs of particles
        if (verbose_) { Info << "cfdemCloudMix::findCell()..." << endl; }
        findCells();
        // locateM().findCell(NULL, positions_, cellIDs_, numberOfParticles());
        if (verbose_) { Info << "cfdemCloudMix::findCell() - done\n" << endl; }

        // // 计算颗粒尺寸与其周围网格平均尺寸的比值, 并将颗粒索引按照颗粒尺寸归类
        // if (verbose_) { Info << "\nvoidFractionModel::getDimensionRatios()..." << endl; }
        // const_cast<voidFractionModel&>(voidFractionM()).getDimensionRatios(dimensionRatios_,
        //                                                                    globalDimensionRatios_,
        //                                                                    sumCellsNumbers_,
        //                                                                    sumCellsVolumes_);
        // if (verbose_) { Info << "voidFractionModel::getDimensionRatios() - done" << endl; }

        // set void fraction field
        if (verbose_) { Info << "cfdemCloudMix::setvoidFraction()..." << endl; }
        voidFractionM().setvoidFraction(NULL,
                                        voidfractions_,
                                        particleWeights_,
                                        particleVolumes_,
                                        particleV_);
        if (verbose_) { Info << "cfdemCloudMix::setvoidFraction() - done\n" << endl; }

        if (verbose_) { Info << "cfdemCloudMix::setInterFace()..." << endl; }
        setInterFace(interFace);
        if (verbose_) { Info << "cfdemCloudMix::setInterFace() - done\n" << endl; }
      }

      // update voidFractionField
      alpha == voidFractionM().voidFractionNext();
      alpha.correctBoundaryConditions();

      // set particles forces
      for (int index = 0; index <  numberOfParticles_; ++index) {
        for (int i = 0; i < 3; i++) {
          impForces_[index][i] = 0;
          expForces_[index][i] = 0;
          DEMForces_[index][i] = 0;
        }
      }
      for (int i = 0; i < cfdemCloud::nrForceModels(); i++) {
        cfdemCloud::forceM(i).setForce();
      }
      // write DEM data
      giveDEMdata();

      dataExchangeM().couple(1);
      haveEvolvedOnce_ = true;

      Info << "\n// * * * * * * * * * * " << "Coupling - done" << " * * * * * * * * * * //\n" << endl;
    }  // End of doCoupleNow

    IOM().dumpDEMdata();
    if (skipAfter_) {
      timeStepsToSkip_ -= 1;
      Info << "Will skip LagrangeToEuler mapping after " << timeStepsToSkip_ << " time steps" <<  endl;
    }
  }  // End of ignore
  Info << "\n// * * * * * * * * * * * * * * * " << "Foam::cfdemCloudMix::evolve() - done" << " * * * * * * * * * * * * * * * //\n" << endl;
  return doCouple;
}

}  // End of namespace Foam
