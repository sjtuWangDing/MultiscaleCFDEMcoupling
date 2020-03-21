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
    // dataExchangeM().destroy(volumefractions_, 1);
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

bool cfdemCloudMix::reAllocArrays() const {
  if (numberOfParticlesChanged_ && !arraysReallocated_) {
    // get arrays of new length
    dataExchangeM().allocateArray(positions_, 0., 3);
    dataExchangeM().allocateArray(velocities_, 0., 3);
    dataExchangeM().allocateArray(fluidVel_, 0., 3);
    dataExchangeM().allocateArray(fAcc_, 0., 3);
    dataExchangeM().allocateArray(impForces_, 0., 3);
    dataExchangeM().allocateArray(expForces_, 0., 3);
    dataExchangeM().allocateArray(DEMForces_, 0., 3);
    dataExchangeM().allocateArray(Cds_, 0., 1);
    dataExchangeM().allocateArray(radii_, 0., 1);
    dataExchangeM().allocateArray(particleV_, 0., 1);
    // dataExchangeM().allocateArray(cellIDs_, -1., voidFractionM().maxCellsPerParticle());

    if (usedForSolverPiso()) {
      size_t width = std::max(voidFractionM().maxCellsNumPerFineParticle(), voidFractionM().maxCellsNumPerMiddleParticle());
      dataExchangeM().allocateArray(cellIDs_, -1., width);
      dataExchangeM().allocateArray(voidfractions_, 1., width);
      dataExchangeM().allocateArray(particleWeights_, 0., width);
      dataExchangeM().allocateArray(particleVolumes_, 0., width);
    }

    if (usedForSolverIB()) {
      dataExchangeM().allocateArray(cellIDs_, -1., voidFractionM().maxCellsPerParticle());
      dataExchangeM().allocateArray(voidfractions_, 1., voidFractionM().maxCellsPerParticle());
      dataExchangeM().allocateArray(particleWeights_, 0., voidFractionM().maxCellsPerParticle());
      dataExchangeM().allocateArray(particleVolumes_, 0., voidFractionM().maxCellsPerParticle());
      dataExchangeM().allocateArray(angularVelocities_, 0, 3);
      dataExchangeM().allocateArray(dragPrev_, 0, 3);
      dataExchangeM().allocateArray(DEMTorques_, 0, 3);
    }

    // Though following arraies only used for coasrse particles, we allocate array for every singal paritcle
    // dataExchangeM().allocateArray(angularVelocities_, 0, 3);
    // dataExchangeM().allocateArray(dragPrev_, 0, 3);
    // dataExchangeM().allocateArray(DEMTorques_, 0, 3);

    if (namesFieldsUserCFDEMToExt.size() != particleDatFieldsUserCFDEMToExt.size()){
      allocateParticleDatFieldsUserCFDEMToExt();
    } else {
      reAllocateParticleDatFieldsUserCFDEMToExt();
    }
    arraysReallocated_ = true;
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
    }

    if (usedForSolverIB()) {
      dataExchangeM().allocateArray(cellIDs_, -1., voidFractionM().maxCellsPerParticle(), nP);
      dataExchangeM().allocateArray(voidfractions_, 1., voidFractionM().maxCellsPerParticle(), nP);
      dataExchangeM().allocateArray(particleWeights_, 0., voidFractionM().maxCellsPerParticle(), nP);
      dataExchangeM().allocateArray(particleVolumes_, 0., voidFractionM().maxCellsPerParticle(), nP);
      dataExchangeM().allocateArray(angularVelocities_, 0, 3, nP);
      dataExchangeM().allocateArray(dragPrev_, 0, 3, nP);
      dataExchangeM().allocateArray(DEMTorques_, 0, 3, nP);
    }

    // Though following arraies only used for coasrse particles, we allocate array for every singal paritcle
    // dataExchangeM().allocateArray(angularVelocities_, 0, 3, nP);
    // dataExchangeM().allocateArray(dragPrev_, 0, 3, nP);
    // dataExchangeM().allocateArray(DEMTorques_, 0, 3, nP);

    if (namesFieldsUserCFDEMToExt.size() != particleDatFieldsUserCFDEMToExt.size()) {
      allocateParticleDatFieldsUserCFDEMToExt();
    } else {
      reAllocateParticleDatFieldsUserCFDEMToExt();
    }
    arraysReallocated_ = true;
    return true;
  }
  return false;
}

void cfdemCloudMix::doCoupleDebug() {
  if (verbose_) {
    int timeIndex = mesh().time().timeIndex();
    int timeIndexOffset = dataExchangeM().timeIndexOffset();
    scalar CFDts = mesh().time().deltaT().value();         // 流体时间步长
    int cfdStep = timeIndex - timeIndexOffset;             // 发生耦合的流体时间步数
    scalar couplingTime = dataExchangeM().couplingTime();  // 耦合时间步长
    int couplingStep = dataExchangeM().couplingStep();     // 耦合时间步数
    Info << "current cfdStep: " << cfdStep << endl;
    Info << "current couplingStep: " << couplingStep << endl;
    Info << "cfdStep * CFDts: " << cfdStep * CFDts - SMALL << endl;
    Info << "couplingStep * couplingTime: " << couplingStep * couplingTime << endl;
  }
}

void cfdemCloudMix::findCellDebug() {
  if (verbose_) {
    for (int index = 0; index < numberOfParticles(); ++index) {
      Info << "particle index: " << index << ", position: " << position(index) << ", cellIDs: "
        << cellIDs()[index][0] << ", cell centre: " << mesh().C()[cellIDs()[index][0]] << endl;
    }
  }
}

// // @brief 更新欧拉场
// // @param alpha  <[in, out] 小颗粒空隙率场
// // @param Us     <[in, out] 局部平均小颗粒速度场
// // @param U      <[in] 流体速度场
// bool cfdemCloudMix::evolve(volScalarField& alpha,
//                            volVectorField& Us,
//                            volVectorField& U) {
//   Info << "\n// * * * * * * * * * * * * * * * " << "Foam::cfdemCloudMix::evolve()" << " * * * * * * * * * * * * * * * //\n" << endl;
//   numberOfParticlesChanged_ = false;
//   arraysReallocated_ = false;
//   bool doCouple = false;

//   if (!ignore()) {
//     if (!writeTimePassed_ && mesh_.time().outputTime()) {
//       writeTimePassed_ = true;
//     }
//     // 因为耦合时间步长 = 流体时间步长的整数倍，所以判断当前流体时间步是否同时也是耦合时间步
//     doCoupleDebug();
//     if (dataExchangeM().doCoupleNow()) {
//       Info << "\n// * * * * * * * * * * " << "Coupling" << " * * * * * * * * * * //\n" << endl;
//       // 当 dataExchangeModel 为 twoWayMPI 时, couple(0) 函数从 DEM 求解器读取颗粒数量, 同时调用 reAllocArray 函数分配内存
//       Info << "number of particles = " << numberOfParticles() << endl;
//       dataExchangeM().couple(0);
//       Info << "number of particles = " << numberOfParticles() << endl;
//       doCouple = true;

//       if (verbose_) {
//         Info << "\n// * * * * * * * * * * " << "ResetVolFields" << " * * * * * * * * * * //\n" << endl;
//         Info << "number of particles = " << numberOfParticles() << endl;
//         // dataExchangeM().couplingStep() 返回当前耦合时间步
//         Info << "current couplingStep: " << dataExchangeM().couplingStep() << endl;
//       }

//       // 重置局部平均颗粒速度
//       // averagingM().UsPrev() == averagingM().UsNext();
//       // averagingM().UsNext() == dimensionedVector("zero",  averagingM().UsNext().dimensions(), vector::zero);
//       averagingM().resetVectorAverage(averagingM().UsPrev(),
//                                       averagingM().UsNext(), false);

//       // 重置小颗粒空隙率场
//       // voidfractionPrev_ == voidfractionNext_;
//       // voidfractionNext_ == dimensionedScalar("one", voidfractionNext_.dimensions(), 1.);
//       resetVoidFraction();

//       // 重置隐式力场
//       // forceM(0).impParticleForces() ==
//       //   dimensionedVector("zero", forceM(0).impParticleForces().dimensions(), vector::zero);
//       averagingM().resetVectorAverage(forceM(0).impParticleForces(),
//                                       forceM(0).impParticleForces(), true);

//       // 重置显式力场
//       // forceM(0).expParticleForces() ==
//       //   dimensionedVector("zero", forceM(0).expParticleForces().dimensions(), vector::zero);
//       averagingM().resetVectorAverage(forceM(0).expParticleForces(),
//                                       forceM(0).expParticleForces(), true);

//       // 重置颗粒速度影响因数场
//       // averagingM().resetWeightFields();
//       averagingM().UsWeightField() == dimensionedScalar("zero", averagingM().UsWeightField().dimensions(), 0.0);

//       // 重置动量交换场
//       for (int i = 0; i < momCoupleModels_.size(); i++) {
//         // momCoupleM(i).KslPrev_ == momCoupleM(i).KslNext_;
//         // momCoupleM(i).KslNext_ ==  dimensionedScalar("zero", momCoupleM(i).KslNext_.dimensions(),  0.0);
//         momCoupleM(i).resetMomSourceField();
//       }
//       if (verbose_) {
//         Info << "\n// * * * * * * * * * * " << "ResetVolFields - done" << " * * * * * * * * * * //\n" << endl;
//       }

//       // End of reset fields, start of evolve fields!

//       // 获取 DEM 数据
//       if (verbose_) { Info << "cfdemCloudMix::getDEMdata()..." << endl; }
//       getDEMdata();
//       if (verbose_) { Info << "cfdemCloudMix::getDEMdata() - done\n" << endl; }

//       // search cellIDs of particles
//       if (verbose_) { Info << "cfdemCloudMix::findCell()..." << endl; }
//       // locateM().findCell(NULL, positions_, cellIDs_, numberOfParticles());
//       findCells();
//       if (verbose_) { Info << "cfdemCloudMix::findCell() - done\n" << endl; }

//       // 计算颗粒尺寸与其周围网格平均尺寸的比值, 并将颗粒索引按照颗粒尺寸归类
//       if (verbose_) { Info << "voidFractionModel::getDimensionRatios()..." << endl; }
//       const_cast<voidFractionModel&>(voidFractionM()).getDimensionRatios(dimensionRatios_,
//                                                                          fineParticleIndexs_,
//                                                                          middleParticleIndexs_,
//                                                                          coarseParticleIndexs_,
//                                                                          missingParticleIndexs_);
//       if (verbose_) { Info << "voidFractionModel::getDimensionRatios() - done\n " << endl; }

//       // 设置小颗粒空隙率场
//       if (verbose_) Info << "cfdemCloudMix::setvoidFraction()..." << endl;
//       voidFractionM().voidFractionModelInit(voidfractions_,
//                                             particleWeights_,
//                                             particleVolumes_,
//                                             particleV_,
//                                             dimensionRatios_);
//       // @brief 设置空隙率和流体体积分数前, 必须先调用 voidFractionModelInit 函数
//       voidFractionM().setVoidFractionAndVolumeFraction(NULL,
//                                                        volumefractions_,
//                                                        voidfractions_,
//                                                        particleWeights_,
//                                                        particleVolumes_,
//                                                        particleV_,
//                                                        dimensionRatios_);
//       // setVoidFraction();
//       if (verbose_) Info << "cfdemCloudMix::setvoidFraction() - done\n" << endl;

//       // 设置局部平均颗粒速度场，将小颗粒的速度 velocities_ 映射到局部平局小颗粒速度，本质上是将拉格朗日场映射到欧拉场
//       // setVectorAverages();
//       if (verbose_) Info << "cfdemCloudMix::setVectorAverage()..." << endl;
//       averagingM().setVectorAverage(averagingM().UsNext(),
//                                     velocities_,
//                                     particleWeights_,
//                                     averagingM().UsWeightField(),
//                                     NULL,
//                                     NULL,
//                                     false);
//       if (verbose_) Info << "cfdemCloudMix::setVectorAverage() - done\n" << endl;

//       // smoothen "next" fields
//       smoothingM().dSmoothing();
//       smoothingM().smoothen(voidFractionM().voidFractionNext());
//       // only smoothen if we use implicit force coupling in cells void of particles
//       // because we need unsmoothened Us field to detect cells for explicit force coupling
//       if (!treatVoidCellsAsExplicitForce()) {
//         smoothingM().smoothenReferenceField(averagingM().UsNext());
//       }
//       Info << "\n// * * * * * * * * * * " << "Coupling - done" << " * * * * * * * * * * //\n" << endl;
//     }  // End of doCoupleNow

//     // 获取 timeStepFraction，以防止使用自适应流体时间步
//     if (verbose_) { Info << "timeStepFraction() = " << dataExchangeM().timeStepFraction() << endl; }
//     if (dataExchangeM().timeStepFraction() > 1.001 ||
//        dataExchangeM().timeStepFraction() < -0.001) {
//       FatalError << "cfdemCloud::dataExchangeM().timeStepFraction() = "<< dataExchangeM().timeStepFraction()
//         << " must be >1 or <0 : This might be due to the fact that you used a adjustable CFD time step. Please use a fixed CFD time step."
//         << abort(FatalError);
//     }

//     // update voidFractionField
//     setAlpha(alpha);
//     if (dataExchangeM().couplingStep() < 2) {
//       alpha.oldTime() = alpha; // supress volume src
//       alpha.oldTime().correctBoundaryConditions();
//     }
//     alpha.correctBoundaryConditions();

//     // calc ddt(voidfraction)
//     calcDdtVoidfraction(alpha, Us);

//     // update mean particle velocity Field
//     Us = averagingM().UsInterp();
//     Us.correctBoundaryConditions();

//     if (doCouple) {
//       // set particles forces
//       if (verbose_) { Info << "cfdemCloudMix::setForce()..." << endl; }
//       setForces();
//       if (verbose_) { Info << "cfdemCloudMix::setForce() - done\n" << endl; }

//       if (verbose_) { Info << "cfdemCloudMix::ncalcMultiphaseTurbulence()..." << endl; }
//       calcMultiphaseTurbulence();
//       if (verbose_) { Info << "cfdemCloudMix::calcMultiphaseTurbulence() - done\n" << endl; }

//       // get next force field
//       if (verbose_) { Info << "cfdemCloudMix::nsetParticleForceField()..." << endl; }
//       setParticleForceField();
//       if (verbose_) { Info << "cfdemCloudMix::setParticleForceField() - done\n" << endl; }

//       // write DEM data
//       if (verbose_) { Info << "cfdemCloudMix::ngiveDEMdata()..." << endl; }
//       giveDEMdata();
//       if (verbose_) { Info << "cfdemCloudMix::giveDEMdata() - done\n" << endl; }

//       dataExchangeM().couple(1);
//     }
//     if (verbose_) {
//       // #include "debugInfo.H"
//     }
//     // do particle IO
//     IOM().dumpDEMdata();
//   }  // End of ignore
//   Info << "\n// * * * * * * * * * * * * * * * " << "Foam::cfdemCloudMix::evolve() - done" << " * * * * * * * * * * * * * * * //\n" << endl;
//   return doCouple;
// }

// @brief 更新函数
// @param alpha  <[in, out] 小颗粒空隙率场
// @param Us     <[in, out] 局部平均小颗粒速度场
// @param U      <[in] 流体速度场
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
      forAll (mesh_.C(),cellI) {
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
        // locateM().findCell(NULL, positions_, cellIDs_, numberOfParticles());
        findCells();
        if (verbose_) { Info << "cfdemCloudMix::findCell() - done\n" << endl; }

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
