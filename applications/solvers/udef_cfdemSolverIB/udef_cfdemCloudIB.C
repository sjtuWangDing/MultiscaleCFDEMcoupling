#include "fileName.H"
#include "udef_cfdemCloudIB.H"
#include "voidFractionModel.H"
#include "forceModel.H"
#include "locateModel.H"
#include "dataExchangeModel.H"
#include "IOModel.H"
#include "mpi.h"
#include "IOmanip.H"
#include "OFversion.H"

namespace Foam {

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

udef_cfdemCloudIB::udef_cfdemCloudIB(const fvMesh& mesh):
  cfdemCloud(mesh),
  angularVelocities_(NULL),
  DEMTorques_(NULL),
  pRefCell_(readLabel(mesh.solutionDict().subDict("PISO").lookup("pRefCell"))),
  pRefValue_(readScalar(mesh.solutionDict().subDict("PISO").lookup("pRefValue"))),
  haveEvolvedOnce_(false),
  skipLagrangeToEulerMapping_(false),
  skipAfter_(false),
  timeStepsToSkip_(0),
  calculateTortuosity_(false),
  frontMeshRefine_(false)
{
  if (this->couplingProperties().found("skipLagrangeToEulerMapping")) {
    Info << "Will skip lagrange-to-Euler mapping..." << endl;
    skipLagrangeToEulerMapping_ = true;
  }

  if (this->couplingProperties().found("timeStepsBeforeSkipping")) {
    skipAfter_ = true;
    timeStepsToSkip_ = readScalar(this->couplingProperties().lookup("timeStepsBeforeSkipping"));
    Info << "Will skip LagrangeToEuler mapping after " << timeStepsToSkip_ << " time steps" <<  endl;
  }

  if (this->couplingProperties().found("tortuosity")) {
    calculateTortuosity_ = true;
    flowDir_ = this->couplingProperties().subDict("tortuosity").lookup("flowDirection");
    flowDir_ = flowDir_ / mag(flowDir_);
    Info << "Will calculate tortuosity in the mean flow direction ("
      << flowDir_[0] << ", " << flowDir_[1] << ", " << flowDir_[2] << ")" << endl;
  }

  // Must check for walls in case of checkPeriodicCells
  // periodic check will mirror particles and probing points to ensure proper behavior near processor bounds
  if (checkPeriodicCells_) {
    // Enforce reading of the blocking for periodic checks
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
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

udef_cfdemCloudIB::~udef_cfdemCloudIB()
{
    dataExchangeM().destroy(angularVelocities_, 3);
    dataExchangeM().destroy(dragPrev_, 3);
    dataExchangeM().destroy(DEMTorques_, 3);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool udef_cfdemCloudIB::evolve(volScalarField& alpha,
                               volScalarField& interFace) {
  numberOfParticlesChanged_ = false;
  arraysReallocated_ = false;
  bool doCouple = false;

  if (skipAfter_ && timeStepsToSkip_ < 1) {
    skipLagrangeToEulerMapping_=true;
  }

  if (!writeTimePassed_ && mesh_.time().outputTime()) {
    writeTimePassed_ = true;
  }

  if (dataExchangeM().doCoupleNow()) {
    Info << "\n timeStepFraction() = " << dataExchangeM().timeStepFraction() << endl;
    dataExchangeM().couple(0);
    doCouple = true;
    // Info << "skipLagrangeToEulerMapping_: " << skipLagrangeToEulerMapping_ 
    //   << " haveEvolvedOnce_: " << haveEvolvedOnce_ << endl;
    if (!skipLagrangeToEulerMapping_ || !haveEvolvedOnce_) {
      // get DEM data
      if (verbose_) { Info << "- getDEMdata()" << endl; }
      getDEMdata();
      Info << "nr particles = " << numberOfParticles() << endl;
      if (verbose_) { Info << "getDEMdata done" << endl; }

      // search cellID of particles
      if (verbose_) { Info << "- findCell()" << endl; }
      locateM().findCell(NULL, positions_, cellIDs_, numberOfParticles());
      if (verbose_) { Info << "findCell done" << endl; }

      // set void fraction field
      if (verbose_) { Info << "- setvoidFraction()" << endl; }
      voidFractionM().setvoidFraction(NULL, voidfractions_, particleWeights_,
                                      particleVolumes_, particleV_);
      if (verbose_) { Info << "setvoidFraction done" << endl; }

      // set interface
      if (verbose_) { Info << "setInterFace" << endl; }
      setInterFace(interFace);
      if (verbose_) { Info << "setInterFace done" << endl; }
    }

    // update voidFractionField
    // there might be a better approach, see cfdemCloud.C
    alpha == voidFractionM().voidFractionNext();
    alpha.correctBoundaryConditions();

    // set particles forces
    if (verbose_) { Info << "- setForce(forces_)" << endl; }
    for(int index = 0; index < numberOfParticles_; ++index){
      for(int i = 0; i < 3; i++){
        impForces_[index][i] = 0;
        expForces_[index][i] = 0;
        DEMForces_[index][i] = 0;
      }
    }
    for (int i = 0; i < nrForceModels(); i++) {
      forceM(i).setForce();
    }
    if (verbose_) { Info << "setForce done." << endl; }

    // write DEM data
    if (verbose_) { Info << " -giveDEMdata()" << endl; }
    giveDEMdata();

    dataExchangeM().couple(1);
    haveEvolvedOnce_=true;
  }
  Info << "evolve done." << endl;
  // do particle IO
  IOM().dumpDEMdata();
  if (skipAfter_) {
    timeStepsToSkip_--;
    Info << "Will skip LagrangeToEuler mapping after " << timeStepsToSkip_ << " time steps" <<  endl;
  }
  return doCouple;
}

void udef_cfdemCloudIB::getDEMdata() {
  cfdemCloud::getDEMdata();
  dataExchangeM().getData("omega", "vector-atom", angularVelocities_);
}

bool udef_cfdemCloudIB::reAllocArrays() const {
  if (cfdemCloud::reAllocArrays()) {
    Info << "Foam::cfdemCloudIB::reAllocArrays()" << endl;
    dataExchangeM().allocateArray(angularVelocities_, 0, 3);
    dataExchangeM().allocateArray(dragPrev_, 0, 3);
    dataExchangeM().allocateArray(DEMTorques_, 0, 3);
    return true;
  }
  return false;
}

void udef_cfdemCloudIB::giveDEMdata() {
  cfdemCloud::giveDEMdata();
  dataExchangeM().giveData("hdtorque", "vector-atom", DEMTorques_);
}

void udef_cfdemCloudIB::setParticleVelocity(volVectorField& U) {
  label cellI = 0;
  vector uParticle(0, 0, 0);
  vector rVec(0, 0, 0);
  vector velRot(0, 0, 0);
  vector angVel(0, 0, 0);
  for(int index = 0; index < numberOfParticles(); index++) {
    for(int subCell = 0; subCell < cellsPerParticle()[index][0]; subCell++) {
      // 获取第 index 个颗粒覆盖第第 subCell 个网格编号
      cellI = cellIDs()[index][subCell];
      if (cellI >= 0) {
        // 计算颗粒中心到网格中心的相对矢量
        for (int i = 0; i < 3; ++i) {
          rVec[i] = U.mesh().C()[cellI][i] - position(index)[i];
        }
        // 获取颗粒角速度
        for (int i = 0; i < 3; ++i) {
          angVel[i] = angularVelocities()[index][i];
        }
        // 计算网格中心处的转动速度
        velRot = angVel ^ rVec;
        // 计算网格中心的颗粒速度 = 平动速度 + 转动速度
        for (int i = 0; i < 3; ++i) {
          uParticle[i] = velocities()[index][i] + velRot[i];
        }
        // 计算网格中的混合流体速度 = (1 - 空隙率) * 颗粒速度 + 空隙率 * 流体速度
        U[cellI] = (1 - voidfractions_[index][subCell]) * uParticle + voidfractions_[index][subCell] * U[cellI];
      }
    }
  }
  U.correctBoundaryConditions();
}

// 因为在速度映射过程中，在混合流体中（颗粒内部和外部），速度是连续的，但是在颗粒边界处，会出现颗粒不连续的情况，所以这里引入 phiIB 这个标量去修正速度场。
// 设：在求解完压力泊松方程后，速度场为 U1，下个时间步的速度场为 U2（即待求速度场），那么 U2 = U1 - div(phiIB)，同时 U2 应该满足连续方程，即 ddt(voidfraction) + div(voidfraction * U2) = 0，则得到关于 phiIB 的修正方程：
// ddt(voidfraction) + div(voidfraction * U1) = div(voidfraction, div(phiIB))
void udef_cfdemCloudIB::calcVelocityCorrection(volScalarField& p,
                                               volVectorField& U,
                                               volScalarField& phiIB,
                                               volScalarField& voidfraction) {
  setParticleVelocity(U);
  // make field divergence free - set reference value in case it is needed
  fvScalarMatrix phiIBEqn
  (
    fvm::laplacian(phiIB) == fvc::div(U) + fvc::ddt(voidfraction)
  );
  if (phiIB.needReference()) {
      phiIBEqn.setReference(pRefCell_, pRefValue_);
  }
  phiIBEqn.solve();

  U = U - fvc::grad(phiIB);
  U.correctBoundaryConditions();

  // correct the pressure as well
  p = p + phiIB / U.mesh().time().deltaT();
  p.correctBoundaryConditions();

  if (couplingProperties_.found("checkinterface")) {
    Info << "checking no-slip on interface..." << endl;
    // #include "checkInterfaceVelocity.H" // TODO: check carefully!
  }
}


// defines the mesh refinement zone around a particle twice the particle size in each direction
void udef_cfdemCloudIB::setInterFace(volScalarField& interFace) {
  interFace == dimensionedScalar("zero", interFace.dimensions(), 0.);
  for(int par = 0; par < numberOfParticles(); par++) {
    // 获取颗粒中心位置
    vector ParPos(positions()[par][0], positions()[par][1], positions()[par][2]);
    const boundBox& globalBb = mesh().bounds();
    double skin = 2.0;
    forAll(mesh_.C(), cellI) {
      // 获取网格中心位置
      vector posC = mesh_.C()[cellI];
      if (checkPeriodicCells_) {
        // Some cells may be located on the other side of a periodic boundary.
        // In this case, the particle center has to be mirrored in order to correctly
        // evaluate the interpolation points.
        vector minPeriodicParticlePos = ParPos;
        voidFractionM().minPeriodicDistance(par, posC, ParPos, globalBb,
                                            minPeriodicParticlePos, 
                                            wall_periodicityCheckRange());
        ParPos = minPeriodicParticlePos;
      }
      // 如果网格在颗粒的两倍直径内，则设置该网格的 interFace
      double value = voidFractionM().pointInParticle(par, ParPos, posC, skin);
      if(value <= 0.0) {
        interFace[cellI] = value + 1.0;
      }
    }
  }
}

}  // End of namespace Foam