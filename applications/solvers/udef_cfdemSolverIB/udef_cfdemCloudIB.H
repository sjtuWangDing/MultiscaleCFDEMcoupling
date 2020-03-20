#ifndef udef_cfdemCloudIB_H
#define udef_cfdemCloudIB_H

#include "cfdemCloud.H"

namespace Foam {

class udef_cfdemCloudIB: public cfdemCloud {

protected:
  // @brief 颗粒角速度
  mutable double** angularVelocities_;

  // @brief 流体对颗粒的力矩
  mutable double** DEMTorques_;

  label pRefCell_;
  scalar pRefValue_;
  mutable bool haveEvolvedOnce_;

  // @brief 是否跳过拉格朗日映射到欧拉场
  mutable bool skipLagrangeToEulerMapping_;

  // @brief 是否在指定时间步后 Skip Lagrange To Euler Mapping.
  mutable bool skipAfter_;

  // @brief 指定在多少时间步后 Skip Lagrange To Euler Mapping.
  mutable int  timeStepsToSkip_;

  mutable bool calculateTortuosity_;

  // @brief Mean flow direction
  vector flowDir_;
  volScalarField* frontMeshRefineField_;
  bool frontMeshRefine_;

public:
  // Constructors
  udef_cfdemCloudIB(const fvMesh& mesh);

  // Destructor
  ~udef_cfdemCloudIB();

  // Member Functions
  void giveDEMdata();

  void getDEMdata();

  bool reAllocArrays() const;

  bool evolve(volScalarField&, volScalarField&);

  void calcVelocityCorrection(volScalarField&, volVectorField&, volScalarField&, volScalarField&);

  virtual void setParticleVelocity(volVectorField&);

  void setInterFace(volScalarField&);

  double getTortuosity(vector dir);

  // @brief Function to allow for setting and activating special refinement operations
  // @param refine_ - pointer is meant for a field that is used refinement outside of the cloud, can point to the same field as interFace.
  void setRefinementField(volScalarField* refine_) {
    frontMeshRefineField_ = refine_;
    frontMeshRefine_ = true;
  }

  // Access Functions.
  double ** angularVelocities() const {
    return angularVelocities_;
  };

  double** DEMTorques() const {
    return DEMTorques_;
  }

  vector angularVelocity(int index) {
    vector vel;
    for (int i = 0; i < 3; ++i) {
      vel[i] = angularVelocities_[index][i];
    }
    return vel;
  }
};  // End of class udef_cfdemCloudIB

}  // End of namespace Foam

#endif  // udef_cfdemCloudIB_H
