// check model type
// referring to: Zhou et al. 2010,JFM

{
  // 根据 Zhou 的论文, 定义了三种模型, setI, setII, setIII
  // "Bfull" 模型对应 setI
  // "A" 模型对应 setII
  // "B" 模型对应 setIII
  auto checkUseForceModel = [&particleCloud] (const std::string& forceModelName) -> bool {
    bool found = false;
    forAll(particleCloud.forceModels(), i) {
      if (particleCloud.forceModels()[i] == forceModelName.c_str()) { found = true; }
    }
    return found;
  };

  word modelType = particleCloud.modelType();
  Info << "\nsolving volume averaged Navier Stokes equations of type " << modelType << endl << endl;

  if (modelType == "Bfull" || modelType == "A") {
    // 在 "Bfull" 以及 "A" 模型中, 流体对颗粒的作用力包含:
    // 阻力 + 压力梯度力 + 粘性应力 + 非主导力(虚拟质量力 + 升力(Saffman force + Magnus force) + Basset force)
    // 但是在 "A" 模型中, 颗粒对流体的作用力为包含: 阻力 + 非主导力, 而在 "Bfull" 模型中, 颗粒对流体对作用力包含上面全部力
    // 所以这里我们检查必须使用压力梯度力, 粘性力 和 不能使用 Archimedes force(压力梯度力中包含力浮力)

    // 不能使用 Archimedes force
    if (checkUseForceModel("Archimedes") || checkUseForceModel("mixArchimedes")) {
      FatalError << "mixCheckModelType.H: can not use Archimedes or mixArchimedes with model type " << modelType << abort(FatalError);
    }

    // 必须使用压力梯度力
    if (!checkUseForceModel("gradPForce")) {
      FatalError << "mixCheckModelType.H: gradPForce not found with model type " << modelType << abort(FatalError);
    }

    // 必须使用粘性力
    if (!checkUseForceModel("viscForce")) {
      FatalError << "mixCheckModelType.H: viscForce not found with model type " << modelType << abort(FatalError);
    }
  } else if (modelType == "B") {
    // 必须使用 Archimedes force, 且只能指定 Archimedes or mixArchimedes 中任意一个
    if (checkUseForceModel("Archimedes") && checkUseForceModel("mixArchimedes")) {
      FatalError << "mixCheckModelType.H: Archimedes or mixArchimedes model both found! You must use only one of them\n" << abort(FatalError);
    } else if ((!checkUseForceModel("Archimedes")) && (!checkUseForceModel("mixArchimedes"))) {
      FatalError << "mixCheckModelType.H: Archimedes and mixArchimedes model both not found! You must use only one of them\n" << abort(FatalError);
    }

    // 不能使用压力梯度力
    if (checkUseForceModel("gradPForce")) {
      FatalError << "mixCheckModelType.H: can not use gradPForce with model type B!\n" << abort(FatalError);
    }

    // 不能使用粘性力
    if (checkUseForceModel("viscForce")) {
      FatalError << "mixCheckModelType.H: can not use viscForce with model type B!\n" << abort(FatalError);
    }
  } else if (modelType == "none") {
    Info << "\nsolving volume averaged Navier Stokes equations of type none\n"<< endl;
  } else {
    FatalError << "\nno suitable model type specified:" << modelType << endl << abort(FatalError);
  }
}
