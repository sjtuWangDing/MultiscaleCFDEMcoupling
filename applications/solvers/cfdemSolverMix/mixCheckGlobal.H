// include this file before the cloud class is created.
{
  IOdictionary dict(IOobject("couplingProperties",
                              mesh.time().constant(),
                              mesh,
                              IOobject::MUST_READ,
                              IOobject::NO_WRITE));

  // check if NO explicitCouple model is declared - otherwise error!
  wordList mc(dict.lookup("momCoupleModels"));
  for (int i = 0; i < mc.size(); ++i) {
    if ("explicitCouple" == mc[i] || "mixExplicitCouple" == mc[i]) {
      FatalError << "mixCheckGlobal.H: "
        << "You specified an explicit momemtum coupling model, which is not suitable for this solver. Please use implicitCouple or mixImplicitCouple model"
        << abort(FatalError);
    }
  }

  // check use mix force model - otherwise error!
  wordList fm(dict.lookup("forceModels"));
  for (int i = 0; i < fm.size(); ++i) {
    if ("mixArchimedes" == fm[i] ||
        "mixDiFeliceDrag" == fm[i] ||
        "mixArchimedesIB" == fm[i] ||
        "mixShirgaonkarIB" == fm[i]) { continue; }
    FatalError << "mixCheckGlobal.H: "
      << "You specified a wrong force submodel, which is not suitable for this solver. Please use mixArchimedes or mixDiFeliceDrag or mixArchimedesIB or ShirgaonkarIB model"
      << abort(FatalError);
  }
}