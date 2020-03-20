// check if NO explicitCouple model is declared - otherwise error!
// include this file before the cloud class is created.
{
  IOdictionary dict(IOobject("couplingProperties",
                              mesh.time().constant(),
                              mesh,
                              IOobject::MUST_READ,
                              IOobject::NO_WRITE));

  wordList mc(dict.lookup("momCoupleModels"));

  for (int i = 0; i < mc.size(); ++i) {

    if ("explicitCouple" == mc[i] || "mixExplicitCouple" == mc[i]) {
      FatalError << "mixCheckCouple.H: You specified an explicit momemtum coupling model, which is not suitable for this solver. Please use implicitCouple or mixImplicitCouple model"
      << abort(FatalError);
    }
  }
}
