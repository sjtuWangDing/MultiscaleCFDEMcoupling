{
  Info << "cfdemSolverMix: mixCheckModelNone..." << endl;
  double err = 1.e-15;
  if (modelType == "none") {
    // check volumefraction == 1
    forAll (volumefraction, cellI) {
      if (fabs(volumefraction[cellI] - 1) > err) {
        Pout << "cellI: " << cellI << " volumefraction: " << volumefraction[cellI] << endl;
      }
    }
  }
  Info << "cfdemSolverMix: mixCheckModelNone - done\n" << endl;
}