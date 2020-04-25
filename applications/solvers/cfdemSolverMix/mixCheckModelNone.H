{
  Info << "cfdemSolverMix: mixCheckModelNone..." << endl;
  double err = 1.e-15;
  if (modelType == "none") {
    for (int cellI = 0; cellI < voidfraction.size(); ++cellI) {
      // (1) check voidfraction == 1
      if (fabs(voidfraction[cellI] - 1) > err) {
        Pout << "cellI: " << cellI << " voidfraction: " << voidfraction[cellI] << endl;
      }

      // (2) check Ksl == 0
      if (fabs(Ksl[cellI] - 0) > err) {
        Pout << "cellI: " << cellI << " Ksl: " << Ksl[cellI] << endl;
      }

      // (3) check Us[cellI] == vector::zero
      if (Us[cellI] != vector::zero) {
        Pout << "cellI: " << cellI << " Us: " << Us[cellI] << endl;
      }
    }
  }
  Info << "cfdemSolverMix: mixCheckModelNone - done\n" << endl;
}