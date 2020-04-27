{
  Info << "cfdemSolverMix: mixCheckModelNone..." << endl;
  double err = 1.e-15;
  if (modelType == "none") {
    for (int cellI = 0; cellI < voidfraction.size(); ++cellI) {
      // (1) check voidfraction == 1
      if (fabs(voidfraction[cellI] - 1) > err) {
        // Pout << "cellI: " << cellI << " voidfraction: " << voidfraction[cellI] << endl;
        FatalError << "cellI: " << cellI << " voidfraction: " << voidfraction[cellI] << abort(FatalError);
      }

      // (2) check Ksl == 0
      if (fabs(Ksl[cellI] - 0) > err) {
        // Pout << "cellI: " << cellI << " Ksl: " << Ksl[cellI] << endl;
        FatalError << "cellI: " << cellI << " Ksl: " << Ksl[cellI] << abort(FatalError);
      }

      // (3) check Us[cellI] == vector::zero
      if (Us[cellI] != vector::zero) {
        // Pout << "cellI: " << cellI << " Us: " << Us[cellI] << endl;
        FatalError << "cellI: " << cellI << " Us: " << Us[cellI] << abort(FatalError);
      }

      // (4) check phi == phiByVoidfraction
      if (fabs(phi[cellI] - phiByVoidfraction[cellI]) > err) {
        // Pout << "cellI: " << cellI << " phi: " << phi[cellI] << " != " << phiByVoidfraction[cellI] << endl;
        FatalError << "cellI: " << cellI << " phi: " << phi[cellI] << " != " << phiByVoidfraction[cellI] << abort(FatalError);
      }
    }
  }
  Info << "cfdemSolverMix: mixCheckModelNone - done\n" << endl;
}