//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
							  G4double Z) {
  G4int verboseLevel = 2;

  if (verboseLevel > 3) {
    G4cout << " >>> G4InuclSpecialFunctions::bindingEnergyAsymptotic" << G4endl;
  }

  // calculates the nuclei binding energy 
  // using smooth liquid high energy formula  
  G4double X = (1.0 - 2.0 * Z / A) * (1.0 - 2.0 * Z / A);
  G4double X1 = pow(A, 0.3333333);
  G4double X2 = X1 * X1;
  G4double X3 = 1.0 / X1;
  G4double X4 = 1.0 / X2;
  G4double DM = 17.035 * (1.0 - 1.846 * X) * A -
    25.8357 * (1.0 - 1.712 * X) * X2 * 
    ((1.0  -0.62025 * X4) * (1.0 - 0.62025 * X4)) -
    0.779 * Z * (Z - 1.0) * X3 * 
    (1.0 - 1.5849 * X4 + 1.2273 / A + 1.5772 * X4 * X4) +
    0.4328 * pow(Z, 1.333333) * X3 * 
    (1.0 - 0.57811 * X3 - 0.14518 * X4 + 0.496 / A);

  return DM;
}
