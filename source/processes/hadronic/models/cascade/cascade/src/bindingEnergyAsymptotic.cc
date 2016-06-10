//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: bindingEnergyAsymptotic.cc 66241 2012-12-13 18:34:42Z gunter $
//
// 20100202  M. Kelsey -- Eliminate unnecessary use of std::pow()
// 20100914  M. Kelsey -- Migrate to integer A and Z; discard unused verbose
// 20150622  M. Kelsey -- For integer G4cbrt, move Z^4 outside.

#include "G4InuclSpecialFunctions.hh"

G4double G4InuclSpecialFunctions::bindingEnergyAsymptotic(G4int A, G4int Z) {
  // calculates the nuclei binding energy 
  // using smooth liquid high energy formula  
  G4double X = (1.0 - 2.0*Z/A); X *= X;
  G4double X1 = G4cbrt(A);
  G4double X2 = X1 * X1;
  G4double X3 = 1.0 / X1;
  G4double X4 = 1.0 / X2;
  G4double X5 = (1.0 - 0.62025 * X4); X5 *= X5;
  G4double Z13 = G4cbrt(Z);

  G4double DM = 17.035 * (1.0 - 1.846 * X) * A -
    25.8357 * (1.0 - 1.712 * X) * X2 * X5 -
    0.779 * Z * (Z - 1) * X3 * 
    (1.0 - 1.5849 * X4 + 1.2273 / A + 1.5772 * X4 * X4) +
    0.4328 * Z13*Z13*Z13*Z13 * X3 * 
    (1.0 - 0.57811 * X3 - 0.14518 * X4 + 0.496 / A);

  return DM;
}
