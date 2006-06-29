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
#include "G4InuclSpecialFunctions.hh"

G4double G4InuclSpecialFunctions::bindingEnergy(G4double A, G4double Z) {
  G4int verboseLevel = 2;

  if (verboseLevel > 3) {
    G4cout << " >>> G4InuclSpecialFunctions::bindingEnergy" << G4endl;
  }

  // calculates the nuclei binding energy using Kummel or exact or asymptotic
  // high temperature 
  G4double DM;
  G4double AN = A - Z;

  if (AN < 0.1 || Z < 0.1) {
    DM = 0.0;

  } else {

    if (A <= 256.0) {

      if (AN >= 20. && Z >= 20) { 

	if (Z < 1.7 * AN && Z > 0.3 * AN) { // standard
	  DM = bindingEnergyKummel(A, Z);

	} else { // bad case
	  DM = bindingEnergyAsymptotic(A, Z);
	}; 

      } else {

	if (A > 60.0 || Z > 21) { // bad case
	  DM = bindingEnergyAsymptotic(A, Z);

	} else { // exact case
	  DM = bindingEnergyExact(A, Z);
	}; 
      }; 

    } else {
      DM = bindingEnergyAsymptotic(A, Z);
    }; 
  };  

  // G4cout << " A " << A << " Z " << Z << " DM " << DM << G4endl;

  return DM;
}
