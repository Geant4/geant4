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

  if(AN < 0.1 || Z < 0.1) {
    DM = 0.0;
  } else {
    if(A <= 256.0) {
      if(AN >= 20. && Z >= 20) { 
	if(Z < 1.7 * AN && Z > 0.3 * AN) { // standard
	  DM = bindingEnergyKummel(A, Z);
	} else { // bad case
	  DM = bindingEnergyAsymptotic(A, Z);
	}; 
      } else {
	if(A > 60.0 || Z > 21) { // bad case
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




