// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4AlphaEvaporationProbability.hh"

G4AlphaEvaporationProbability::G4AlphaEvaporationProbability() :
		G4EvaporationProbability(4,2,4) // A,Z,Gamma
{
	const G4int NumExcitedStates = 31+1;
	ExcitEnergies.reshape(NumExcitedStates);
	ExcitSpins.reshape(NumExcitedStates);
	for (G4int i = 0; i < NumExcitedStates; i++) {
		ExcitEnergies(i) = 0.0;
		ExcitSpins(i) = 0;
	}

	ExcitEnergies(18) = 7.98*MeV;
	ExcitEnergies(20) = 6.90*MeV;
	ExcitEnergies(25) = 5.83*MeV;
	ExcitEnergies(26) = 8.57*MeV;
	ExcitEnergies(31) = 5.33*MeV;

	ExcitSpins(18) = 4;
	ExcitSpins(20) = 6;
	ExcitSpins(25) = 7;
	ExcitSpins(26) = 4;
	ExcitSpins(31) = 13;
	
	SetExcitationEnergiesPtr(&ExcitEnergies);
	SetExcitationSpinsPtr(&ExcitSpins);		
}


G4AlphaEvaporationProbability::G4AlphaEvaporationProbability(const G4AlphaEvaporationProbability &right)
{
 G4Exception("G4AlphaEvaporationProbability::copy_constructor meant to not be accessable");
}




const G4AlphaEvaporationProbability & G4AlphaEvaporationProbability::
operator=(const G4AlphaEvaporationProbability &right)
{
  G4Exception("G4AlphaEvaporationProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4AlphaEvaporationProbability::operator==(const G4AlphaEvaporationProbability &right) const
{
  return false;
}

G4bool G4AlphaEvaporationProbability::operator!=(const G4AlphaEvaporationProbability &right) const
{
  return true;
}

G4double G4AlphaEvaporationProbability::CCoeficient(const G4double aZ) const
{
	// Data comes from 
	// Dostrovsky, Fraenkel and Friedlander
	// Physical Review, vol 116, num. 3 1959
	// 
	// const G4int size = 5;
	// G4double Zlist[5] = { 10.0, 20.0, 30.0, 50.0, 70.0};
	//	G4double Calpha[5] = { 0.10, 0.10, 0.10, 0.08, 0.06};
	G4double C = 0.0;
	
	
	if (aZ <= 30) {
		C = 0.10;
	} else if (aZ <= 50) {
		C = 0.1 + -((aZ-50.)/20.)*0.02;
	} else if (aZ < 70) {
		C = 0.08 + -((aZ-70.)/20.)*0.02;
	} else {
		C = 0.06;
	}
	return C;
}
