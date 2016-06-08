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


#include "G4He3EvaporationProbability.hh"

G4He3EvaporationProbability::G4He3EvaporationProbability() :
		G4EvaporationProbability(3,2,6) // A,Z,Gamma
{
	const G4int NumExcitedStates = 31+1;
	ExcitEnergies.reshape(NumExcitedStates);
	ExcitSpins.reshape(NumExcitedStates);
	for (G4int i = 0; i < NumExcitedStates; i++) {
		ExcitEnergies(i) = 0.0;
		ExcitSpins(i) = 0;
	}
	
	ExcitEnergies(18) = 7.29*MeV;
	ExcitEnergies(20) = 6.48*MeV;
	ExcitEnergies(25) = 5.69*MeV;
	ExcitEnergies(26) = 8.31*MeV;
	ExcitEnergies(31) = 5.10*MeV;

	ExcitSpins(18) = 6;
	ExcitSpins(20) = 8;
	ExcitSpins(25) = 3;
	ExcitSpins(26) = 2;
	ExcitSpins(31) = 7;	
	
	SetExcitationEnergiesPtr(&ExcitEnergies);
	SetExcitationSpinsPtr(&ExcitSpins);	
}

G4He3EvaporationProbability::G4He3EvaporationProbability(const G4He3EvaporationProbability &right)
{
 G4Exception("G4He3EvaporationProbability::copy_constructor meant to not be accessable");
}




const G4He3EvaporationProbability & G4He3EvaporationProbability::
operator=(const G4He3EvaporationProbability &right)
{
  G4Exception("G4He3EvaporationProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4He3EvaporationProbability::operator==(const G4He3EvaporationProbability &right) const
{
  return false;
}

G4bool G4He3EvaporationProbability::operator!=(const G4He3EvaporationProbability &right) const
{
  return true;
}

G4double G4He3EvaporationProbability::CCoeficient(const G4double aZ) const
{
	// Data comes from 
	// Dostrovsky, Fraenkel and Friedlander
	// Physical Review, vol 116, num. 3 1959
	// 
	// const G4int size = 5;
	// G4double Zlist[5] = { 10.0, 20.0, 30.0, 50.0, 70.0};
	//	G4double Calpha[5] = { 0.10, 0.10, 0.10, 0.08, 0.06};
	// C for He3 is equal to C for alpha times 4/3
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
	return C*(4.0/3.0);
}
