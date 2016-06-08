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


#include "G4ProtonEvaporationProbability.hh"

G4ProtonEvaporationProbability::G4ProtonEvaporationProbability() :
		G4EvaporationProbability(1,1,2) // A,Z,Gamma
{
	const G4int NumExcitedStates = 31+1;
	ExcitEnergies.reshape(NumExcitedStates);
	ExcitSpins.reshape(NumExcitedStates);
	for (G4int i = 0; i < NumExcitedStates; i++) {
		ExcitEnergies(i) = 0.0;
		ExcitSpins(i) = 0;
	}
	ExcitEnergies(15) = 5.96*MeV;
	ExcitEnergies(17) = 1.74*MeV;
	ExcitEnergies(18) = 4.44*MeV;
	ExcitEnergies(19) = 1.67*MeV;
	ExcitEnergies(20) = 4.32*MeV;
	ExcitEnergies(22) = 3.68*MeV;
	ExcitEnergies(23) = 6.69*MeV;
	ExcitEnergies(25) = 3.95*MeV;
	ExcitEnergies(26) = 6.32*MeV;
	ExcitEnergies(27) = 0.30*MeV;
	ExcitEnergies(28) = 6.18*MeV;
	ExcitEnergies(29) = 6.92*MeV;
	ExcitEnergies(30) = 3.06*MeV;
	ExcitEnergies(31) = 3.57*MeV;


	ExcitSpins(15) = 8;
	ExcitSpins(17) = 1;
	ExcitSpins(18) = 6;
	ExcitSpins(19) = 5;
	ExcitSpins(20) = 6;
	ExcitSpins(22) = 4;
	ExcitSpins(23) = 8;
	ExcitSpins(25) = 3;
	ExcitSpins(26) = 4;
	ExcitSpins(27) = 7;
	ExcitSpins(28) = 4;
	ExcitSpins(29) = 5;
	ExcitSpins(30) = 2;
	ExcitSpins(31) = 10;
	
	SetExcitationEnergiesPtr(&ExcitEnergies);
	SetExcitationSpinsPtr(&ExcitSpins);	
	
}

G4ProtonEvaporationProbability::G4ProtonEvaporationProbability(const G4ProtonEvaporationProbability &right)
{
 G4Exception("G4ProtonEvaporationProbability::copy_constructor meant to not be accessable");
}

const G4ProtonEvaporationProbability & G4ProtonEvaporationProbability::
operator=(const G4ProtonEvaporationProbability &right)
{
  G4Exception("G4ProtonEvaporationProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4ProtonEvaporationProbability::operator==(const G4ProtonEvaporationProbability &right) const
{
  return false;
}

G4bool G4ProtonEvaporationProbability::operator!=(const G4ProtonEvaporationProbability &right) const
{
  return true;
}

G4double G4ProtonEvaporationProbability::CCoeficient(const G4double aZ) const
{
	// Data comes from 
	// Dostrovsky, Fraenkel and Friedlander
	// Physical Review, vol 116, num. 3 1959
	// 
	// const G4int size = 5;
	// G4double Zlist[5] = { 10.0, 20.0, 30.0, 50.0, 70.0};
	// G4double Cp[5] = { 0.50, 0.28, 0.20, 0.15, 0.10};
	G4double C = 0.0;
	
	if (aZ >= 70) {
		C = 0.10;
	} else {
		C = ((((0.15417e-06*aZ) - 0.29875e-04)*aZ + 0.21071e-02)*aZ - 0.66612e-01)*aZ + 0.98375;
	}
	
	return C;
	
}
