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


#include "G4TritonEvaporationProbability.hh"

G4TritonEvaporationProbability::G4TritonEvaporationProbability() :
		G4EvaporationProbability(3,1,6) // A,Z,Gamma
{
	const G4int NumExcitedStates = 31+1;
	ExcitEnergies.reshape(NumExcitedStates);
	ExcitSpins.reshape(NumExcitedStates);
	for (G4int i = 0; i < NumExcitedStates; i++) {
		ExcitEnergies(i) = 0.0;
		ExcitSpins(i) = 0;
	}
	
	ExcitEnergies(15) = 6.26*MeV;
	ExcitEnergies(17) = 3.59*MeV;
	ExcitEnergies(18) = 6.76*MeV;
	ExcitEnergies(20) = 6.34*MeV;
	ExcitEnergies(23) = 7.34*MeV;
	ExcitEnergies(25) = 5.11*MeV;
	ExcitEnergies(26) = 7.57*MeV;
	ExcitEnergies(28) = 7.28*MeV;
	ExcitEnergies(31) = 4.46*MeV;

	ExcitSpins(15) = 5;
	ExcitSpins(17) = 5;
	ExcitSpins(18) = 10;
	ExcitSpins(20) = 2;
	ExcitSpins(23) = 5;
	ExcitSpins(25) = 5;
	ExcitSpins(26) = 8;
	ExcitSpins(28) = 8;
	ExcitSpins(31) = 3;

	SetExcitationEnergiesPtr(&ExcitEnergies);
	SetExcitationSpinsPtr(&ExcitSpins);	
}

G4TritonEvaporationProbability::G4TritonEvaporationProbability(const G4TritonEvaporationProbability &right)
{
 G4Exception("G4TritonEvaporationProbability::copy_constructor meant to not be accessable");
}




const G4TritonEvaporationProbability & G4TritonEvaporationProbability::
operator=(const G4TritonEvaporationProbability &right)
{
  G4Exception("G4TritonEvaporationProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4TritonEvaporationProbability::operator==(const G4TritonEvaporationProbability &right) const
{
  return false;
}

G4bool G4TritonEvaporationProbability::operator!=(const G4TritonEvaporationProbability &right) const
{
  return true;
}

G4double G4TritonEvaporationProbability::CCoeficient(const G4double aZ) const
{
	// Data comes from 
	// Dostrovsky, Fraenkel and Friedlander
	// Physical Review, vol 116, num. 3 1959
	// 
	// const G4int size = 5;
	// G4double Zlist[5] = { 10.0, 20.0, 30.0, 50.0, 70.0};
	// G4double Cp[5] = { 0.50, 0.28, 0.20, 0.15, 0.10};
	// C for triton is equal to C for protons divided by 3
	G4double C = 0.0;
	
	if (aZ >= 70) {
		C = 0.10;
	} else {
		C = ((((0.15417e-06*aZ) - 0.29875e-04)*aZ + 0.21071e-02)*aZ - 0.66612e-01)*aZ + 0.98375;
	}
	
	return C/3.0;
	
}
