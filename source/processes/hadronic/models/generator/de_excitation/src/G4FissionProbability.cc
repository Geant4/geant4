// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//


#include "G4FissionProbability.hh"




G4FissionProbability::G4FissionProbability(const G4FissionProbability &right)
{
 G4Exception("G4FissionProbability::copy_constructor meant to not be accessable");
}




const G4FissionProbability & G4FissionProbability::operator=(const G4FissionProbability &right)
{
  G4Exception("G4FissionProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4FissionProbability::operator==(const G4FissionProbability &right) const
{
  return false;
}

G4bool G4FissionProbability::operator!=(const G4FissionProbability &right) const
{
  return true;
}


G4double G4FissionProbability::EmissionProbability(const G4Fragment & fragment, const G4double MaximalKineticEnergy)
	// Compute integrated probability of fission channel
{
	if (MaximalKineticEnergy <= 0.0) return 0.0;
	G4double A = fragment.GetA();
	G4double Z = fragment.GetZ();
	G4double U = fragment.GetExcitationEnergy();
  
	G4double Ucompound = U - EvaporationPairingCorrection(G4int(A),G4int(Z));
	G4double Ufission = U - FissionPairingCorrection(G4int(A),G4int(Z));
  
	G4double SystemEntropy = 2.0*sqrt(theEvapLDP.LevelDensityParameter(A,Z,Ucompound)*Ucompound);
	
	G4double afission = theFissLDP.LevelDensityParameter(A,Z,Ufission);

	G4double Cf = 2.0*sqrt(afission*MaximalKineticEnergy);

	G4double Q1 = 1.0 + (Cf - 1.0)*exp(Cf);
	G4double Q2 = 4.0*pi*afission*exp(SystemEntropy);
	
	G4double probability = Q1/Q2;
 
	return probability;
}

G4double G4FissionProbability::EvaporationPairingCorrection(const G4int A, const G4int Z) const
{
	const G4double PairingConstant = 12.0*MeV;
	const G4int N = A - Z;
	G4double Pair = (1.0 - G4double(Z) + 2.0*(Z/2)) + (1.0 - G4double(N) + 2.0*(N/2));
	G4double PCorrection = Pair*PairingConstant/sqrt(G4double(A));
	return PCorrection;
}

G4double G4FissionProbability::FissionPairingCorrection(const G4int A, const G4int Z) const
{
	const G4double PairingConstant = 14.0*MeV;
	const G4int N = A - Z;
	G4double Pair = (1.0 - G4double(Z) + 2.0*(Z/2)) + (1.0 - G4double(N) + 2.0*(N/2));
	G4double PCorrection = Pair*PairingConstant/sqrt(G4double(A));
	return PCorrection;
}
