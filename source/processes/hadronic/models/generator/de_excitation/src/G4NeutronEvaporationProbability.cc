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


#include "G4NeutronEvaporationProbability.hh"

G4NeutronEvaporationProbability::G4NeutronEvaporationProbability() :
		G4EvaporationProbability(1,0,2) // A,Z,Gamma
{
	const G4int NumExcitedStates = 31+1;
	ExcitEnergies.reshape(NumExcitedStates);
	ExcitSpins.reshape(NumExcitedStates);
	for (G4int i = 0; i < NumExcitedStates; i++) {
		ExcitEnergies(i) = 0.0;
		ExcitSpins(i) = 0;
	}
	
	ExcitEnergies( 9) = 3.56*MeV;
	ExcitEnergies(10) = 0.48*MeV;
	ExcitEnergies(11) = 0.98*MeV;
	ExcitEnergies(12) = 0.43*MeV;
	ExcitEnergies(15) = 3.37*MeV;
	ExcitEnergies(17) = 0.72*MeV;
	ExcitEnergies(18) = 2.13*MeV;
	ExcitEnergies(19) = 0.95*MeV;
	ExcitEnergies(20) = 2.00*MeV;
	ExcitEnergies(21) = 4.44*MeV;
	ExcitEnergies(22) = 3.09*MeV;
	ExcitEnergies(23) = 6.09*MeV;
	ExcitEnergies(25) = 2.31*MeV;
	ExcitEnergies(26) = 5.28*MeV;
	ExcitEnergies(27) = 0.12*MeV;
	ExcitEnergies(28) = 5.22*MeV;
	ExcitEnergies(29) = 6.10*MeV;
	ExcitEnergies(30) = 0.87*MeV;
	ExcitEnergies(31) = 1.98*MeV;


	ExcitSpins( 9) = 1;
	ExcitSpins(10) = 2;
	ExcitSpins(11) = 3;
	ExcitSpins(12) = 2;
	ExcitSpins(15) = 5;
	ExcitSpins(17) = 3;
	ExcitSpins(18) = 2;
	ExcitSpins(19) = 5;
	ExcitSpins(20) = 2;
	ExcitSpins(21) = 5;
	ExcitSpins(22) = 2;
	ExcitSpins(23) = 3;
	ExcitSpins(25) = 1;
	ExcitSpins(26) = 8;
	ExcitSpins(27) = 1;
	ExcitSpins(28) = 8;
	ExcitSpins(29) = 8;
	ExcitSpins(30) = 2;
	ExcitSpins(31) = 5;

	SetExcitationEnergiesPtr(&ExcitEnergies);
	SetExcitationSpinsPtr(&ExcitSpins);
}


G4NeutronEvaporationProbability::G4NeutronEvaporationProbability(const G4NeutronEvaporationProbability &right)
{
 G4Exception("G4NeutronEvaporationProbability::copy_constructor meant to not be accessable");
}




const G4NeutronEvaporationProbability & G4NeutronEvaporationProbability::
operator=(const G4NeutronEvaporationProbability &right)
{
  G4Exception("G4NeutronEvaporationProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4NeutronEvaporationProbability::operator==(const G4NeutronEvaporationProbability &right) const
{
  return false;
}

G4bool G4NeutronEvaporationProbability::operator!=(const G4NeutronEvaporationProbability &right) const
{
  return true;
}
