// This code implementation is the intellectual property of
// the GEANT4 collaboration.
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


G4double G4FissionProbability::EmissionProbability(const G4Fragment & fragment, const G4double photonExcitation)
  //
{
  G4double A = fragment.GetA();
  G4double Z = fragment.GetZ();
  G4double U = fragment.GetExcitationEnergy();
  G4double SystemEntropy = 2.0*sqrt((theEvapLDP.LevelDensityParameter(A,Z,U)/(1./MeV))*A*U/MeV);

  // Compute integrated probability of fission channel
  if (theChannel->GetMaximalKineticEnergy() <= 0.0) return 0.0;

  G4double Q1 = 2.0*sqrt((theFissLDP.LevelDensityParameter(A,Z,U)/(1./MeV))*A*
			 theChannel->GetMaximalKineticEnergy()/MeV);

  G4double Q2 = 1./(4.0*pi);

  //  G4double Tfis = 21.e-6*940.0;


  //return min(Tfis,(Q2/((theFissLDP.LevelDensityParameter(A,Z,U)/(1./MeV))*A))*
  //	     ((Q1-1.0)*exp(Q1-SystemEntropy)+exp(-SystemEntropy)));


   return (Q2/((theFissLDP.LevelDensityParameter(A,Z,U)/(1./MeV))*A))*
     ((Q1-1.0)*exp(Q1-SystemEntropy)+exp(-SystemEntropy));
}


