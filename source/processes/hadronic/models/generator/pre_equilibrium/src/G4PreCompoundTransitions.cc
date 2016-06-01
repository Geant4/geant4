// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// by V. Lara

#include "G4PreCompoundTransitions.hh"

G4PreCompoundTransitions::
G4PreCompoundTransitions(const G4Fragment & aFragment)
{
  // Fermi energy; method internal units are MeV
  const G4double FermiEnergy = 45.0;
  
  // Number of holes
  G4double H = aFragment.GetNumberOfHoles();
  // Number of Particles 
  G4double P = aFragment.GetNumberOfParticles();

  // Relative Energy (T_{rel})
  G4double RelativeEnergy = (8.0/5.0)*FermiEnergy + (aFragment.GetExcitationEnergy()/MeV)/(H+P);
  //    (aFragment.GetExcitationEnergy()/MeV)/aFragment.GetNumberOfExcitons();
  
  // Relative Velocity: 
  // <V_{rel}>^2
  G4double RelativeVelocitySqr = 2.0*RelativeEnergy/(G4Proton::Proton()->GetPDGMass()/MeV);
  // <V_{rel}>
  G4double RelativeVelocity = sqrt(RelativeVelocitySqr);

  // Proton-Proton Cross Section (in mbarn)
  G4double ppXSection = 10.63/RelativeVelocitySqr - 29.93/RelativeVelocity + 42.9;
  // Proton-Neutron Cross Section (in mbarn)
  G4double npXSection = 34.10/RelativeVelocitySqr - 82.20/RelativeVelocity + 82.2;

  // Averaged Cross Section: \sigma(V_{rel})
  G4double AveragedXSection = (ppXSection+npXSection)/2.0;

  // Fermi energy Relative energy ratio
  G4double FermiRelRatio = FermiEnergy/RelativeEnergy;

  // This factor is introduced to take into account the Pauli principle
  G4double PauliFactor = 1.0 - (7.0/5.0)*FermiRelRatio;
  if (FermiRelRatio > 0.5) PauliFactor += (2.0/5.0)*FermiRelRatio*pow(2.0 - (1.0/FermiRelRatio), 5.0/2.0);

  // Transition probability for \Delta n = +2
  TransitionProb1 = 0.00332*AveragedXSection*PauliFactor*sqrt(RelativeEnergy)/
    pow(1.2 + 1.0/(4.7*RelativeVelocity), 3.0);

  
  G4double GE = G4PreCompoundParameters::GetAddress()->GetLevelDensity()*
                aFragment.GetA()*aFragment.GetExcitationEnergy()/MeV;

  // Transition probability for \Delta n = -2 (at F(p,h) = 0)
  //  TransitionProb2 = max(0, (TransitionProb1*P*H*(P+H+1.0)*(P+H-2.0))/(GE*GE));
  TransitionProb2 = (TransitionProb1*P*H*(P+H+1.0)*(P+H-2.0))/(GE*GE);
  if (TransitionProb2 < 0.0) TransitionProb2 = 0.0; 

  // Transition probability for \Delta n = 0 (at F(p,h) = 0)
  TransitionProb3 = TransitionProb1*(P+H+1.0)*(P*(P-1.0)+4.0*P*H+H*(H-1.0))/((P+H)*GE);
  
  return;
}

const G4PreCompoundTransitions & G4PreCompoundTransitions::operator=(const G4PreCompoundTransitions &right)
{
    G4Exception("G4PreCompoundTransitions::operator= meant to not be accessable");
    return *this;
}


G4bool G4PreCompoundTransitions::operator==(const G4PreCompoundTransitions &right) const
{
  return false;
}

G4bool G4PreCompoundTransitions::operator!=(const G4PreCompoundTransitions &right) const
{
  return true;
}




G4int G4PreCompoundTransitions::GetDeltaNExciton()
{
  G4int result = 0;
  G4double ChosenTransition = G4UniformRand()*GetTotalProbability();
  if (ChosenTransition <= TransitionProb1) 
  {
    // Number of excitons is increased on \Delta n = +2
    result = 2;
  } 
  else if (ChosenTransition <= TransitionProb1+TransitionProb2) 
  {
    // Number of excitons is increased on \Delta n = -2
    result = -2;
  }
  return result;
}
