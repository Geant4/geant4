//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4PreCompoundTransitions.cc,v 1.9 2001/08/01 17:08:32 hpw Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
// by V. Lara

#include "G4PreCompoundTransitions.hh"

#include "G4EvaporationLevelDensityParameter.hh"
#include "G4PairingCorrection.hh"


G4PreCompoundTransitions::
G4PreCompoundTransitions(const G4Fragment & aFragment)
{
  // Fermi energy
  const G4double FermiEnergy = G4PreCompoundParameters::GetAddress()->GetFermiEnergy();
  
	// Nuclear radius
  const G4double r0 = G4PreCompoundParameters::GetAddress()->GetTransitionsr0();
   
  // In order to calculate the level density parameter
  G4EvaporationLevelDensityParameter theLDP;

  // Number of holes
  G4double H = aFragment.GetNumberOfHoles();
  // Number of Particles 
  G4double P = aFragment.GetNumberOfParticles();
  // Number of Excitons 
  G4double N = P+H;

  // Nucleus 
  G4double A = aFragment.GetA();
  G4double Z = aFragment.GetZ();
  G4double U = aFragment.GetExcitationEnergy();

  // Relative Energy (T_{rel})
  G4double RelativeEnergy = (8.0/5.0)*FermiEnergy + U/N;
  
  // Relative Velocity: 
  // <V_{rel}>^2
  G4double RelativeVelocitySqr = 2.0*RelativeEnergy/proton_mass_c2;
  // <V_{rel}>
  G4double RelativeVelocity = sqrt(RelativeVelocitySqr);

  // Proton-Proton Cross Section
  G4double ppXSection = (10.63/RelativeVelocitySqr - 29.93/RelativeVelocity + 42.9)*millibarn;
  // Proton-Neutron Cross Section
  G4double npXSection = (34.10/RelativeVelocitySqr - 82.20/RelativeVelocity + 82.2)*millibarn;

  // Averaged Cross Section: \sigma(V_{rel})
  G4double AveragedXSection = (ppXSection+npXSection)/2.0;

  // Fermi relative energy ratio
  G4double FermiRelRatio = FermiEnergy/RelativeEnergy;

  // This factor is introduced to take into account the Pauli principle
  G4double PauliFactor = 1.0 - (7.0/5.0)*FermiRelRatio;
  if (FermiRelRatio > 0.5) PauliFactor += (2.0/5.0)*FermiRelRatio*pow(2.0 - (1.0/FermiRelRatio), 5.0/2.0);

  // Interaction volume 
  G4double Vint = (4.0/3.0)*pi*pow(2.0*r0 + hbarc/(proton_mass_c2*RelativeVelocity) , 3.0);

  // Transition probability for \Delta n = +2
  //  TransitionProb1 = 0.00332*AveragedXSection*PauliFactor*sqrt(RelativeEnergy)/
  //    pow(1.2 + 1.0/(4.7*RelativeVelocity), 3.0);
  TransitionProb1 = AveragedXSection*PauliFactor*sqrt(2.0*RelativeEnergy/proton_mass_c2)/Vint;
  if (TransitionProb1 < 0.0) TransitionProb1 = 0.0; 

  // g = 0.595aA; GE = g*E where E is Excitation Energy
  G4double a = theLDP.LevelDensityParameter(A,Z,U-G4PairingCorrection::GetPairingCorrection(A,Z));
  //  G4double a = G4PreCompoundParameters::GetAddress()->GetLevelDensity();
  G4double GE = 0.595*a*A*U;


  // F(p,h) = 0.25*(p^2 + h^2 + p - h) - 0.5*h
  G4double Fph = (P*P+H*H+P-H)/4.0 - H/2.0;
  // F(p+1,h+1)
  G4double Fph1 = Fph + N/2.0;
  // (n+1)/n ((g*E - F(p,h))/(g*E - F(p+1,h+1)))^(n+1)
  G4double ProbFactor = ((N+1.0)/N) * pow((GE-Fph)/(GE-Fph1),N+1.0);


  // Transition probability for \Delta n = -2 (at F(p,h) = 0)
  //  TransitionProb2 = max(0, (TransitionProb1*P*H*(P+H+1.0)*(P+H-2.0))/(GE*GE));
  //  TransitionProb2 = (TransitionProb1*P*H*(P+H+1.0)*(P+H-2.0))/(GE*GE);
  TransitionProb2 = TransitionProb1 * ProbFactor * (P*H*(N+1.0)*(N-2.0))/((GE-Fph)*(GE-Fph));
  if (TransitionProb2 < 0.0) TransitionProb2 = 0.0; 

  // Transition probability for \Delta n = 0 (at F(p,h) = 0)
  //  TransitionProb3 = TransitionProb1*(P+H+1.0)*(P*(P-1.0)+4.0*P*H+H*(H-1.0))/((P+H)*GE);
  TransitionProb3 = TransitionProb1 * ProbFactor * (P*(P-1.0) + 4.0*P*H + H*(H-1.0))/(GE-Fph);
  if (TransitionProb3 < 0.0) TransitionProb3 = 0.0; 
  

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




G4Fragment G4PreCompoundTransitions::PerformTransition(const G4Fragment & aFragment)
{
  G4Fragment result(aFragment);
  G4double ChosenTransition = G4UniformRand()*this->GetTotalProbability();
  G4int deltaN = 0;
  if (ChosenTransition <= TransitionProb1) 
    {
      // Number of excitons is increased on \Delta n = +2
      deltaN = 2;
    } 
  else if (ChosenTransition <= TransitionProb1+TransitionProb2) 
    {
      // Number of excitons is increased on \Delta n = -2
      deltaN = -2;
    }
  result.SetNumberOfParticles(result.GetNumberOfParticles()+deltaN/2);
  result.SetNumberOfHoles(result.GetNumberOfHoles()+deltaN/2); 
  // With weight Z/A, number of charged particles is decreased on +1
  if ((deltaN > 0 || result.GetNumberOfCharged() > 0) &&
      (G4UniformRand() <= result.GetZ()/result.GetA())){
    result.SetNumberOfCharged(result.GetNumberOfCharged()+deltaN/2);
  }

  // Number of charged can not be greater that number of particles
  if ( result.GetNumberOfParticles() < result.GetNumberOfCharged() ) {
    result.SetNumberOfCharged(result.GetNumberOfParticles());
  }
  
  return result;
}

