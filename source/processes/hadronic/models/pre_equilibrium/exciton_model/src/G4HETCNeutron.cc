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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#include "G4HETCNeutron.hh"

G4double G4HETCNeutron::K(const G4Fragment & aFragment)
{
  if (GetStage() != 1) return 1.0;
  // Number of protons in projectile
  G4double Pa = static_cast<G4int>(aFragment.GetParticleDefinition()->GetPDGCharge());
  // Number of neutrons in projectile 
  G4double Na = aFragment.GetParticleDefinition()->GetBaryonNumber();
  G4double TargetA = aFragment.GetA() - Na;
  G4double TargetZ = aFragment.GetZ() - Pa;
  Na -= Pa;
  G4double r = TargetZ/TargetA;
  
  G4double P = aFragment.GetNumberOfParticles();
  G4double H = aFragment.GetNumberOfHoles();
  
  G4double result = 0.0;
  if (P > 0)
    {
      result = (H*(1.0-r)+Na)/P;

      result /= (TargetA-TargetZ)/TargetA;
    }
  
  return std::max(0.0,result);
}


G4double G4HETCNeutron::GetKineticEnergy(const G4Fragment & aFragment)
{
  G4double H = aFragment.GetNumberOfHoles();
  G4double Pb = aFragment.GetNumberOfParticles() - GetA();
  G4double Nb = Pb + H;
  
  G4double Ab = std::max(0.0,(Pb*Pb+H*H+Pb-3*H)/4.0);
  G4double Emax = GetMaximalKineticEnergy() - Ab;
  
  G4double cut = GetBeta() / (GetBeta()+Emax/(Nb+1));
  G4double x(0.0);
  if (G4UniformRand() <= cut)
    {
      x = BetaRand(static_cast<G4int>(Nb),1);
    }
  else 
    {
      x = BetaRand(static_cast<G4int>(Nb),2);
    }

  return Emax * (1.0 - x);
}
