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
// by V. Lara

#include "G4HETCChargedFragment.hh"
#include "G4PreCompoundParameters.hh"


G4double G4HETCChargedFragment::
GetKineticEnergy(const G4Fragment & aFragment)
{
  // Number of protons in projectile
  G4double Pa = aFragment.GetParticleDefinition()->GetPDGCharge();
  // Number of neutrons in projectile 
  G4double Na = aFragment.GetParticleDefinition()->GetBaryonNumber();
  Na -= Pa;
  
  G4double Pb = aFragment.GetNumberOfParticles();
  G4double H = aFragment.GetNumberOfHoles();

  G4double Ab = std::max(0.0,(Pb*Pb+H*H+Pb-3*H)/4.0);
  G4double Emax = GetMaximalKineticEnergy() - Ab;

  G4double x = BetaRand(static_cast<G4int>(Pb+H),2);
  
  return Emax - (Emax-GetCoulombBarrier())*x;
}
