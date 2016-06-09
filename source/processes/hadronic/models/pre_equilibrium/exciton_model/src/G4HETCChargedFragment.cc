//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
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
