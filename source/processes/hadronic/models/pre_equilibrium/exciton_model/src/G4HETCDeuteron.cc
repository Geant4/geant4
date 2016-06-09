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
#include "G4HETCDeuteron.hh"

G4double G4HETCDeuteron::K(const G4Fragment & aFragment)
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
  if (P > 1)
    {
      result = 2.0* (H*(H-1.0)*r*(r-1.0)+H*(Na*r+Pa*(1.0-r)) + Pa*Na)/(P*(P-1.0));

      result /= (TargetZ/TargetA)*((TargetA-TargetZ)/TargetA);
    }

  return std::max(0.0,result);

}
