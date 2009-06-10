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
// $Id: G4CrossSectionPsCreationChampion.cc,v 1.2 2009-06-10 13:32:36 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// -------------------------------------------------------------------

#include "G4CrossSectionPsCreationChampion.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CrossSectionPsCreationChampion::G4CrossSectionPsCreationChampion()
{
  lowEnergyLimit = 10. * eV;
  highEnergyLimit = 5. * keV;

  G4ParticleDefinition* positronDef = G4Positron::PositronDefinition();
  if (positronDef ==0 ) G4Exception("G4CrossSectionPsCreationchampion constructor: positron is not defined");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4CrossSectionPsCreationChampion::~G4CrossSectionPsCreationChampion()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4CrossSectionPsCreationChampion::CrossSection(const G4Track& track )
{
  const G4DynamicParticle* particle = track.GetDynamicParticle();
  G4double k = particle->GetKineticEnergy();
  const G4ParticleDefinition* particleDefinition = track.GetDefinition();
  
  if (particleDefinition != G4Positron::PositronDefinition())
     G4Exception("G4CrossSectionPsCreationChampion::CrossSection: attempting to calculate cross section for wrong particle");

  G4double totalCrossSection = 0.;

  if (k > lowEnergyLimit && k < highEnergyLimit)
  {
    totalCrossSection = partialCrossSection.Sum(k,particleDefinition);
  }
  return totalCrossSection;
}
