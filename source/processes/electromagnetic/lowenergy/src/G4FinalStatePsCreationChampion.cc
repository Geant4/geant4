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
// $Id: G4FinalStatePsCreationChampion.cc,v 1.2 2009-06-10 13:32:36 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// -------------------------------------------------------------------

#include "G4FinalStatePsCreationChampion.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4FinalStatePsCreationChampion::G4FinalStatePsCreationChampion()
{
  lowEnergyLimit = 10. * eV;
  highEnergyLimit = 5. * keV;

  G4ParticleDefinition* positronDef = G4Positron::PositronDefinition();

  if (positronDef == 0) G4Exception("G4FinalStatePsCreationChampion constructor: positron is not defined");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4FinalStatePsCreationChampion::~G4FinalStatePsCreationChampion()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4FinalStateProduct& G4FinalStatePsCreationChampion::GenerateFinalState(const G4Track& track, const G4Step& /* step */)
{
  product.Clear();

  const G4DynamicParticle* particle = track.GetDynamicParticle();

  G4double k = particle->GetKineticEnergy();

  if (k > lowEnergyLimit && k < highEnergyLimit)
  {
      G4ParticleDefinition* definition = particle->GetDefinition();
      G4ParticleMomentum primaryDirection = particle->GetMomentumDirection();

      G4int positroniumState = cross.RandomSelectState(k,definition);
      
      G4int ionizationShell = cross.RandomSelectShell(k,definition,positroniumState);
        
      G4double bindingEnergy = waterStructure.IonisationEnergy(ionizationShell);
      
      G4double creationEnergy = 13.6*eV/(2*(positroniumState+1)*(positroniumState+1)) ; // Ry = 13.6 * eV;
      
      G4double secondaryKinetic = k + creationEnergy - bindingEnergy; 

      product.KillPrimaryParticle();
      product.AddEnergyDeposit(bindingEnergy-creationEnergy);

      G4DNAGenericIonsManager * instance(G4DNAGenericIonsManager::Instance());
      
      if (secondaryKinetic<0) 
      {
        G4cout << "***** WARNING in G4FinalStatePsCreationChampion: Ps secondaryKinetic < 0 --> set to zero " << G4endl; 
        secondaryKinetic=0; 
      }
      
      if (positroniumState == 0) 
      product.AddSecondary 
       (new G4DynamicParticle(instance->GetIon("Ps-1s"),track.GetDynamicParticle()->GetMomentumDirection(),secondaryKinetic));
      
      if (positroniumState == 1) 
      product.AddSecondary 
      (new G4DynamicParticle(instance->GetIon("Ps-2s"),track.GetDynamicParticle()->GetMomentumDirection(),secondaryKinetic));

  }
  
  return product;
}

