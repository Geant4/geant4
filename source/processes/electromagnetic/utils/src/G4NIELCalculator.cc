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
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4NIELCalculator
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 30.05.2019
//
// Modifications:
//
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4NIELCalculator.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4LossTableManager.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VProcess.hh"
#include "G4Track.hh"
#include "Randomize.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4NIELCalculator::G4NIELCalculator(G4VEmModel* mod, G4int verb)
  : fModel(mod), fVerbose(verb)
{
  G4LossTableManager::Instance()->SetNIELCalculator(this);
  if(fVerbose > 0) {
    G4cout << "G4NIELCalculator: is created with the model <" 
           << fModel->GetName() << ">" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4NIELCalculator::AddEmModel(G4VEmModel* mod) 
{
  if(mod && mod != fModel) { 
    fModel = mod; 
    if(fVerbose > 0) {
      G4cout << "G4NIELCalculator: new model <" << fModel->GetName()
             << "> is added" << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4NIELCalculator::Initialise() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4NIELCalculator::ComputeNIEL(const G4Step* step)
{
  G4double niel = 0.0;
  G4double T2 = step->GetPostStepPoint()->GetKineticEnergy();
  if(fModel && T2 > 0.) {
    const G4Track* track = step->GetTrack();
    const G4ParticleDefinition* part = track->GetParticleDefinition();
    G4double length = step->GetStepLength(); 

    if(length > 0.0 && part->GetPDGMass() > 100*CLHEP::MeV) {

      // primary
      G4double T1= step->GetPreStepPoint()->GetKineticEnergy();
      G4double T = 0.5*(T1 + T2);
      const G4MaterialCutsCouple* couple = 
	step->GetPreStepPoint()->GetMaterialCutsCouple(); 
      niel = length*fModel->ComputeDEDXPerVolume(couple->GetMaterial(),part,T);
      niel = std::min(niel, T1);
    }
  }
  return niel;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4NIELCalculator::RecoilEnergy(const G4Step* step)
{
  G4double erec = 0.0;
  const std::vector<const G4Track*>* sec = step->GetSecondaryInCurrentStep(); 

  if(sec) {
    for(auto track : *sec) {
      const G4ParticleDefinition* part = track->GetParticleDefinition();
      if(part->IsGeneralIon()) {
        erec += track->GetKineticEnergy();
      }
    }
  }
  return erec;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
