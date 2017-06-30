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
/// \file electromagnetic/TestEm5/src/StackingAction.cc
/// \brief Implementation of the StackingAction class
//
// $Id: StackingAction.cc 104417 2017-05-30 08:30:48Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StackingAction.hh"

#include "Run.hh"
#include "EventAction.hh"
#include "HistoManager.hh"
#include "StackingMessenger.hh"

#include "G4RunManager.hh"
#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction(EventAction* EA)
 : G4UserStackingAction(), fEventAction(EA),
   fKillSecondary(0),fStackMessenger(0),fPhotoGamma(-1),fComptGamma(-1),
   fPhotoAuger(-1),fComptAuger(-1),fPixeGamma(-1),fPixeAuger(-1),
   fElectronDNAGamma(-1),fElectronDNAAuger(-1),fProtonDNAGamma(-1),
   fProtonDNAAuger(-1),fHydrogenDNAGamma(-1),fHydrogenDNAAuger(-1),
   fAlphaDNAGamma(-1),fAlphaDNAAuger(-1),fAlphaPlusDNAGamma(-1),
   fAlphaPlusDNAAuger(-1),fHeliumDNAGamma(-1),fHeliumDNAAuger(-1),
   fGenericIonDNAGamma(-1),fGenericIonDNAAuger(-1),fIDdefined(false)
{
  fStackMessenger = new StackingMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::~StackingAction()
{
  delete fStackMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  //keep primary particle
  if (aTrack->GetParentID() == 0) { return fUrgent; }

  if(!fIDdefined) {
    fIDdefined = true;
    fPhotoGamma = G4PhysicsModelCatalog::GetIndex("phot_fluo");
    fComptGamma = G4PhysicsModelCatalog::GetIndex("compt_fluo");
    fPhotoAuger = G4PhysicsModelCatalog::GetIndex("phot_auger");
    fComptAuger = G4PhysicsModelCatalog::GetIndex("compt_auger");
    fPixeGamma = G4PhysicsModelCatalog::GetIndex("gammaPIXE");
    fPixeAuger = G4PhysicsModelCatalog::GetIndex("e-PIXE");
    fElectronDNAGamma = 
     G4PhysicsModelCatalog::GetIndex("e-_G4DNAIonisation_fluo");
    fElectronDNAAuger = 
     G4PhysicsModelCatalog::GetIndex("e-_G4DNAIonisation_auger");
    fProtonDNAGamma = 
     G4PhysicsModelCatalog::GetIndex("proton_G4DNAIonisation_fluo");
    fProtonDNAAuger = 
     G4PhysicsModelCatalog::GetIndex("proton_G4DNAIonisation_auger");
    fHydrogenDNAGamma = 
     G4PhysicsModelCatalog::GetIndex("hydrogen_G4DNAIonisation_fluo");
    fHydrogenDNAAuger = 
     G4PhysicsModelCatalog::GetIndex("hydrogen_G4DNAIonisation_auger");
    fAlphaDNAGamma = 
     G4PhysicsModelCatalog::GetIndex("alpha_G4DNAIonisation_fluo");
    fAlphaDNAAuger = 
     G4PhysicsModelCatalog::GetIndex("alpha_G4DNAIonisation_auger");
    fAlphaPlusDNAGamma = 
     G4PhysicsModelCatalog::GetIndex("alpha+_G4DNAIonisation_fluo");
    fAlphaPlusDNAAuger = 
     G4PhysicsModelCatalog::GetIndex("alpha+_G4DNAIonisation_auger");
    fHeliumDNAGamma = 
     G4PhysicsModelCatalog::GetIndex("helium_G4DNAIonisation_fluo");
    fHeliumDNAAuger = 
     G4PhysicsModelCatalog::GetIndex("helium_G4DNAIonisation_auger");
    fGenericIonDNAGamma = 
     G4PhysicsModelCatalog::GetIndex("GenericIon_G4DNAIonisation_fluo");
    fGenericIonDNAAuger = 
     G4PhysicsModelCatalog::GetIndex("GenericIon_G4DNAIonisation_auger");
  }
  G4int idx = aTrack->GetCreatorModelID();

  //count secondary particles
    
  Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun()); 
  run->CountParticles(aTrack->GetDefinition());
  /*
  G4cout << "###StackingAction: new " 
         << aTrack->GetDefinition()->GetParticleName()
         << " E(MeV)= " << aTrack->GetKineticEnergy()
         << "  " << aTrack->GetMomentumDirection() << G4endl;
  */
  //
  //energy spectrum of secondaries
  //
  G4double energy = aTrack->GetKineticEnergy();
  G4double charge = aTrack->GetDefinition()->GetPDGCharge();

  if (charge != 0.) {
    analysisManager->FillH1(2,energy);
    analysisManager->FillH1(4,energy);
    if(idx == fPhotoAuger || idx == fComptAuger) {
      analysisManager->FillH1(50,energy);
      analysisManager->FillH1(52,energy);
    } else if(idx == fPixeAuger) {
      analysisManager->FillH1(54,energy);
      analysisManager->FillH1(56,energy);
    } else if(idx == fElectronDNAAuger || 
              idx == fProtonDNAAuger || 
              idx == fHydrogenDNAAuger || 
              idx == fAlphaDNAAuger || 
              idx == fAlphaPlusDNAAuger || 
              idx == fHeliumDNAAuger || 
              idx == fGenericIonDNAAuger) {
      analysisManager->FillH1(58,energy);
      analysisManager->FillH1(60,energy);
    }
  }

  if (aTrack->GetDefinition() == G4Gamma::Gamma()) {
    analysisManager->FillH1(3,energy);
    analysisManager->FillH1(5,energy);
    if(idx == fPhotoGamma || idx == fComptGamma) {
      analysisManager->FillH1(51,energy);
      analysisManager->FillH1(53,energy);
    } else if(idx == fPixeGamma) {
      analysisManager->FillH1(55,energy);
      analysisManager->FillH1(57,energy);
    } else if(idx == fElectronDNAGamma || 
              idx == fProtonDNAGamma || 
              idx == fHydrogenDNAGamma || 
              idx == fAlphaDNAGamma || 
              idx == fAlphaPlusDNAGamma || 
              idx == fHeliumDNAGamma || 
              idx == fGenericIonDNAGamma) {
      analysisManager->FillH1(59,energy);
      analysisManager->FillH1(61,energy);
    }
  }  

  //stack or delete secondaries
  G4ClassificationOfNewTrack status = fUrgent;
  if (fKillSecondary) {
    if (fKillSecondary == 1) {
     fEventAction->AddEnergy(energy);
     status = fKill;
    }  
    if (aTrack->GetDefinition() == G4Gamma::Gamma()) {
      status = fKill;
    }
  }
    
  return status;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
