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
///////////////////////////////////////////////////////////////////////////////
// File: CCalStackingAction.cc
// Description: Stacking action needed for the application
///////////////////////////////////////////////////////////////////////////////
#include "CCalStackingAction.hh"
#include "G4StackManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "CCaloSD.hh"
#include "CCalSDList.hh"
#include "G4RunManager.hh"
#include "G4Navigator.hh"

//#define debug
//#define ddebug

CCalStackingAction::CCalStackingAction()
  : fTimeLimit(10000*CLHEP::ns),isInitialized(false) 
{}

CCalStackingAction::~CCalStackingAction(){}

void CCalStackingAction::PrepareNewEvent(){

  if(!isInitialized)   initialize();
  stage = firstStage;
  nurgent = 0;
  acceptSecondaries = 1;
}

void CCalStackingAction::initialize(){

  isInitialized = true;
 
  numberOfSD = CCalSDList::getInstance()->getNumberOfCaloSD();
#ifdef debug
  G4cout << "CCalStackingAction look for " << numberOfSD 
         << " calorimeter-like SD" << G4endl;
#endif
  G4int i = 0;
  for (i=0; i<numberOfSD; i++) {
    G4String theName(CCalSDList::getInstance()->getCaloSDName(i));
    SDName[i] = theName;
#ifdef debug
    G4cout << "Found SD  name " << theName << G4endl;
#endif
    theCaloSD[i] = 0;
  }   

  G4SDManager* sd = G4SDManager::GetSDMpointerIfExist();
  if (sd != 0) {
    
    for (i=0; i<numberOfSD; i++){

      G4VSensitiveDetector* aSD = sd->FindSensitiveDetector(SDName[i]);
      if (aSD==0) {
#ifdef debug
        G4cout << "CCalStackingAction::initialize: No SD with name " << SDName[i]
               << " in this Setup " << G4endl;
#endif
      } else {
        theCaloSD[i] = dynamic_cast<CCaloSD*>(aSD);
        theCaloSD[i]->SetPrimaryID(0);
      }           
    }
#ifdef debug
    G4cout << "CCalStackingAction::initialize: Could not get SD Manager !" 
           << G4endl;
#endif
  }   
}

G4ClassificationOfNewTrack CCalStackingAction::ClassifyNewTrack(const G4Track* aTrack){

  G4ClassificationOfNewTrack classification=fKill;
  G4int parentID = aTrack->GetParentID();
#ifdef ddebug
  G4TrackStatus status = aTrack->GetTrackStatus();
  G4cout << "Classifying track " << aTrack->GetTrackID()
         << " with status " << aTrack->GetTrackStatus() << G4endl;  
#endif
    
  if (aTrack->GetGlobalTime() > fTimeLimit) {
#ifdef debug
    G4cout << "Kills particle " << aTrack->GetDefinition()->GetParticleName() 
           << " of energy " << aTrack->GetKineticEnergy()/MeV << " MeV" 
           << G4endl;
#endif
    return classification = fKill;
  }
    
  if (stage<end) {
    /////////////////
    /// PRIMARIES ///
    /////////////////
    if (parentID == 0 ) {
      if ( nurgent == 0) {
        nurgent++;
        classification = fUrgent;
        setPrimaryID(aTrack->GetTrackID());
      }
      else  classification = fWaiting;   
    }

    ///////////////////
    /// SECONDARIES ///
    ///////////////////
       
    if (parentID > 0) {
      if (acceptSecondaries == 1) {
        if (trackStartsInCalo(const_cast<G4Track *>(aTrack))!=0 )
          classification = fUrgent;
        else
          classification = fWaiting; 
      } else {
        if(nurgent == 0){                     
          nurgent++;
          classification = fUrgent;
          setPrimaryID(aTrack->GetTrackID());
        } else
          classification = fWaiting;        
      }       
    }
       
       
  } else 
    classification = G4UserStackingAction::ClassifyNewTrack(aTrack);

#ifdef ddebug
  G4cout << " returning classification= " << classification
         << " for track "<< aTrack->GetTrackID() << G4endl;
#endif
  return classification;

}


void CCalStackingAction::NewStage(){

#ifdef ddebug
  G4cout << "In NewStage with stage = " << stage << G4endl;
#endif
  if (stage <end) {
    nurgent = 0;                    
    setPrimaryID(0);
    acceptSecondaries = 0;
    stackManager->ReClassify();
    acceptSecondaries = 1;
    if (stackManager->GetNUrgentTrack() == 0) {
      stage = stageLevel(stage+1);
    }
        
  }
}

G4bool CCalStackingAction::trackStartsInCalo(const G4Track* ){

 /// This method should check that the secondary particle
 /// was produced inside the detector calorimeter and 
 /// really is part of the shower.
 /// If it has been produced before the calorimeter 
 /// for ex. Bremsstrahlung, it should be treated as a new
 /// particle producing a new shower.

 return true;
}

void CCalStackingAction::setPrimaryID(G4int id){
  
  for (G4int i=0; i<numberOfSD; i++){
    if(theCaloSD[i] != 0)theCaloSD[i]->SetPrimaryID(id);
  }

}
