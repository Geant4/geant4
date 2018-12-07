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
/// \file RE05/src/RE05StackingAction.cc
/// \brief Implementation of the RE05StackingAction class
//

#include "RE05StackingAction.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "RE05StackingActionMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RE05StackingAction::RE05StackingAction()
 : G4UserStackingAction(),
   fTrkHits(0), fMuonHits(0), fMessenger(0), 
   fStage(0),
   fReqMuon(2),
   fReqIso(10),
   fAngRoI(30.0*deg)
{ 
  fMessenger = new RE05StackingActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RE05StackingAction::~RE05StackingAction()
{ delete fMessenger; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack 
RE05StackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  G4ClassificationOfNewTrack classification = fWaiting;
  switch(fStage)
  {
  case 0: // Stage 0 : Primary muons only
    if(aTrack->GetParentID()==0)
    {
      G4ParticleDefinition * particleType = aTrack->GetDefinition();
      if((particleType==G4MuonPlus::MuonPlusDefinition())
       ||(particleType==G4MuonMinus::MuonMinusDefinition()))
      { classification = fUrgent; }
    }
    break;

  case 1: // Stage 1 : Charged primaries only
          //           Suspended tracks will be sent to the waiting stack
    if(aTrack->GetParentID()!=0) { break; }
    if(aTrack->GetTrackStatus()==fSuspend) { break; }
    if(aTrack->GetDefinition()->GetPDGCharge()==0.) { break; }
    classification = fUrgent;
    break;

  default: // Stage 2 : Accept all primaries
           //           Accept all secondaries in RoI
           //           Kill secondaries outside RoI
    if(aTrack->GetParentID()==0)
    { 
      classification = fUrgent;
      break;
    }
    if((fAngRoI<0.)||InsideRoI(aTrack,fAngRoI))
    { 
      classification = fUrgent;
      break;
    }
    classification = fKill;
  }
  return classification;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool RE05StackingAction::InsideRoI(const G4Track * aTrack,G4double ang)
{
  if(!fMuonHits)
  { fMuonHits = (RE05MuonHitsCollection*)GetCollection("muonCollection"); }
  if(!fMuonHits)
  { G4cerr << "muonCollection NOT FOUND" << G4endl;
    return true; }

  G4int nhits = fMuonHits->entries();

  const G4ThreeVector trPos = aTrack->GetPosition();
  for(G4int i=0;i<nhits;i++)
  {
    G4ThreeVector muHitPos = (*fMuonHits)[i]->GetPos();
    G4double angl = muHitPos.angle(trPos);
    if(angl<ang) { return true; }
  }

  return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VHitsCollection* RE05StackingAction::GetCollection(G4String colName)
{
  G4SDManager* SDMan = G4SDManager::GetSDMpointer();
  G4RunManager* runMan = G4RunManager::GetRunManager();
  int colID = SDMan->GetCollectionID(colName);
  if(colID>=0)
  {
    const G4Event* currentEvent = runMan->GetCurrentEvent();
    G4HCofThisEvent* HCE = currentEvent->GetHCofThisEvent();
    return HCE->GetHC(colID);
  }
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE05StackingAction::NewStage()
{
  fStage++;
  G4int nhits;
  if(fStage==1)
  {
  // Stage 0->1 : check if at least "fReqMuon" hits on muon chamber
  //              otherwise abort current event
    if(!fMuonHits)
    { fMuonHits = (RE05MuonHitsCollection*)GetCollection("muonCollection"); }
    if(!fMuonHits)
    { G4cerr << "muonCollection NOT FOUND" << G4endl;
      return; }
    nhits = fMuonHits->entries();
////    G4cout << "Stage 0->1 : " << nhits << " hits found in the muon chamber."
////         << G4endl;
    if(nhits<fReqMuon)
    { 
      stackManager->clear();
////      G4cout << "++++++++ event aborted" << G4endl;
      return;
    }
    stackManager->ReClassify();
    return;
  }

  else if(fStage==2)
  {
  // Stage 1->2 : check the isolation of muon tracks
  //              at least "fReqIsoMuon" isolated muons
  //              otherwise abort current event.
  //              Isolation requires "fReqIso" or less hits
  //              (including own hits) in the RoI region
  //              in the tracker layers.
    nhits = fMuonHits->entries();
    if(!fTrkHits)
    { fTrkHits = (RE05TrackerHitsCollection*)GetCollection("trackerCollection"); }
    if(!fTrkHits)
    { G4cerr << "trackerCollection NOT FOUND" << G4endl;
      return; }
    G4int nTrkhits = fTrkHits->entries();
    G4int isoMuon = 0;
    for(G4int j=0;j<nhits;j++)
    {
      G4ThreeVector hitPos = (*fMuonHits)[j]->GetPos();
      G4int nhitIn = 0;
      for(G4int jj=0;(jj<nTrkhits)&&(nhitIn<=fReqIso);jj++)
      {
        G4ThreeVector trkhitPos = (*fTrkHits)[jj]->GetPos();
        if(trkhitPos.angle(hitPos)<fAngRoI) nhitIn++;
      }
      if(nhitIn<=fReqIso) isoMuon++;
    }
////    G4cout << "Stage 1->2 : " << isoMuon << " isolated muon found." << G4endl;
    if(isoMuon<fReqIsoMuon)
    {
      stackManager->clear();
////      G4cout << "++++++++ event aborted" << G4endl;
      return;
    }
    stackManager->ReClassify();
    return;
  }

  else
  {
  // Other fStage change : just re-classify
    stackManager->ReClassify();
  }
}
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE05StackingAction::PrepareNewEvent()
{ 
  fStage = 0; 
  fTrkHits = 0;
  fMuonHits = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
