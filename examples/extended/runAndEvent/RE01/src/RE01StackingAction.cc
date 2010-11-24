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
// $Id: RE01StackingAction.cc,v 1.3 2010-11-24 22:47:23 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//


#include "RE01StackingAction.hh"

#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ios.hh"

#include "RE01TrackInformation.hh"
#include "RE01CalorimeterHit.hh"

RE01StackingAction::RE01StackingAction()
:stage(0),trackerHitsColID(-1),calorimeterHitsColID(-1)
{ ; }

RE01StackingAction::~RE01StackingAction()
{ ; }

G4ClassificationOfNewTrack 
RE01StackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  G4ClassificationOfNewTrack classification = fUrgent;

  if(stage==0)
  {
    RE01TrackInformation* trackInfo;
    if(aTrack->GetTrackStatus()==fSuspend) // Track reached to calorimeter
    {
      trackInfo = (RE01TrackInformation*)(aTrack->GetUserInformation());
      trackInfo->SetTrackingStatus(0);
      trackInfo->SetSourceTrackInformation(aTrack);
//      G4cout << "Track " << aTrack->GetTrackID() << " (parentID " << aTrack->GetParentID() 
//             << ") has reached to calorimeter and has been suspended at " << aTrack->GetPosition() << G4endl;
      classification = fWaiting;
    }
    else if(aTrack->GetParentID()==0) // Primary particle
    {
      trackInfo = new RE01TrackInformation(aTrack);
      trackInfo->SetTrackingStatus(1);
      G4Track* theTrack = (G4Track*)aTrack;
      theTrack->SetUserInformation(trackInfo);
    }
  }
  return classification;
}

G4VHitsCollection* RE01StackingAction::GetCalCollection()
{
  G4SDManager* SDMan = G4SDManager::GetSDMpointer();
  G4RunManager* runMan = G4RunManager::GetRunManager();
  if(calorimeterHitsColID<0)
  { calorimeterHitsColID = SDMan->GetCollectionID("calCollection"); }
  if(calorimeterHitsColID>=0)
  {
    const G4Event* currentEvent = runMan->GetCurrentEvent();
    G4HCofThisEvent* HCE = currentEvent->GetHCofThisEvent();
    return HCE->GetHC(calorimeterHitsColID);
  }
  return 0;
}

void RE01StackingAction::NewStage()
{
G4cout << "+++++++++++ Stage " << stage << G4endl;
  if(stage==0)
  {
    // display trajetory information in the tracking region
    G4cout << G4endl;
    G4cout << "Tracks in tracking region have been processed. -- Stage 0 over." << G4endl;
    G4cout << G4endl;
  }
  else
  {
    // display calorimeter information caused by a "source" track in the tracker region
//    G4cout << G4endl;
//    G4cout << "Processing one shower originated by a source track in tracker region is over." << G4endl;
//    G4cout << G4endl;

    RE01CalorimeterHitsCollection* CHC = (RE01CalorimeterHitsCollection*)GetCalCollection();
    if(CHC)
    { 
      int n_hit = CHC->entries();
      G4double totE = 0;
      G4int n_hitByATrack = 0;
      for(int i=0;i<n_hit;i++)
      {
        G4double edepByATrack = (*CHC)[i]->GetEdepByATrack();
        if(edepByATrack>0.)
        {
          totE += edepByATrack;
          if(n_hitByATrack==0)
          { (*CHC)[i]->GetTrackInformation()->Print(); }
          n_hitByATrack++;
          G4cout << "Cell[" << (*CHC)[i]->GetZ() << "," << (*CHC)[i]->GetPhi() << "]    " 
                 << edepByATrack/GeV << " [GeV]" << G4endl;
          (*CHC)[i]->ClearEdepByATrack();
        }
      }
      if(n_hitByATrack>0)
      {
        G4cout << "###  Total energy deposition in calorimeter by a source track in "
               << n_hitByATrack << " cells : " << totE / GeV << " (GeV)" << G4endl;
        G4cout << G4endl;
      }
    }
  }

  if(stackManager->GetNUrgentTrack())
  {
    // Transfer all tracks in Urgent stack to Waiting stack, since all tracks
    // in Waiting stack have already been transfered to Urgent stack before
    // invokation of this method.
    stackManager->TransferStackedTracks(fUrgent,fWaiting);

    // Then, transfer only one track to Urgent stack.
    stackManager->TransferOneStackedTrack(fWaiting,fUrgent);

    stage++;
  }
}
    
void RE01StackingAction::PrepareNewEvent()
{ 
  stage = 0; 
}


