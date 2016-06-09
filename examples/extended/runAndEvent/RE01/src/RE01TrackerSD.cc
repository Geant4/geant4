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
// $Id: RE01TrackerSD.cc,v 1.1 2004/11/26 07:37:43 asaim Exp $
// GEANT4 tag $Name: geant4-07-01 $
//



#include "RE01TrackerSD.hh"
#include "RE01TrackerHit.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "RE01TrackInformation.hh"

RE01TrackerSD::RE01TrackerSD(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="trackerCollection");
}

RE01TrackerSD::~RE01TrackerSD(){;}

void RE01TrackerSD::Initialize(G4HCofThisEvent* HCE)
{
  static int HCID = -1;
  trackerCollection = new RE01TrackerHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 
  if(HCID<0)
  { HCID = GetCollectionID(0); }
  HCE->AddHitsCollection(HCID,trackerCollection);
}

G4bool RE01TrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep==0.) return false;

  RE01TrackerHit* newHit = new RE01TrackerHit();
  newHit->SetEdep( edep );
  newHit->SetPos( aStep->GetPreStepPoint()->GetPosition() );
  RE01TrackInformation* trackInfo = (RE01TrackInformation*)(aStep->GetTrack()->GetUserInformation());
  if(trackInfo->GetTrackingStatus()>0)
  { newHit->SetTrackID( aStep->GetTrack()->GetTrackID() ); }
  else
  { newHit->SetTrackID( -1 ); }
  trackerCollection->insert( newHit );

  return true;
}

void RE01TrackerSD::EndOfEvent(G4HCofThisEvent*)
{
}

void RE01TrackerSD::clear()
{
} 

void RE01TrackerSD::DrawAll()
{
} 

void RE01TrackerSD::PrintAll()
{
} 
