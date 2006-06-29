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
// $Id: RE01TrackerSD.cc,v 1.2 2006-06-29 17:44:32 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
