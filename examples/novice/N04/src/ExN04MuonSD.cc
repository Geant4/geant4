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

#include "ExN04MuonSD.hh"
#include "ExN04MuonHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ios.hh"

ExN04MuonSD::ExN04MuonSD(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="muonCollection");
  positionResolution = 5*cm;
}

ExN04MuonSD::~ExN04MuonSD(){;}

void ExN04MuonSD::Initialize(G4HCofThisEvent*HCE)
{
  static int HCID = -1;
  muonCollection = new ExN04MuonHitsCollection
                   (SensitiveDetectorName,collectionName[0]); 
  if(HCID<0)
  { HCID = GetCollectionID(0); }
  HCE->AddHitsCollection(HCID,muonCollection);
}

G4bool ExN04MuonSD::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist)
{
  G4Track* aTrack = aStep->GetTrack();
  G4double edep = aStep->GetTotalEnergyDeposit();

  if(edep==0.) return true;

  ExN04MuonHit* aHit;
  int nHit = muonCollection->entries();
  G4ThreeVector hitpos = aStep->GetPreStepPoint()->GetPosition();
  for(int i=0;i<nHit;i++)
  {
    aHit = (*muonCollection)[i];
    G4ThreeVector pos = aHit->GetPos();
    G4double dist2 = sqr(pos.x()-hitpos.x())
                    +sqr(pos.y()-hitpos.y())+sqr(pos.z()-hitpos.z());
    if(dist2<=sqr(positionResolution))
    aHit->AddEdep(edep);
    return true;
  }

  aHit = new ExN04MuonHit();
  aHit->SetEdep( edep );
  aHit->SetPos( aStep->GetPreStepPoint()->GetPosition() );
  muonCollection->insert( aHit );

  return true;
}

void ExN04MuonSD::EndOfEvent(G4HCofThisEvent*HCE)
{;}

void ExN04MuonSD::clear()
{
} 

void ExN04MuonSD::DrawAll()
{
} 

void ExN04MuonSD::PrintAll()
{
} 

