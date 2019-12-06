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
/// \file ExN04MuonSD.cc
/// \brief Implementation of the ExN04MuonSD class
//

#include "ExN04MuonSD.hh"
#include "ExN04MuonHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
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

G4bool ExN04MuonSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
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

void ExN04MuonSD::EndOfEvent(G4HCofThisEvent*)
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

