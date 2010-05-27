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
// $Id: MyTrackerSD.cc,v 1.6 2010-05-27 15:00:18 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"

#include "MyTrackerSD.hh"
#include "MyTrackerHit.hh"

MyTrackerSD::MyTrackerSD(G4String name)
:G4VSensitiveDetector(name)
{
  collectionName.insert("TrackerCollection");
}

MyTrackerSD::~MyTrackerSD(){;}

void MyTrackerSD::Initialize(G4HCofThisEvent*)
{
  TrackerCollection = new MyTrackerHitsCollection(SensitiveDetectorName,
						  collectionName[0]); 
}

G4bool MyTrackerSD::ProcessHits(G4Step*aStep,G4TouchableHistory*)
{
  G4ThreeVector hitPoint
   = ( aStep->GetPreStepPoint()->GetPosition()
     + aStep->GetPostStepPoint()->GetPosition() ) * 0.5;
  G4double edep = 1.0;

  if(verboseLevel>0)
  { G4cout << " New Tracker Hit at " << hitPoint << G4endl; }

  MyTrackerHit* newHit = new MyTrackerHit;
  newHit->SetEdep( edep );
  newHit->SetPos( hitPoint);
  TrackerCollection->insert(newHit);

  return true;
}

void MyTrackerSD::EndOfEvent(G4HCofThisEvent*HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection(HCID, TrackerCollection );
}

void MyTrackerSD::clear()
{
} 

void MyTrackerSD::DrawAll()
{
} 

void MyTrackerSD::PrintAll()
{
} 

