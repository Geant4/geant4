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
//
// $Id: Tst50TrackerSD.cc,v 1.4 2003-05-17 18:11:54 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 17 May     2003 SG      first implemntation
// -------------------------------------------------------------------

#include "Tst50TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"

Tst50TrackerSD::Tst50TrackerSD(G4String name)
  :G4VSensitiveDetector(name)
{ 
  G4String HCname;
  collectionName.insert(HCname="trackerCollection");
  hitID = new G4int[500];
}

Tst50TrackerSD::~Tst50TrackerSD()
{
  delete [] hitID;
}

void Tst50TrackerSD::Initialize(G4HCofThisEvent* HCE)
{
  trackerCollection = new Tst50TrackerHitsCollection
    (SensitiveDetectorName,collectionName[0]); 
  for (G4int j=0;j<1;j++)
    {hitID [j]= -1;}; 
}

G4bool Tst50TrackerSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ROhist)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep==0.) return false;

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    
  G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 
  theTouchable->MoveUpHistory();     
    
  if ( hitID[0]==-1)
    { 
      Tst50TrackerHit* newHit = new Tst50TrackerHit();
      newHit->SetEdep(edep); 
      hitID[0] = trackerCollection->insert(newHit) - 1;
    }
  else
    { (*trackerCollection)[hitID[0]]->AddEnergy(edep); 
    }
 
  return true;
}

void Tst50TrackerSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
    { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection(HCID,trackerCollection);
}


