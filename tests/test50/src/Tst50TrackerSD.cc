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
// $Id: Tst50TrackerSD.cc,v 1.6 2003-07-03 13:43:10 guatelli Exp $
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
#include "Tst50DetectorConstruction.hh"

Tst50TrackerSD::Tst50TrackerSD(G4String name, Tst50DetectorConstruction* det)
  :G4VSensitiveDetector(name),Detector(det)
{ 
  G4String HCname;
  collectionName.insert(HCname="Tst50Collection");
  tst50Collection = 0;
}

Tst50TrackerSD::~Tst50TrackerSD()
{
}

void Tst50TrackerSD::Initialize(G4HCofThisEvent*)
{
  tst50Collection = new Tst50TrackerHitsCollection
    (SensitiveDetectorName,collectionName[0]); 
}

G4bool Tst50TrackerSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep==0.) return false;

  Tst50TrackerHit* newHit = new Tst50TrackerHit();
 
  newHit->SetEdep(edep);
  tst50Collection->insert( newHit );
 
  return true;
}

void Tst50TrackerSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
    { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection(HCID,tst50Collection);
}


