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
// $Id: G4PSPassageCellCurrent.cc 81087 2014-05-20 15:44:27Z gcosmo $
//
// G4PSPassageCellCurrent
#include "G4PSPassageCellCurrent.hh"
#include "G4StepStatus.hh"
#include "G4Track.hh"
#include "G4VSolid.hh"
#include "G4UnitsTable.hh"
////////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring number of tracks,
//   where only tracks passing through the geometry are taken 
//  into account.
//
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura.
// 2010-07-22   Introduce Unit specification.
// 2010-07-22   Add weighted option
// 
///////////////////////////////////////////////////////////////////////////////

G4PSPassageCellCurrent::G4PSPassageCellCurrent(G4String name, G4int depth)
    :G4VPrimitiveScorer(name,depth),HCID(-1),fCurrentTrkID(-1),fCurrent(0),
     EvtMap(0),weighted(true)
{
    SetUnit("");
}

G4PSPassageCellCurrent::~G4PSPassageCellCurrent()
{;}

G4bool G4PSPassageCellCurrent::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{

  if ( IsPassed(aStep) ) {
    if(weighted) fCurrent = aStep->GetPreStepPoint()->GetWeight();
    G4int index = GetIndex(aStep);
    EvtMap->add(index,fCurrent);
  }

  return TRUE;
}

G4bool G4PSPassageCellCurrent::IsPassed(G4Step* aStep){
  G4bool Passed = FALSE;

  G4bool IsEnter = aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary;
  G4bool IsExit  = aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary;

  G4int  trkid  = aStep->GetTrack()->GetTrackID();

  if ( IsEnter &&IsExit ){         // Passed at one step
    Passed = TRUE;                 
  }else if ( IsEnter ){            // Enter a new geometry
    fCurrentTrkID = trkid;         // Resetting the current track.
  }else if ( IsExit ){             // Exit a current geometry
    if ( fCurrentTrkID == trkid ) {
      Passed = TRUE;               // if the track is same as entered.
    }
  }else{                           // Inside geometry
    if ( fCurrentTrkID == trkid ){ // Adding the track length to current one ,
    }
  }
  return Passed;
}

void G4PSPassageCellCurrent::Initialize(G4HCofThisEvent* HCE)
{
  fCurrentTrkID = -1;

  EvtMap = new G4THitsMap<G4double>(detector->GetName(),
				    GetName());
  if ( HCID < 0 ) HCID = GetCollectionID(0);
  HCE->AddHitsCollection(HCID,EvtMap);

}

void G4PSPassageCellCurrent::EndOfEvent(G4HCofThisEvent*)
{
}

void G4PSPassageCellCurrent::clear(){
  EvtMap->clear();
}

void G4PSPassageCellCurrent::DrawAll()
{;}

void G4PSPassageCellCurrent::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() <<G4endl; 
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
  for(; itr != EvtMap->GetMap()->end(); itr++) {
    G4cout << "  copy no.: " << itr->first
	   << "  cell current : " << *(itr->second)
	   << " [tracks] "
	   << G4endl;
  }
}

void G4PSPassageCellCurrent::SetUnit(const G4String& unit)
{
  if (unit == "" ){
    unitName = unit;
    unitValue = 1.0;
  }else{
      G4String msg = "Invalid unit ["+unit+"] (Current  unit is [" +GetUnit()+"] ) for " + GetName();
    G4Exception("G4PSPassageCellCurrent::SetUnit","DetPS0012",JustWarning,msg);
  }
}


