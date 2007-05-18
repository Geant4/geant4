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
// $Id: G4PSTrackCounter.cc,v 1.1 2007-05-18 00:00:38 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4PSTrackCounter
#include "G4PSTrackCounter.hh"


///////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring number of tracks in a cell.
//
// Created: 2007-02-02  Tsukasa ASO, Akinori Kimura.
//
///////////////////////////////////////////////////////////////////////////////

G4PSTrackCounter::G4PSTrackCounter(G4String name, G4int direction, G4int depth)
  :G4VPrimitiveScorer(name,depth),HCID(-1),fDirection(direction),weighted(false)
{;}

G4PSTrackCounter::~G4PSTrackCounter()
{;}

G4bool G4PSTrackCounter::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{

  G4StepPoint* preStep = aStep->GetPreStepPoint();
  G4StepPoint* posStep = aStep->GetPostStepPoint();
  G4bool IsEnter = preStep->GetStepStatus()==fGeomBoundary;
  G4bool IsExit  = posStep->GetStepStatus()==fGeomBoundary;

  // Regular : count in prestep volume.
  G4int index = GetIndex(aStep);

  //  G4cout << " trk " << aStep->GetTrack()->GetTrackID()
  //	 << " index " << index << " In " << IsEnter << " Out " <<IsExit
  //	 << " PreStatus " << preStep->GetStepStatus() 
  //	 << " PostStatus " << posStep->GetStepStatus() 
  //	 << " w " << aStep->GetPreStepPoint()->GetWeight()
  // <<    G4endl;

  G4bool flag = FALSE;

  if ( IsEnter && fDirection == fCurrent_In ) flag = TRUE;
  else if ( IsExit  && fDirection == fCurrent_Out ) flag = TRUE;
  else if ( (IsExit||IsEnter) && fDirection == fCurrent_InOut  ) flag = TRUE;

  if ( flag ){
    //G4cout << " Count " << G4endl;
      G4double val = 1.0;
      if(weighted) val *= aStep->GetPreStepPoint()->GetWeight();
      EvtMap->add(index,val);  
  }

  return TRUE;
}

void G4PSTrackCounter::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(detector->GetName(),GetName());
  if(HCID < 0) {HCID = GetCollectionID(0);}
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
  //G4cout <<   "----- PS Initialize " << G4endl;
}

void G4PSTrackCounter::EndOfEvent(G4HCofThisEvent*)
{;}

void G4PSTrackCounter::clear(){
  EvtMap->clear();
}

void G4PSTrackCounter::DrawAll()
{;}

void G4PSTrackCounter::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() << G4endl;
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
  for(; itr != EvtMap->GetMap()->end(); itr++) {
    G4cout << "  copy no.: " << itr->first
	   << "  track count: " << *(itr->second)
	   << G4endl;
  }
}
