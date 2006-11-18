//
// ********************************************************************
// * Disclaimer                                                       *
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
// $Id: G4PSCellCharge.cc,v 1.1 2006-11-18 18:20:40 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4PSCellCharge
#include "G4PSCellCharge.hh"
#include "G4Track.hh"


///////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring cell charge.
//   The Cell Charge is defined by  a sum of deposited charge inside the cell.
//
// Created: 2006-08-20  Tsukasa ASO
// 
///////////////////////////////////////////////////////////////////////////////

G4PSCellCharge::G4PSCellCharge(G4String name, G4int depth)
    :G4VPrimitiveScorer(name,depth),HCID(-1)
{;}

G4PSCellCharge::~G4PSCellCharge()
{;}

G4bool G4PSCellCharge::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{

    // Enter or First step of primary.
    if( aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary 
	|| ( aStep->GetTrack()->GetParentID() == 0 &&
	     aStep->GetTrack()->GetCurrentStepNumber() == 1 ) ){
	G4double CellCharge = aStep->GetPreStepPoint()->GetCharge();
	CellCharge *= aStep->GetPreStepPoint()->GetWeight();
	G4int index = GetIndex(aStep);
	EvtMap->add(index,CellCharge);
    }

    // Exit
    if( aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary){
	G4double CellCharge = aStep->GetPreStepPoint()->GetCharge();
	CellCharge *= aStep->GetPreStepPoint()->GetWeight();
	G4int index = GetIndex(aStep);
	CellCharge *= -1.0;
	EvtMap->add(index,CellCharge);
    }

  return TRUE;
}

void G4PSCellCharge::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(detector->GetName(),
				    GetName());
  if ( HCID < 0 ) HCID = GetCollectionID(0);
  HCE->AddHitsCollection(HCID,EvtMap);
}

void G4PSCellCharge::EndOfEvent(G4HCofThisEvent*)
{;}

void G4PSCellCharge::clear(){
  EvtMap->clear();
}

void G4PSCellCharge::DrawAll()
{;}

void G4PSCellCharge::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() <<G4endl; 
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
  for(; itr != EvtMap->GetMap()->end(); itr++) {
    G4cout << "  copy no.: " << itr->first
	   << "  cell charge : " << *(itr->second) << " [e]"
	   << G4endl;
  }
}

