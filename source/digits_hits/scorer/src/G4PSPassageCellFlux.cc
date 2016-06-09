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
// $Id: G4PSPassageCellFlux.cc,v 1.1 2007/07/11 01:31:03 asaim Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// G4PSPassageCellFlux
#include "G4PSPassageCellFlux.hh"
#include "G4StepStatus.hh"
#include "G4Track.hh"
#include "G4VSolid.hh"
#include "G4UnitsTable.hh"
////////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring cell flux.
//   The Cell Flux is defined by  a track length divided by a geometry
//   volume, where only tracks passing through the geometry are taken 
//  into account. e.g. the unit of Cell Flux is mm/mm3.
//
//   If you want to score all tracks in the geometry volume,
//  please use G4PSCellFlux.
//
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura.
// 
///////////////////////////////////////////////////////////////////////////////

G4PSPassageCellFlux::G4PSPassageCellFlux(G4String name, G4int depth)
  :G4VPrimitiveScorer(name,depth),HCID(-1),fCurrentTrkID(-1),fCellFlux(0)
{;}

G4PSPassageCellFlux::~G4PSPassageCellFlux()
{;}

G4bool G4PSPassageCellFlux::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{

  if ( IsPassed(aStep) ) {
    fCellFlux /= 
      aStep->GetPreStepPoint()->GetPhysicalVolume()
      ->GetLogicalVolume()->GetSolid()->GetCubicVolume();
    G4int index = GetIndex(aStep);
    EvtMap->add(index,fCellFlux);
  }

  return TRUE;
}

G4bool G4PSPassageCellFlux::IsPassed(G4Step* aStep){
  G4bool Passed = FALSE;

  G4bool IsEnter = aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary;
  G4bool IsExit  = aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary;

  G4int  trkid  = aStep->GetTrack()->GetTrackID();
  G4double trklength  = aStep->GetStepLength();
  trklength *= aStep->GetPreStepPoint()->GetWeight();

  if ( IsEnter &&IsExit ){         // Passed at one step
    fCellFlux = trklength;         // Track length is absolutely given.
    Passed = TRUE;                 
  }else if ( IsEnter ){            // Enter a new geometry
    fCurrentTrkID = trkid;         // Resetting the current track.
    fCellFlux  = trklength;     
  }else if ( IsExit ){             // Exit a current geometry
    if ( fCurrentTrkID == trkid ) {// if the track is same as entered,
      fCellFlux  += trklength;     // add the track length to current one.
      Passed = TRUE;               
    }
  }else{                           // Inside geometry
    if ( fCurrentTrkID == trkid ){ // if the track is same as entered,
      fCellFlux  += trklength;     // adding the track length to current one.
    }
  }

  return Passed;
}

void G4PSPassageCellFlux::Initialize(G4HCofThisEvent* HCE)
{
  fCurrentTrkID = -1;

  EvtMap = new G4THitsMap<G4double>(detector->GetName(),
				    GetName());
  if ( HCID < 0 ) HCID = GetCollectionID(0);
  HCE->AddHitsCollection(HCID,EvtMap);

}

void G4PSPassageCellFlux::EndOfEvent(G4HCofThisEvent*)
{;}

void G4PSPassageCellFlux::clear(){
  EvtMap->clear();
}

void G4PSPassageCellFlux::DrawAll()
{;}

void G4PSPassageCellFlux::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() <<G4endl; 
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
  for(; itr != EvtMap->GetMap()->end(); itr++) {
    G4cout << "  copy no.: " << itr->first
	   << "  cell flux : " << *(itr->second)*cm*cm << " [cm^-2]"
	   << G4endl;
  }
}

