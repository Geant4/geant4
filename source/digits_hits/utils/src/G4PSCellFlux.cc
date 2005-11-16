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
// $Id: G4PSCellFlux.cc,v 1.1 2005-11-16 23:12:42 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4PSCellFlux
#include "G4PSCellFlux.hh"
#include "G4Track.hh"
#include "G4VSolid.hh"

///////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring cell flux.
//   The Cell Flux is defined by  a track length divided by a geometry
//   volume, where all of the tracks in the geometry are taken 
//  into account. 
//    Please use G4PSPassageCellFlux, if you want to score only 
//  tracks which passed through the volume.
//
//
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura.
// 
///////////////////////////////////////////////////////////////////////////////

G4PSCellFlux::G4PSCellFlux(G4String name, G4int depth)
    :G4VPrimitiveScorer(name,depth),HCID(-1)
{;}

G4PSCellFlux::~G4PSCellFlux()
{;}

G4bool G4PSCellFlux::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
    G4double stepLength = aStep->GetStepLength();
    if ( stepLength == 0. ) return FALSE;
    G4double volume     = aStep->GetPreStepPoint()->GetPhysicalVolume()
                         ->GetLogicalVolume()->GetSolid()->GetCubicVolume();
    G4double CellFlux = stepLength / volume ;
    CellFlux *= aStep->GetPreStepPoint()->GetWeight(); 
    G4int index = GetIndex(aStep);
    EvtMap->add(index,CellFlux);

  return TRUE;
}

void G4PSCellFlux::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(detector->GetName(),
				    GetName());
  if ( HCID < 0 ) HCID = GetCollectionID(0);
  HCE->AddHitsCollection(HCID,EvtMap);
}

void G4PSCellFlux::EndOfEvent(G4HCofThisEvent*)
{;}

void G4PSCellFlux::clear(){
  EvtMap->clear();
}

void G4PSCellFlux::DrawAll()
{;}

void G4PSCellFlux::PrintAll()
{
  G4cout << " PrimitiveSenstivity " << GetName() <<G4endl; 
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
  for(; itr != EvtMap->GetMap()->end(); itr++) {
    G4cout << "  copy no.: " << itr->first
	   << "  cell flux : " << *(itr->second) /(mm/mm3)
	   << G4endl;
  }
}

