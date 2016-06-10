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
// $Id: G4PSCellFlux.cc 81087 2014-05-20 15:44:27Z gcosmo $
// GEANT4 tag $Name: geant4-09-04 $
//
// G4PSCellFlux
#include "G4PSCellFlux.hh"

#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4UnitsTable.hh"
///////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring cell flux.
//   The Cell Flux is defined by  a sum of track length divided 
//   by the geometry volume, where all of the tracks in the geometry 
//   are taken into account. 
//
//   If you want to score only tracks passing through the geometry volume,
//  please use G4PSPassageCellFlux.
//
//
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura.
// 2010-07-22   Introduce Unit specification.
// 2010-07-22   Add weighted option
// 
///////////////////////////////////////////////////////////////////////////////

G4PSCellFlux::G4PSCellFlux(G4String name, G4int depth)
    :G4VPrimitiveScorer(name,depth),HCID(-1),EvtMap(0),weighted(true)
{
    DefineUnitAndCategory();
    SetUnit("percm2");
    //verboseLevel = 10;
}

G4PSCellFlux::G4PSCellFlux(G4String name, const G4String& unit, G4int depth)
    :G4VPrimitiveScorer(name,depth),HCID(-1),EvtMap(0),weighted(true)
{
    DefineUnitAndCategory();
    SetUnit(unit);
}

G4PSCellFlux::~G4PSCellFlux()
{;}

G4bool G4PSCellFlux::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double stepLength = aStep->GetStepLength();
  if ( stepLength == 0. ) return FALSE;

  G4int idx = ((G4TouchableHistory*)
	       (aStep->GetPreStepPoint()->GetTouchable()))
               ->GetReplicaNumber(indexDepth);
  G4double cubicVolume = ComputeVolume(aStep, idx);

  G4double CellFlux = stepLength / cubicVolume;
  if (weighted) CellFlux *= aStep->GetPreStepPoint()->GetWeight(); 
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
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() <<G4endl; 
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
  for(; itr != EvtMap->GetMap()->end(); itr++) {
    G4cout << "  copy no.: " << itr->first
	   << "  cell flux : " << *(itr->second)/GetUnitValue() 
	   << " [" << GetUnit() << "]"
	   << G4endl;
  }
}

void G4PSCellFlux::SetUnit(const G4String& unit)
{
    CheckAndSetUnit(unit,"Per Unit Surface");
}

void G4PSCellFlux::DefineUnitAndCategory(){
   // Per Unit Surface
   new G4UnitDefinition("percentimeter2","percm2","Per Unit Surface",(1./cm2));
   new G4UnitDefinition("permillimeter2","permm2","Per Unit Surface",(1./mm2));
   new G4UnitDefinition("permeter2","perm2","Per Unit Surface",(1./m2));
}

G4double G4PSCellFlux::ComputeVolume(G4Step* aStep, G4int idx){

  G4VPhysicalVolume* physVol = aStep->GetPreStepPoint()->GetPhysicalVolume();
  G4VPVParameterisation* physParam = physVol->GetParameterisation();
  G4VSolid* solid = 0;
  if(physParam)
  { // for parameterized volume
    if(idx<0)
    {
      G4ExceptionDescription ED;
      ED << "Incorrect replica number --- GetReplicaNumber : " << idx << G4endl;
      G4Exception("G4PSCellFlux::ComputeVolume","DetPS0001",JustWarning,ED);
    }
    solid = physParam->ComputeSolid(idx, physVol);
    solid->ComputeDimensions(physParam,idx,physVol);
  }
  else
  { // for ordinary volume
    solid = physVol->GetLogicalVolume()->GetSolid();
  }
  
  return solid->GetCubicVolume();
}
