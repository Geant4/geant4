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
// $Id: G4PSPassageCellFlux.cc 81087 2014-05-20 15:44:27Z gcosmo $
// GEANT4 tag $Name: geant4-09-04 $
//
// G4PSPassageCellFlux
#include "G4PSPassageCellFlux.hh"

#include "G4SystemOfUnits.hh"
#include "G4StepStatus.hh"
#include "G4Track.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPVParameterisation.hh"
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
// 2010-07-22   Introduce Unit specification.
// 2010-07-22   Add weighted option
// 
///////////////////////////////////////////////////////////////////////////////

G4PSPassageCellFlux::G4PSPassageCellFlux(G4String name, G4int depth)
  : G4VPrimitiveScorer(name,depth),HCID(-1),fCurrentTrkID(-1),fCellFlux(0),
    EvtMap(0),weighted(true)
{
    DefineUnitAndCategory();
    SetUnit("percm2");
}

G4PSPassageCellFlux::G4PSPassageCellFlux(G4String name, const G4String& unit,
					 G4int depth)
  : G4VPrimitiveScorer(name,depth),HCID(-1),fCurrentTrkID(-1),fCellFlux(0),
    EvtMap(0),weighted(true)
{
    DefineUnitAndCategory();
    SetUnit(unit);
}

G4PSPassageCellFlux::~G4PSPassageCellFlux()
{;}

G4bool G4PSPassageCellFlux::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{

  if ( IsPassed(aStep) ) {
    G4int idx = ((G4TouchableHistory*)
		 (aStep->GetPreStepPoint()->GetTouchable()))
                 ->GetReplicaNumber(indexDepth);
    G4double cubicVolume = ComputeVolume(aStep, idx);

    fCellFlux /= cubicVolume;
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
  if ( weighted ) trklength *= aStep->GetPreStepPoint()->GetWeight();

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
	   << "  cell flux : " << *(itr->second)/GetUnitValue()
	   << " [" << GetUnit()
	   << G4endl;
  }
}

void G4PSPassageCellFlux::SetUnit(const G4String& unit)
{
    CheckAndSetUnit(unit,"Per Unit Surface");
}

void G4PSPassageCellFlux::DefineUnitAndCategory(){
   // Per Unit Surface
   new G4UnitDefinition("percentimeter2","percm2","Per Unit Surface",(1./cm2));
   new G4UnitDefinition("permillimeter2","permm2","Per Unit Surface",(1./mm2));
   new G4UnitDefinition("permeter2","perm2","Per Unit Surface",(1./m2));
}


G4double G4PSPassageCellFlux::ComputeVolume(G4Step* aStep, G4int idx){

  G4VPhysicalVolume* physVol = aStep->GetPreStepPoint()->GetPhysicalVolume();
  G4VPVParameterisation* physParam = physVol->GetParameterisation();
  G4VSolid* solid = 0;
  if(physParam)
  { // for parameterized volume
    if(idx<0)
    {
      G4ExceptionDescription ED;
      ED << "Incorrect replica number --- GetReplicaNumber : " << idx << G4endl;
      G4Exception("G4PSPassageCellFlux::ComputeVolume","DetPS0013",JustWarning,ED);
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
