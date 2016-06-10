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
// $Id: G4PSDoseDeposit.cc 81087 2014-05-20 15:44:27Z gcosmo $
// GEANT4 tag $Name: geant4-09-04 $
//
// G4PSDoseDeposit
#include "G4PSDoseDeposit.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4UnitsTable.hh"

////////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring only energy deposit.
// 
//
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura.
// 2010-07-22   Introduce Unit specification.
// 
///////////////////////////////////////////////////////////////////////////////

G4PSDoseDeposit::G4PSDoseDeposit(G4String name, G4int depth)
  :G4VPrimitiveScorer(name,depth),HCID(-1),EvtMap(0)
{
    SetUnit("Gy");
}

G4PSDoseDeposit::G4PSDoseDeposit(G4String name, const G4String& unit,
				 G4int depth)
  :G4VPrimitiveScorer(name,depth),HCID(-1),EvtMap(0)
{
    SetUnit(unit);
}

G4PSDoseDeposit::~G4PSDoseDeposit()
{;}

G4bool G4PSDoseDeposit::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if ( edep == 0. ) return FALSE;

  G4int idx = ((G4TouchableHistory*)
	       (aStep->GetPreStepPoint()->GetTouchable()))
               ->GetReplicaNumber(indexDepth);
  G4double cubicVolume = ComputeVolume(aStep, idx);


  G4double density = aStep->GetTrack()->GetStep()->GetPreStepPoint()->GetMaterial()->GetDensity();
  G4double dose    = edep / ( density * cubicVolume );
  dose *= aStep->GetPreStepPoint()->GetWeight(); 
  G4int  index = GetIndex(aStep);
  EvtMap->add(index,dose);  
  return TRUE;
}

void G4PSDoseDeposit::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(GetMultiFunctionalDetector()->GetName(),
				    GetName());
  if(HCID < 0) {HCID = GetCollectionID(0);}
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void G4PSDoseDeposit::EndOfEvent(G4HCofThisEvent*)
{;}

void G4PSDoseDeposit::clear()
{
  EvtMap->clear();
}

void G4PSDoseDeposit::DrawAll()
{;}

void G4PSDoseDeposit::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() << G4endl;
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
  for(; itr != EvtMap->GetMap()->end(); itr++) {
    G4cout << "  copy no.: " << itr->first
	   << "  dose deposit: " 
	   << *(itr->second)/GetUnitValue()
	   << " ["<<GetUnit() <<"]"
	   << G4endl;
  }
}

void G4PSDoseDeposit::SetUnit(const G4String& unit)
{
	CheckAndSetUnit(unit,"Dose");
}

G4double G4PSDoseDeposit::ComputeVolume(G4Step* aStep, G4int idx){

  G4VPhysicalVolume* physVol = aStep->GetPreStepPoint()->GetPhysicalVolume();
  G4VPVParameterisation* physParam = physVol->GetParameterisation();
  G4VSolid* solid = 0;
  if(physParam)
  { // for parameterized volume
    if(idx<0)
    {
      G4ExceptionDescription ED;
      ED << "Incorrect replica number --- GetReplicaNumber : " << idx << G4endl;
      G4Exception("G4PSDoseDeposit::ComputeVolume","DetPS0004",JustWarning,ED);
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



