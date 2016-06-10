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
// $Id: G4PSNofStep.cc 81087 2014-05-20 15:44:27Z gcosmo $
//
// G4PSNofStep
#include "G4PSNofStep.hh"
#include "G4UnitsTable.hh"

// (Description)
//   This is a primitive scorer class for scoring number of steps in the
//  Cell.
//
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura.
// 2010-07-22   Introduce Unit specification.
//

G4PSNofStep::G4PSNofStep(G4String name, G4int depth)
    :G4VPrimitiveScorer(name,depth),HCID(-1),EvtMap(0),boundaryFlag(false)
{
    SetUnit("");
}

G4PSNofStep::~G4PSNofStep()
{;}

G4bool G4PSNofStep::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  if ( boundaryFlag ) {
      if ( aStep->GetStepLength() == 0. ) return FALSE;
  }
  G4int  index = GetIndex(aStep);
  G4double val = 1.0;
  EvtMap->add(index,val);  
  return TRUE;
}

void G4PSNofStep::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(detector->GetName(),GetName());
  if(HCID < 0) {HCID = GetCollectionID(0);}
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void G4PSNofStep::EndOfEvent(G4HCofThisEvent*)
{;}

void G4PSNofStep::clear(){
  EvtMap->clear();
}

void G4PSNofStep::DrawAll()
{;}

void G4PSNofStep::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() << G4endl;
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
  for(; itr != EvtMap->GetMap()->end(); itr++) {
    G4cout << "  copy no.: " << itr->first
	   << "  num of step: " << *(itr->second)
	   << " [steps] "
	   << G4endl;
  }
}

void G4PSNofStep::SetUnit(const G4String& unit)
{
  if (unit == "" ){
    unitName = unit;
    unitValue = 1.0;
  }else{
      G4String msg = "Invalid unit ["+unit+"] (Current  unit is [" +GetUnit()+"] ) for " + GetName();
      G4Exception("G4PSNofStep::SetUnit","DetPS0011",JustWarning,msg);
  }

}


