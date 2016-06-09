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
// $Id: G4PSEnergyDeposit.cc,v 1.3 2005/12/13 08:43:41 gunter Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// G4PSEnergyDeposit
#include "G4PSEnergyDeposit.hh"
#include "G4UnitsTable.hh"
////////////////////////////////////////////////////////////////////////////////
// Description:
//   This is a primitive scorer class for scoring energy deposit.
// 
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura.
// 
///////////////////////////////////////////////////////////////////////////////

G4PSEnergyDeposit::G4PSEnergyDeposit(G4String name, G4int depth)
  :G4VPrimitiveScorer(name,depth),HCID(-1)
{;}

G4PSEnergyDeposit::~G4PSEnergyDeposit()
{;}

G4bool G4PSEnergyDeposit::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if ( edep == 0. ) return FALSE;
  edep *= aStep->GetPreStepPoint()->GetWeight(); // (Particle Weight)
  G4int  index = GetIndex(aStep);
  EvtMap->add(index,edep);  
  return TRUE;
}

void G4PSEnergyDeposit::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(GetMultiFunctionalDetector()->GetName(),
				    GetName());
  if(HCID < 0) {HCID = GetCollectionID(0);}
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void G4PSEnergyDeposit::EndOfEvent(G4HCofThisEvent*)
{;}

void G4PSEnergyDeposit::clear()
{
  EvtMap->clear();
}

void G4PSEnergyDeposit::DrawAll()
{;}

void G4PSEnergyDeposit::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() << G4endl;
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
  for(; itr != EvtMap->GetMap()->end(); itr++) {
    G4cout << "  copy no.: " << itr->first
	   << "  energy deposit: " << G4BestUnit(*(itr->second),"Energy") 
	   << G4endl;
  }
}
