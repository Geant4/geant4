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
// $Id: G4PSTrackLength.cc,v 1.6 2007/05/18 00:00:38 asaim Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
// G4PSTrackLength
#include "G4PSTrackLength.hh"
#include "G4UnitsTable.hh"
////////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring sum of track length.
// 
//
// Created: 2007-02-02  Tsukasa ASO, Akinori Kimura.
// 
///////////////////////////////////////////////////////////////////////////////

G4PSTrackLength::G4PSTrackLength(G4String name, G4int depth)
    :G4VPrimitiveScorer(name,depth),HCID(-1),weighted(false),multiplyKinE(false),
     divideByVelocity(false)
{;}

G4PSTrackLength::~G4PSTrackLength()
{;}

G4bool G4PSTrackLength::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double trklength  = aStep->GetStepLength();
  if ( trklength == 0. ) return FALSE;
  if(weighted) trklength *= aStep->GetPreStepPoint()->GetWeight();
  if(multiplyKinE) trklength *= aStep->GetPreStepPoint()->GetKineticEnergy();
  if(divideByVelocity) trklength /= aStep->GetPreStepPoint()->GetVelocity();
  G4int  index = GetIndex(aStep);
  EvtMap->add(index,trklength);  
  return TRUE;
}

void G4PSTrackLength::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(detector->GetName(),GetName());
  if(HCID < 0) {HCID = GetCollectionID(0);}
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*)EvtMap);
}

void G4PSTrackLength::EndOfEvent(G4HCofThisEvent*)
{;}

void G4PSTrackLength::clear(){
  EvtMap->clear();
}

void G4PSTrackLength::DrawAll()
{;}

void G4PSTrackLength::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() << G4endl;
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  std::map<G4int,G4double*>::iterator itr = EvtMap->GetMap()->begin();
  for(; itr != EvtMap->GetMap()->end(); itr++) {
      G4cout << "  copy no.: " << itr->first << "  track length: " ;
      if ( multiplyKinE && !divideByVelocity ){
	  G4cout << *(itr->second)/(mm*MeV) <<" mm*MeV";
      } else if ( !multiplyKinE && divideByVelocity ){
	  G4cout << *(itr->second)*second <<" /second";
      } else if ( multiplyKinE && divideByVelocity) {
          G4cout << *(itr->second)/MeV*second <<" MeV/second";
      } else {
	  G4cout  << G4BestUnit(*(itr->second),"Length");
      }
      G4cout << G4endl;
  }
}

