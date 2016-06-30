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
// $Id: G4VPrimitiveScorer.cc 94771 2015-12-09 09:44:05Z gcosmo $
//
// G4VPrimitiveScorer
#include "G4VPrimitiveScorer.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4UnitsTable.hh"

G4VPrimitiveScorer::G4VPrimitiveScorer(G4String name, G4int depth)
  :primitiveName(name),detector(nullptr),filter(nullptr),verboseLevel(0),indexDepth(depth),
   unitName("NoUnit"),unitValue(1.0),fNi(0),fNj(0),fNk(0)
{;} 

G4VPrimitiveScorer::~G4VPrimitiveScorer()
{;}

G4int G4VPrimitiveScorer::GetCollectionID(G4int)
{
  if(detector)
   return G4SDManager::GetSDMpointer()
    ->GetCollectionID(detector->GetName()+"/"+primitiveName); 
  else
   return -1;
}

void G4VPrimitiveScorer::Initialize(G4HCofThisEvent*)
{;}

void G4VPrimitiveScorer::EndOfEvent(G4HCofThisEvent*)
{;}

void G4VPrimitiveScorer::clear()
{;}

void G4VPrimitiveScorer::DrawAll()
{;}

void G4VPrimitiveScorer::PrintAll()
{;}

G4int G4VPrimitiveScorer::GetIndex(G4Step* aStep)
{
  G4StepPoint* preStep = aStep->GetPreStepPoint();
  G4TouchableHistory* th = (G4TouchableHistory*)(preStep->GetTouchable());
  return th->GetReplicaNumber(indexDepth);
}

void G4VPrimitiveScorer::CheckAndSetUnit(const G4String& unit,
					 const G4String& category){
    if ( G4UnitDefinition::GetCategory(unit) == category){
	unitName = unit;
	unitValue = G4UnitDefinition::GetValueOf(unit);
    } else {
	G4String msg = "Invalid unit ["+unit+"] (Current  unit is [" +GetUnit()+"] ) requested for " + GetName();
	G4Exception("G4VPrimitiveScorer::CheckAndSetUnit","Det0151",JustWarning,msg);
    }
}
