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
// $Id: G4ScoringManager.cc,v 1.4 2007-07-12 07:57:16 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4ScoringManager.hh"
#include "G4ScoringMessenger.hh"
#include "G4THitsMap.hh"

G4ScoringManager* G4ScoringManager::fSManager = 0;

G4ScoringManager* G4ScoringManager::GetScoringManager()
{
  if(!fSManager)
  {
    fSManager = new G4ScoringManager;
  }
  return fSManager;
}

G4ScoringManager* G4ScoringManager::GetScoringManagerIfExist()
{ return fSManager; }

G4ScoringManager::G4ScoringManager():verboseLevel(0)
{
  fMessenger = new G4ScoringMessenger(this);
}

G4ScoringManager::~G4ScoringManager()
{
  delete fMessenger;
  fSManager = 0;
}

void G4ScoringManager::Accumulate(G4VHitsCollection* map)
{
  G4String wName = map->GetName();
  G4VScoringMesh* sm = FindMesh(wName);
  if(sm == NULL) return;
  sm->Accumulate(static_cast<G4THitsMap<double>*>(map));
}

G4VScoringMesh* G4ScoringManager::FindMesh(G4String wName)
{
  G4VScoringMesh* sm = 0;
  for(MeshVecItr itr = fMeshVec.begin(); itr != fMeshVec.end(); itr++) {
    G4String smName = (*itr)->GetWorldName();
    if(wName == smName) {
      sm = *itr;
      break;
    }
  }
  return sm;
}
void G4ScoringManager::List() const
{
  G4cout << "G4ScoringManager" << G4endl;
}
