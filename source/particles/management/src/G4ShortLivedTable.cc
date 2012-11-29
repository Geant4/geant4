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
// $Id$
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	History: first implementation, based on object model of
//	27 June 1998  H.Kurashige
// ------------------------------------------------------------
//      added Remove()                  06 Nov.,98 H.Kurashige


#include "G4ShortLivedTable.hh"
#include "G4ParticleTable.hh"
#include "G4StateManager.hh"

#include "G4ios.hh"

G4ShortLivedTable::G4ShortLivedTable()
{
  fShortLivedList = new G4ShortLivedList();
}

G4ShortLivedTable::~G4ShortLivedTable()
{
  if (fShortLivedList ==0) return;

  //  No need to delete here because all particles are dynamic objects

  fShortLivedList->clear();
  delete fShortLivedList;
  fShortLivedList =0;
}

G4ShortLivedTable::G4ShortLivedTable(const G4ShortLivedTable & right)
{
  fShortLivedList = new G4ShortLivedList(*(right.fShortLivedList)); 
}

G4ShortLivedTable & G4ShortLivedTable::operator=(const G4ShortLivedTable &right)
{
  if (this != & right) {
    if (fShortLivedList !=0){
      fShortLivedList->clear();
      delete fShortLivedList;
      fShortLivedList = new G4ShortLivedList(*(right.fShortLivedList));
    } else {
      fShortLivedList = new G4ShortLivedList();
    }
  }
  return *this;
}

G4int G4ShortLivedTable::GetVerboseLevel() const
{
  return G4ParticleTable::GetParticleTable()->GetVerboseLevel();
}

G4bool G4ShortLivedTable::IsShortLived(const G4ParticleDefinition* particle) const
{
  return particle->IsShortLived();
}

void G4ShortLivedTable::clear()
{
  if (G4ParticleTable::GetParticleTable()->GetReadiness()) {
    G4Exception("G4ShortLivedTable::clear()",
		"PART116", JustWarning,
		"No effects because readyToUse is true.");    
    return;
  }

  fShortLivedList->clear();
}

void G4ShortLivedTable::Insert(const G4ParticleDefinition* particle)
{
  if (IsShortLived(particle)) {
    fShortLivedList->push_back(particle);
  } 
}

void G4ShortLivedTable::Remove(const G4ParticleDefinition* particle)
{
  if (G4ParticleTable::GetParticleTable()->GetReadiness()) {
    G4StateManager* pStateManager = G4StateManager::GetStateManager();
    G4ApplicationState currentState = pStateManager->GetCurrentState();
    if (currentState != G4State_PreInit) {
      G4String msg = "Request of removing ";
      msg += particle->GetParticleName();  
      msg += " has No effects other than Pre_Init";
      G4Exception("G4ShortLivedTable::Remove()",
		"PART117", JustWarning, msg);
      return;
    } else {
#ifdef G4VERBOSE
      if (GetVerboseLevel()>0){
	G4cout << particle->GetParticleName()
	       << " will be removed from the ShortLivedTable " << G4endl;
      }
#endif
    }
  }

  if (IsShortLived(particle)) {
    G4ShortLivedList::iterator idx;
    for (idx = fShortLivedList->begin(); idx!= fShortLivedList->end(); ++idx) {
      if ( particle == *idx) {
        fShortLivedList->erase(idx);
        break;
      }
    }
  } else {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1) {
      G4cout << "G4ShortLivedTable::Remove :" << particle->GetParticleName() ;
      G4cout << " is not short lived" << G4endl; 
    }
#endif
  }
}


void G4ShortLivedTable::DumpTable(const G4String &particle_name) const
{
  const G4ParticleDefinition* particle;

  G4ShortLivedList::iterator idx;
  for (idx = fShortLivedList->begin(); idx!= fShortLivedList->end(); ++idx) {
    particle = *idx;
    if (( particle_name == "ALL" ) || (particle_name == "all")){
      particle->DumpTable();
    } else if ( particle_name == particle->GetParticleName() ) {
      particle->DumpTable();
    }
  }
}












