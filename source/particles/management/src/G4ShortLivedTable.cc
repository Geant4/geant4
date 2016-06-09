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
// $Id: G4ShortLivedTable.cc,v 1.14.2.1 2008/04/25 12:21:52 kurasige Exp $
// GEANT4 tag $Name: geant4-09-01-patch-03 $
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

#include "G4ios.hh"

G4ShortLivedTable::G4ShortLivedTable()
{
  fShortLivedList = new G4ShortLivedList();
}

G4ShortLivedTable::~G4ShortLivedTable()
{
  if (fShortLivedList ==0) return;

  //  No need to delete here because all particles are dynamic objects
  //   
  // remove all contents in the short lived List and delete all particles  
  //G4ShortLivedList::iterator i;
  //for (i = fShortLivedList->begin(); i!= fShortLivedList->end(); ++i) {
  //  delete (*i);
  //}

  fShortLivedList->clear();
  delete fShortLivedList;
  fShortLivedList =0;
}

G4int G4ShortLivedTable::GetVerboseLevel() const
{
  return G4ParticleTable::GetParticleTable()->GetVerboseLevel();
}

G4bool G4ShortLivedTable::IsShortLived(G4ParticleDefinition* particle) const
{
  return particle->IsShortLived();
}

void G4ShortLivedTable::Insert(G4ParticleDefinition* particle)
{
  if (IsShortLived(particle)) {
    fShortLivedList->push_back(particle);
  } 
}

void G4ShortLivedTable::Remove(G4ParticleDefinition* particle)
{
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
    if (GetVerboseLevel()>0) {
      G4cout << "G4ShortLivedTable::Remove :" << particle->GetParticleName() ;
      G4cout << " is not short lived" << G4endl; 
    }
#endif
  }
}


void G4ShortLivedTable::DumpTable(const G4String &particle_name) const
{
  G4ParticleDefinition* particle;

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












