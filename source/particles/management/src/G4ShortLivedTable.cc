// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ShortLivedTable.cc,v 1.5 1999-10-29 05:34:28 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	For information related to this code contact:
//	CERN, IT Division, ASD Group
//	History: first implementation, based on object model of
//	27 June 1998  H.Kurashige
// ------------------------------------------------------------
//      added Remove()                  06 Nov.,98 H.Kurashige


#include "G4ShortLivedTable.hh"
#include "G4ParticleTable.hh"

#include "G4ios.hh"

#ifdef WIN32
#  include <Strstrea.h>
#else
#  include <strstream.h>
#endif


G4ShortLivedTable::G4ShortLivedTable()
{
  fShortLivedList = new G4ShortLivedList();
}

G4ShortLivedTable::~G4ShortLivedTable()
{
  // remove all contents in the short lived List and delete all particles  
#ifdef G4USE_STL_MAP
  G4ShortLivedList::iterator i;
  for (i = fShortLivedList->begin(); i!= fShortLivedList->end(); ++i) {
    delete (*i);
  }
  fShortLivedList->clear();
#else
  fShortLivedList->clearAndDestroy();
#endif
  delete fShortLivedList;
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
#ifdef G4USE_STL_MAP
    fShortLivedList->push_back(particle);
#else
    fShortLivedList->insert(particle);
#endif
  } else {
    //#ifdef G4VERBOSE
    //if (GetVerboseLevel()>0) {
    //  G4cout << "G4ShortLivedTable::Insert :" << particle->GetParticleName() ;
    //  G4cout << " is not short lived" << endl; 
    //}
    //#endif
  }
}

void G4ShortLivedTable::Remove(G4ParticleDefinition* particle)
{
  if (IsShortLived(particle)) {
#ifdef G4USE_STL_MAP
    G4ShortLivedList::iterator idx;
    for (idx = fShortLivedList->begin(); idx!= fShortLivedList->end(); ++idx) {
      if ( particle == *idx) {
        fShortLivedList->erase(idx);
      }
    }
#else
    fShortLivedList->remove(particle);
#endif

  } else {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4ShortLivedTable::Remove :" << particle->GetParticleName() ;
      G4cout << " is not short lived" << endl; 
    }
#endif
  }
}


void G4ShortLivedTable::DumpTable(const G4String &particle_name) const
{
  G4ParticleDefinition* particle;

#ifdef G4USE_STL_MAP
  G4ShortLivedList::iterator idx;
  for (idx = fShortLivedList->begin(); idx!= fShortLivedList->end(); ++idx) {
    particle = *idx;
#else
  for (G4int idx= 0; idx < fShortLivedList->entries() ; idx++) {
    particle = (*fShortLivedList)(idx);
#endif 

    if (( particle_name == "ALL" ) || (particle_name == "all")){
      particle->DumpTable();
    } else if ( particle_name == particle->GetParticleName() ) {
      particle->DumpTable();
    }
  }
}

