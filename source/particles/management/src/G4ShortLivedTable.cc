// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ShortLivedTable.cc,v 1.2 1999-04-14 10:28:31 kurasige Exp $
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
  fShortLivedList->clearAndDestroy();
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
    fShortLivedList->insert(particle);
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
    fShortLivedList->remove(particle);
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
  for (G4int idx= 0; idx < fShortLivedList->entries() ; idx++) {
    if (( particle_name == "ALL" ) || (particle_name == "all")){
      ((*fShortLivedList)(idx))->DumpTable();
    } else if ( particle_name == ((*fShortLivedList)(idx))->GetParticleName() ) {
      ((*fShortLivedList)(idx))->DumpTable();
    }
  }
}

