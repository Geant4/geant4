// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ShortLivedTable.hh,v 1.4 1999-10-29 05:34:22 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, IT Division, ASD group
//	History: first implementation, 
//      based on object model of June 27, 98 H.Kurashige
// ------------------------------------------------------------
//      added Remove()                  06 Nov.,98 H.Kurashige

#ifndef G4ShortLivedTable_h
#define G4ShortLivedTable_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"


#ifdef G4USE_STL_MAP
#include "g4std/vector"
#else
#include <rw/tpordvec.h>
#endif 

class G4ParticleTable;

class G4ShortLivedTable
{
 // Class Description
 //   G4ShortLivedTable is the table of pointer to G4ParticleDefinition
 //   In G4ShortLivedTable, each G4ParticleDefinition pointer is stored
 //

 public:
#ifdef G4USE_STL_MAP
   // Use STL Vector as list of shortlives
   typedef vector<G4ParticleDefinition*>  G4ShortLivedList;
#else
   // Use  RWTPtrOrderedVector as list of shortlives
   typedef RWTPtrOrderedVector<G4ParticleDefinition> G4ShortLivedList;
#endif

 public:
   G4ShortLivedTable();

 protected:
   G4ShortLivedTable(const  G4ShortLivedTable &right);

 public:
   virtual ~G4ShortLivedTable();

   G4bool                IsShortLived(G4ParticleDefinition*) const;
   // return true if the particle is shortlived particle
  
   void DumpTable(const G4String &particle_name = "ALL") const;
   // dump information of particles specified by name 

   G4int                 Entries() const;
   G4bool                Contains(const G4ParticleDefinition *particle) const;
   void                  Insert(G4ParticleDefinition* particle);
   void                  Remove(G4ParticleDefinition* particle);
   G4ParticleDefinition* GetParticle(G4int index) const;

 protected:
   G4int                GetVerboseLevel() const;

 private:
   G4ShortLivedList*                  fShortLivedList;

};

inline G4bool  G4ShortLivedTable::Contains(const G4ParticleDefinition* particle) const
{
#ifdef G4USE_STL_MAP
  G4ShortLivedList::iterator i;
  for (i = fShortLivedList->begin(); i!= fShortLivedList->end(); ++i) {
    if (**i==*particle) return true;
  }
  return false;
#else
  return fShortLivedList->contains(particle);
#endif
}

inline G4int G4ShortLivedTable::Entries() const
{
#ifdef G4USE_STL_MAP
  return fShortLivedList->size();
#else
  return fShortLivedList->entries();
#endif
}


inline G4ParticleDefinition*  G4ShortLivedTable::GetParticle(G4int index) const
{
  if ( (index >=0 ) && (index < Entries()) ) {
    return (*fShortLivedList)[index];
  } else {
    return 0; 
  } 
}



#endif
