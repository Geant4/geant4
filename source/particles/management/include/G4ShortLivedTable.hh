// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ShortLivedTable.hh,v 1.1 1999-01-07 16:10:31 gunter Exp $
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
#include <rw/tpordvec.h>
#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include <rw/cstring.h>

class G4ParticleTable;

class G4ShortLivedTable
{
 //   G4ShortLivedTable is the table of pointer to G4ParticleDefinition
 //   In G4ShortLivedTable, each G4ParticleDefinition pointer is stored

 public:
   typedef RWTPtrOrderedVector<G4ParticleDefinition> G4ShortLivedList;

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
  return fShortLivedList->contains(particle);
}

inline G4int G4ShortLivedTable::Entries() const
{
  return fShortLivedList->entries();
}


inline G4ParticleDefinition*  G4ShortLivedTable::GetParticle(G4int index) const
{
  return (*fShortLivedList)[index];
}



#endif
