// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IonTable.hh,v 1.1 1999-01-07 16:10:28 gunter Exp $
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
//      modified GetIon                 02 Aug., 98 H.Kurashige
//      added Remove()                  06 Nov.,98 H.Kurashige

#ifndef G4IonTable_h
#define G4IonTable_h 1

#include "G4ios.hh"
#include <rw/tpordvec.h>
#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include <rw/cstring.h>

class G4ParticleTable;

class G4IonTable
{
 //   G4IonTable is the table of pointer to G4ParticleDefinition
 //   In G4IonTable, each G4ParticleDefinition pointer is stored

 public:
   typedef RWTPtrOrderedVector<G4ParticleDefinition> G4IonList;

 public:
   G4IonTable();

 protected:
   G4IonTable(const  G4IonTable &right);

 public:
   virtual ~G4IonTable();

   G4ParticleDefinition* GetIon(G4int Z, G4int A, G4int J, G4int Q);
   // get a pointer to the ion with A,Z,J,Q 
   // The ion will be created if not exist yet
   
   G4bool                IsIon(G4ParticleDefinition*) const;
   // return true if the particle is ion
  
   void DumpTable(const G4String &particle_name = "ALL") const;
   // dump information of particles specified by name 

   G4String             GetIonName(G4int Z, G4int A, G4int J, G4int Q) const;
   // get ion name
   G4double             GetIonMass(G4int Z, G4int A) const;

   G4int                 Entries() const;
   G4bool                Contains(const G4ParticleDefinition *particle) const;
   void                  Insert(G4ParticleDefinition* particle);
   void                  Remove(G4ParticleDefinition* particle);
   G4ParticleDefinition* GetParticle(G4int index) const;


 protected:
   G4int                GetVerboseLevel() const;

 private:
   G4IonList*                  fIonList;
   G4double                    electronMass, protonMass, neutronMass;
      
   enum { numberOfElements = 110};
   static const G4String       elementName[numberOfElements];

};

inline G4bool  G4IonTable::Contains(const G4ParticleDefinition* particle) const
{
  return fIonList->contains(particle);
}

inline G4int G4IonTable::Entries() const
{
  return fIonList->entries();
}


inline G4ParticleDefinition*  G4IonTable::GetParticle(G4int index) const
{
  return (*fIonList)[index];
}

#endif










