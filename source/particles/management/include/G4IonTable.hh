// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IonTable.hh,v 1.6 1999-08-18 09:15:12 kurasige Exp $
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
//      add GetNucleusMass              15 Mar. 99  H.Kurashige
//          -----
//      Modified GetIon methods         17 Aug. 99 H.Kurashige

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
   // hide copy construictor as protected 
   G4IonTable(const  G4IonTable &right);

 public:
  // destructor
   virtual ~G4IonTable();


 public: 
   G4ParticleDefinition* FindIon(G4int Z, G4int A, G4int L);
   // return pointer of ion if it exists 
   G4ParticleDefinition* GetIon(G4int Z, G4int A, G4int L);
   // get a pointer to the ion with A,Z,J 
   // The ion will be created if not exist yet
   //   Z: Atomic Number
   //   A: Atomic Mass
   //   L: Excitation level (0=stable state)

   G4ParticleDefinition* FindIon(G4int Z, G4int A, G4double E);
   // return pointer of ion if it exists 
   G4ParticleDefinition* GetIon(G4int Z, G4int A, G4double E);
   // get a pointer to the ion with A,Z,J 
   // The ion will be created if not exist yet
   //   Z: Atomic Number
   //   A: Atomic Mass
   //   E: Excitaion energy
   // !!  This method is a dummy now !!
   // !! Implementation will be later !!
 

   G4ParticleDefinition* GetIon(G4int Z, G4int A, G4int J, G4int Q);
   // This method is provided for compatibilties 
   // The last argument of "Q" gives no effect
   // The third argument of J corresponds the excitaion level
   
   // !! PDGCharge inG4ParticleDefinition of ions is           !!
   // !! electric charge of nucleus (i.e. fully ionized ions)  !!

   static G4bool        IsIon(G4ParticleDefinition*);
   // return true if the particle is ion

   G4String             GetIonName(G4int Z, G4int A, G4int L) const;
   // get ion name

   G4double             GetIonMass(G4int Z, G4int A) const;
   G4double             GetNucleusMass(G4int Z, G4int A) const;
   // These two methods returns Nucleus (i.e. full ionized atom) mass 
   // ,where Z is Atomic Number (number of protons) and
   //  A is Atomic Number (number of nucleons)


   G4int                 Entries() const;
   G4bool                Contains(const G4ParticleDefinition *particle) const;
   void                  Insert(G4ParticleDefinition* particle);
   void                  Remove(G4ParticleDefinition* particle);
   G4ParticleDefinition* GetParticle(G4int index) const;

    void DumpTable(const G4String &particle_name = "ALL") const;
   // dump information of particles specified by name 


 protected:

   G4ParticleDefinition* GetLightIon(G4int Z, G4int A) const;
   G4bool                IsLightIon(G4ParticleDefinition*) const;
   // return true if the particle is pre-defined ion
 
   // Utilities
   G4double             GetProtonMass() const;
   G4double             GetNeutronMass() const;
   G4double             GetElectronMass() const;

   // 
   G4int                GetVerboseLevel() const;

 private:
   G4IonList*                  fIonList; 

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










