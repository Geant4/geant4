// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParticleTable.hh,v 1.10 1999-10-29 08:03:52 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, IT Division, ASD group
//	History: first implementation, based on object model of
//	27 June 1996, H.Kurashige
// ------------------------------------------------------------
//      added fParticleMessenger         14 Nov., 97 H.Kurashige
//      added Create/DeleteMessenger     06 Jul., 98 H.Kurashige
//      modified FindIon                 02 Aug., 98 H.Kurashige
//      added dictionary for encoding    24 Sep., 98 H.Kurashige
//      added RemoveAllParticles()        8 Nov., 98 H.Kurashige
//         --------------------------------
//      fixed  some improper codings     08 Apr., 99 H.Kurashige
//      modified FindIon/GetIon methods  17 AUg., 99 H.Kurashige
//      implement new version for using STL map instaed of RW PtrHashedDictionary
//                                       28 ct., 99  H.Kurashige

#ifndef G4ParticleTable_h
#define G4ParticleTable_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"

#ifdef G4USE_STL_MAP
#include "g4std/map"
#include "G4ParticleTableIterator.hh"
#else
#include <rw/tphdict.h>
#endif 

class G4UImessenger;
class G4ParticleMessenger;
class G4IonTable;
class G4ShortLivedTable;

class G4ParticleTable
{
 // Class Description
 //   G4ParticleTable is the table of pointer to G4ParticleDefinition
 //   G4ParticleTable is a "singleton" (only one and staic object)
 //   In G4ParticleTable, each G4ParticleDefinition pointer is stored
 //   with its name as a key to itself. So, each  G4ParticleDefinition 
 //   object must have unique name for itself. 
 //   

 public:
#ifdef G4USE_STL_MAP
   typedef G4ParticleTableIterator<G4String, G4ParticleDefinition*>::Map G4PTblDictionary;
   typedef G4ParticleTableIterator<G4String, G4ParticleDefinition*> G4PTblDicIterator;
   typedef G4ParticleTableIterator<G4int, G4ParticleDefinition*>::Map G4PTblEncodingDictionary;
   typedef G4ParticleTableIterator<G4int, G4ParticleDefinition*> G4PTblEncodingDicIterator;
#else
   typedef RWTPtrHashDictionary<G4String,G4ParticleDefinition> G4PTblDictionary;
   typedef RWTPtrHashDictionaryIterator<G4String,G4ParticleDefinition> G4PTblDicIterator;
   typedef RWTPtrHashDictionary<G4int,G4ParticleDefinition> G4PTblEncodingDictionary;
   typedef RWTPtrHashDictionaryIterator<G4int,G4ParticleDefinition> G4PTblEncodingDicIterator;
#endif

 protected:
   G4ParticleTable();
   G4ParticleTable(const  G4ParticleTable &right);

 public:
   virtual ~G4ParticleTable();
  
 public: // With Description
   static G4ParticleTable* GetParticleTable();
   // return the pointer to G4ParticleTable object
   //   G4ParticleTable is a "singleton" and can get its pointer by this function
   //   At the first time of calling this function, the G4ParticleTable object
   //   is instantiated 

   G4bool   contains(const G4ParticleDefinition *particle);
   G4bool   contains(const G4String &particle_name);
   // returns TRUE if the ParticleTable contains 

   G4int    entries() const;
   G4int    size() const;
   // returns the number of Particles in the ParticleTable
   
   G4ParticleDefinition* GetParticle(G4int index);
   // returns a pointer to i-th particles in the ParticleTable
   //    0<= index < entries()

   const G4String& GetParticleName(G4int index);
   // returns name of i-th particles in the ParticleTable

   G4ParticleDefinition* FindParticle(G4int  PDGEncoding );
   G4ParticleDefinition* FindParticle(const G4String &particle_name);
   G4ParticleDefinition* FindParticle(const G4ParticleDefinition *particle);
   // returns a pointer to the particle (0 if not contained)

   G4ParticleDefinition* FindAntiParticle(G4int  PDGEncoding );
   G4ParticleDefinition* FindAntiParticle(const G4String &particle_name);
   G4ParticleDefinition* FindAntiParticle(const G4ParticleDefinition *particle);
   // returns a pointer to its anti-particle (0 if not contained)

   G4ParticleDefinition* FindIon( G4int    atomicNumber, 
				  G4int    atomicMass, 
				  G4double excitationEnergy );
   //  return the pointer to an ion  (returns 0 if the ion does not exist)
   //  the ion has excitation energy nearest to given excitationEnergy  (0: ground state)

   G4ParticleDefinition* GetIon(  G4int    atomicNumber, 
				  G4int    atomicMass, 
				  G4double   excitationEnergy);
   //  return the pointer to an ion ( create ion if the ion does not exist)
   //  It has excitation energy nearest to given excitationEnergy  (0: ground state)
  
   G4ParticleDefinition* FindIon( G4int atomicNumber, 
				  G4int atomicMass, 
				  G4int dummy1,
				  G4int dummy2 );
   //  return the pointer to an ion
   //  !! This routine behaves same as GetIon( atomicNumber, atomicMass, 0) 
   //  !! The third and fourth arguments are meaningless
   //  !! This routine is provided for compatibility to old version

   G4PTblDicIterator* GetIterator();
   // return the pointer of Iterator (RW compatible)
   
   void DumpTable(const G4String &particle_name = "ALL");
   // dump information of particles specified by name 

 public: 
   G4ParticleDefinition* Insert(G4ParticleDefinition *particle);
   // insert the particle into ParticleTable 
   // return value is same as particle if successfully inserted
   //              or pointer to another G4ParticleDefinition object 
   //                   which has same name of particle
   //              or 0 if fail to insert by another reason

 protected:
   G4ParticleDefinition* Remove(G4ParticleDefinition *particle);
   // Remove Particle

   G4PTblDictionary* GetDictionary();

   const G4String& GetKey(const G4ParticleDefinition *particle) const;
   // return key value of the particle (i.e. particle name)

   const G4PTblEncodingDictionary* GetEncodingDictionary();
   // return the pointer to EncodingDictionary

 public: //With Description

   const G4IonTable* GetIonTable();
   // return the pointer to G4IonTable object

   const G4ShortLivedTable* GetShortLivedTable();
   // return the pointer to G4ShortLivedTable object
 
 public:
   G4UImessenger*       CreateMessenger();
   void                 DeleteMessenger();
  // create/delete messenger for the particle table 
 
 protected:  
   void RemoveAllParticles();
   // remove all particles from G4ParticleTable and 
   // delete them if they were created dynamically  (i.e. not static objects) 

#ifndef G4USE_STL_MAP
 protected:
   static unsigned HashFun(const G4String& particle_name);
   static unsigned EncodingHashFun(const G4int& aEndcoding);
  // hash functions  
#endif

 private:
   G4int verboseLevel;
   // controle flag for output message
   //  0: Silent
   //  1: Warning message
   //  2: More

 public:
   void  SetVerboseLevel(G4int value);
   G4int GetVerboseLevel() const;

 private:
   G4ParticleMessenger* fParticleMessenger;
   G4PTblDictionary*  fDictionary;
   G4PTblDicIterator* fIterator;
   G4PTblEncodingDictionary* fEncodingDictionary;
   G4int DictionaryBucketSize;

   static G4ParticleTable*  fgParticleTable;

   G4IonTable*            fIonTable;
   G4ShortLivedTable*     fShortLivedTable;

   G4String               noName;
};
#include "G4ParticleTable.icc"

#endif






