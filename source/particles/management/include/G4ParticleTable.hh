// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParticleTable.hh,v 1.1 1999-01-07 16:10:30 gunter Exp $
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


#ifndef G4ParticleTable_h
#define G4ParticleTable_h 1

#include "G4ios.hh"
#include <rw/tphdict.h>
#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include <rw/cstring.h>
class G4UImessenger;
class G4ParticleMessenger;
class G4IonTable;
class G4ShortLivedTable;

class G4ParticleTable
{
 //   G4ParticleTable is the table of pointer to G4ParticleDefinition
 //   G4ParticleTable is a "singleton" (only one and staic object)
 //   In G4ParticleTable, each G4ParticleDefinition pointer is stored
 //   with its name as a key to itself. So, each  G4ParticleDefinition 
 //   object must have unique name for itself. 
 //   

 public:
   typedef RWTPtrHashDictionary<G4String,G4ParticleDefinition> G4PTblDictionary;
   typedef RWTPtrHashDictionaryIterator<G4String,G4ParticleDefinition> G4PTblDicIterator;
   typedef RWTPtrHashDictionary<G4int,G4ParticleDefinition> G4PTblEncodingDictionary;
   typedef RWTPtrHashDictionaryIterator<G4int,G4ParticleDefinition> G4PTblEncodingDicIterator;

 protected:
   G4ParticleTable();
   G4ParticleTable(const  G4ParticleTable &right);

 public:
   virtual ~G4ParticleTable();

   void RemoveAllParticles();
   // remove all particles from G4ParticleTable and 
   // delete them if they were created dynamically  (i.e. not static objects) 

   static G4ParticleTable* GetParticleTable();
   // return the pointer to G4ParticleTable object
   //   G4ParticleTable is a "singleton" and can get its pointer by this function
   //   At the first time of calling this function, the G4ParticleTable object
   //   is instantiated 

   G4bool   contains(const G4ParticleDefinition *particle) const;
   G4bool   contains(const G4String &particle_name) const;
   // returns TRUE if the ParticleTable contains 

   G4int    entries() const;
   // returns the number of Particles in the ParticleTable
   
   G4ParticleDefinition* GetParticle(G4int index);
   // returns a pointer to i-th particles in the ParticleTable
   //    0<= index < entries()

   G4String GetParticleName(G4int index);
   // returns name of i-th particles in the ParticleTable

   G4ParticleDefinition* FindParticle(G4int  PDGEncoding )  const;
   G4ParticleDefinition* FindParticle(const G4String &particle_name)  const;
   G4ParticleDefinition* FindParticle(const G4ParticleDefinition *particle) const;
   // returns a pointer to the particle (NULL if not contained)

   G4ParticleDefinition* FindAntiParticle(G4int  PDGEncoding )  const;
   G4ParticleDefinition* FindAntiParticle(const G4String &particle_name) const;
   G4ParticleDefinition* FindAntiParticle(const G4ParticleDefinition *particle) const;
   // returns a pointer to its anti-particle (NULL if not contained)

   G4ParticleDefinition* FindIon( G4int atomicNumber, 
				  G4int atomicMass, 
				  G4int iAngularMomentum = 0, 
				  G4int iCharge = -1);
   //  return the pointer to an ion
   //    iAngularMomentum is integer and should be given in unit of 1/2 
   //    (i.e. you should set iAngularMomentum = 2 if you want an ion with total angular momentum = 1)
   //    iCharge is initeger and should be given in unit of eplus
   //    you can omit spin (default value = 0) and charge (default value = -1 : fully ionized)

   G4ParticleDefinition* Insert(G4ParticleDefinition *particle);
   // insert the particle into ParticleTable 
   // return value is same as particle if successfully inserted
   //              or pointer to another G4ParticleDefinition object which has same name of particle
   //              or NULL if fail to insert by another reason

   G4ParticleDefinition* Remove(G4ParticleDefinition *particle);
   // Remove Particle

   void DumpTable(const G4String &particle_name = "ALL") const;
   // dump information of particles specified by name 

 public:
   G4String GetKey(const G4ParticleDefinition *particle) const;
   G4PTblDictionary* GetDictionary() const;
   G4PTblDicIterator* GetIterator() const; 
   G4PTblEncodingDictionary* GetEncodingDictionary() const;

   G4IonTable* GetIonTable() const;
   // return the pointer to G4IonTable object

   G4ShortLivedTable* GetShortLivedTable() const;
   // return the pointer to G4ShortLivedTable object
 
 public:
   G4UImessenger* CreateMessenger();
   void           DeleteMessenger();

 protected:
   static unsigned HashFun(const G4String& particle_name);
   static unsigned EncodingHashFun(const G4int& aEndcoding);
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
};

inline G4ShortLivedTable*  G4ParticleTable::GetShortLivedTable() const 
{
  return fShortLivedTable;
}

inline G4IonTable*  G4ParticleTable::GetIonTable() const 
{
  return fIonTable;
}

inline  void G4ParticleTable::SetVerboseLevel(G4int value )
{ 
  verboseLevel = value; 
}
inline G4int G4ParticleTable::GetVerboseLevel() const 
{ 
  return verboseLevel; 
}

inline G4ParticleTable::G4PTblDictionary* G4ParticleTable::GetDictionary() const
{
  return fDictionary;
}

inline G4ParticleTable::G4PTblDicIterator* G4ParticleTable::GetIterator() const
{
  return fIterator;
}

inline G4ParticleTable::G4PTblEncodingDictionary* G4ParticleTable::GetEncodingDictionary() const
{
  return fEncodingDictionary;
}

inline G4String G4ParticleTable::GetKey(const G4ParticleDefinition *particle) const
{
  return particle->GetParticleName();
}

inline G4ParticleDefinition* G4ParticleTable::FindParticle(const G4String &particle_name) const
{
  return GetDictionary() -> findValue(&particle_name);
}

inline G4ParticleDefinition* G4ParticleTable::FindParticle(const G4ParticleDefinition *particle) const
{
  G4String key = GetKey(particle);
  return GetDictionary() -> findValue( &key );
}

inline G4ParticleDefinition* G4ParticleTable::FindAntiParticle(G4int aPDGEncoding) const
{
  return FindParticle( FindParticle(aPDGEncoding)->GetAntiPDGEncoding() );
}

inline G4ParticleDefinition* G4ParticleTable::FindAntiParticle(const G4String &particle_name) const
{
  G4int pcode = FindParticle(particle_name) -> GetAntiPDGEncoding();
  return FindParticle(pcode);
}

inline G4ParticleDefinition* G4ParticleTable::FindAntiParticle(const G4ParticleDefinition *particle) const
{
  G4int pcode = particle -> GetAntiPDGEncoding();
  return FindParticle(pcode);
}

inline G4bool  G4ParticleTable::contains(const G4String &particle_name) const
{
   return GetDictionary() -> contains(&particle_name);
}

inline G4bool  G4ParticleTable::contains(const G4ParticleDefinition *particle) const
{
 G4String key = GetKey(particle);
 return GetDictionary() -> contains(&key);
}

inline G4String G4ParticleTable::GetParticleName(G4int index)
{
  G4ParticleDefinition* aParticle =GetParticle(index);
  G4String name = "";
  if (aParticle != NULL) name = aParticle->GetParticleName();
  return name;
}

inline G4int G4ParticleTable::entries() const
{
  return GetDictionary() -> entries();
}


inline unsigned G4ParticleTable::HashFun(const G4String& particle_name)
{
  return particle_name.hash(); 
}

inline unsigned G4ParticleTable::EncodingHashFun(const G4int& aEncoding)
{
  if (aEncoding >0 ) return aEncoding*2;
  else return abs(aEncoding)*2+1;
}
#endif






