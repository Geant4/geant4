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
// $Id: G4ParticleTable.hh,v 1.22 2010-10-30 07:55:00 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
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
//      modified implementation of Remove 21 Mar.,08  H.Kurashige

#ifndef G4ParticleTable_h
#define G4ParticleTable_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"


#include <map>
#include "G4ParticleTableIterator.hh"

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

   typedef G4ParticleTableIterator<G4String, G4ParticleDefinition*>::Map G4PTblDictionary;
   typedef G4ParticleTableIterator<G4String, G4ParticleDefinition*> G4PTblDicIterator;
   typedef G4ParticleTableIterator<G4int, G4ParticleDefinition*>::Map G4PTblEncodingDictionary;
   typedef G4ParticleTableIterator<G4int, G4ParticleDefinition*> G4PTblEncodingDicIterator;

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
   G4ParticleDefinition* FindIon( G4int    atomicNumber, 
				  G4int    atomicMass, 
				  G4int    numberOfLambda, 
				  G4double excitationEnergy );
   //  return the pointer to an ion  (returns 0 if the ion does not exist)
   //  the ion has excitation energy nearest to given excitationEnergy  (0: ground state)

   G4ParticleDefinition* GetIon(  G4int    atomicNumber, 
				  G4int    atomicMass, 
				  G4double   excitationEnergy);
   G4ParticleDefinition* GetIon(  G4int    atomicNumber, 
				  G4int    atomicMass, 
				  G4int    numberOfLambda, 
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

 public: //With Description

   G4IonTable* GetIonTable();
   // return the pointer to G4IonTable object

   const G4ShortLivedTable* GetShortLivedTable();
   // return the pointer to G4ShortLivedTable object
 
 public: // With Description
   G4ParticleDefinition* Insert(G4ParticleDefinition *particle);
   // insert the particle into ParticleTable 
   // return value is same as particle if successfully inserted
   //              or pointer to another G4ParticleDefinition object 
   //                   which has same name of particle
   //              or 0 if fail to insert by another reason

   G4ParticleDefinition* Remove(G4ParticleDefinition *particle);
   // Remove the particle from the table (not delete)

   void RemoveAllParticles();
   // remove all particles from G4ParticleTable 

   void DeleteAllParticles();
   // remove and delete all particles from G4ParticleTable  

 public:
   G4UImessenger*       CreateMessenger();
   void                 DeleteMessenger();
  // create/delete messenger for the particle table 
 
 protected:
   G4PTblDictionary* GetDictionary();

   const G4String& GetKey(const G4ParticleDefinition *particle) const;
   // return key value of the particle (i.e. particle name)

   const G4PTblEncodingDictionary* GetEncodingDictionary();
   // return the pointer to EncodingDictionary

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

   static G4ParticleTable*  fgParticleTable;

   G4IonTable*            fIonTable;
   G4ShortLivedTable*     fShortLivedTable;

   G4String               noName;

   G4bool  readyToUse;
 
 public:
   void SetReadiness();
   G4bool GetReadiness() const;
 private:
   void CheckReadiness();
};
#include "G4ParticleTable.icc"

#endif






