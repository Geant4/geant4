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
// G4ParticleTable
//
// Class description:
//
// G4ParticleTable is the table of pointers to G4ParticleDefinition.
// It is a "singleton" (only one static object).
// Each G4ParticleDefinition pointer is stored with its name as a key
// to itself. So, each G4ParticleDefinition object must have unique
// name.

// Authors: G.Cosmo, 2 December 1995 - Design, based on object model
//          H.Kurashige, 27 June 1996 - First implementation
// History:
// - 14 Nov 1997, H.Kurashige - Added messenger
// - 24 Sep 1998, H.Kurashige - Added dictionary for encoding
// - 28 Oct 1999, H.Kurashige - Migration to STL maps
// - 15 Sep 2017, K.L.Genser - Added support for MuonicAtom
// --------------------------------------------------------------------
#ifndef G4ParticleTable_hh
#define G4ParticleTable_hh 1

#include "G4ParticleDefinition.hh"
#include "G4ParticleTableIterator.hh"
#include "G4Threading.hh"
#include "G4ios.hh"
#include "globals.hh"

#include <map>

class G4UImessenger;
class G4ParticleMessenger;
class G4IonTable;

class G4ParticleTable
{
  public:
    using G4PTblDictionary = G4ParticleTableIterator<G4String, G4ParticleDefinition*>::Map;
    using G4PTblDicIterator = G4ParticleTableIterator<G4String, G4ParticleDefinition*>;
    using G4PTblEncodingDictionary = G4ParticleTableIterator<G4int, G4ParticleDefinition*>::Map;
    using G4PTblEncodingDicIterator = G4ParticleTableIterator<G4int, G4ParticleDefinition*>;

    virtual ~G4ParticleTable();

    // Copy constructor and assignment operator not allowed
    G4ParticleTable(const G4ParticleTable&) = delete;
    G4ParticleTable& operator=(const G4ParticleTable&) = delete;

    // This method is similar to the constructor. It is used by each worker
    // thread to achieve the partial effect as that of the master thread
    void WorkerG4ParticleTable();

    // This method is similar to the destructor. It is used by each worker
    // thread to achieve the partial effect as that of the master thread
    void DestroyWorkerG4ParticleTable();

    // Return the pointer to the G4ParticleTable object
    // G4ParticleTable is a "singleton" and can get its pointer by this
    // function. At the first time of calling this function, the
    // G4ParticleTable object is instantiated
    static G4ParticleTable* GetParticleTable();

    // Returns TRUE if the ParticleTable contains the particle's pointer
    inline G4bool contains(const G4ParticleDefinition* particle) const;
    G4bool contains(const G4String& particle_name) const;

    // Returns the number of particles in the ParticleTable
    G4int entries() const;
    G4int size() const;

    // Returns a pointer to the i-th particle in the ParticleTable
    // 0 <= index < entries()
    G4ParticleDefinition* GetParticle(G4int index) const;

    // Returns the name of i-th particle in the ParticleTable
    const G4String& GetParticleName(G4int index) const;

    // Returns a pointer to the particle (0 if not contained)
    G4ParticleDefinition* FindParticle(G4int PDGEncoding);
    G4ParticleDefinition* FindParticle(const G4String& particle_name);
    G4ParticleDefinition* FindParticle(const G4ParticleDefinition* particle);

    // Returns a pointer to its anti-particle (0 if not contained)
    inline G4ParticleDefinition* FindAntiParticle(G4int PDGEncoding);
    inline G4ParticleDefinition* FindAntiParticle(const G4String& p_name);
    inline G4ParticleDefinition* FindAntiParticle(const G4ParticleDefinition* p);

    // Returns the pointer to the Iterator
    G4PTblDicIterator* GetIterator() const;

    // Dumps information of particles specified by name
    void DumpTable(const G4String& particle_name = "ALL");

    // Returns the pointer to the G4IonTable object
    G4IonTable* GetIonTable() const;

    // Inserts the particle into ParticleTable.
    // Returned value is the same as particle if successfully inserted
    //   or the pointer to another G4ParticleDefinition object
    //      which has same particle name
    //   or nullptr if failing to insert by other reason
    G4ParticleDefinition* Insert(G4ParticleDefinition* particle);

    // Removes the particle from the table (not delete)
    G4ParticleDefinition* Remove(G4ParticleDefinition* particle);

    // Removes all particles from G4ParticleTable
    void RemoveAllParticles();

    // Removes and deletes all particles from G4ParticleTable
    void DeleteAllParticles();

    // Creates messenger
    G4UImessenger* CreateMessenger();

    void SelectParticle(const G4String& name);

    inline const G4ParticleDefinition* GetSelectedParticle() const;

    inline void SetVerboseLevel(G4int value);
    inline G4int GetVerboseLevel() const;

    inline void SetReadiness(G4bool val = true);
    inline G4bool GetReadiness() const;

    inline G4ParticleDefinition* GetGenericIon() const;
    inline void SetGenericIon(G4ParticleDefinition*);

    inline G4ParticleDefinition* GetGenericMuonicAtom() const;
    inline void SetGenericMuonicAtom(G4ParticleDefinition*);

    // Public data ----------------------------------------------------

    // These fields should be thread local or thread private. For a singleton
    // class, we can change any member field as static without any problem
    // because there is only one instance. Then we are allowed to add
    // "G4ThreadLocal"
    G4ParticleMessenger* fParticleMessenger = nullptr;
    static G4ThreadLocal G4PTblDictionary* fDictionary;
    static G4ThreadLocal G4PTblDicIterator* fIterator;
    static G4ThreadLocal G4PTblEncodingDictionary* fEncodingDictionary;

    // Particle table is being shared
    static G4ParticleTable* fgParticleTable;

    // This field should be thread private. However, we have to keep one copy
    // of the ion table pointer. So we change all important fields of
    // G4IonTable to be thread local
    G4IonTable* fIonTable = nullptr;

    // These shadow pointers are used by each worker thread to copy the content
    // from the master thread
    //
    static G4PTblDictionary* fDictionaryShadow;
    static G4PTblDicIterator* fIteratorShadow;
    static G4PTblEncodingDictionary* fEncodingDictionaryShadow;

#ifdef G4MULTITHREADED
    // Shared instance of a mutex
    static G4GLOB_DLL G4Mutex& particleTableMutex();
    static G4GLOB_DLL G4int& lockCount();
#endif

  protected:
    const G4PTblDictionary* GetDictionary() const;

    // Returns key value of the particle (i.e. particle name)
    inline const G4String& GetKey(const G4ParticleDefinition* particle) const;

    // Returns the pointer to EncodingDictionary
    const G4PTblEncodingDictionary* GetEncodingDictionary() const;

  private:
    // Provate default constructor
    G4ParticleTable();

    void CheckReadiness() const;

    // Private data ---------------------------------------------------

    G4ParticleDefinition* genericIon = nullptr;
    G4ParticleDefinition* genericMuonicAtom = nullptr;
    const G4ParticleDefinition* selectedParticle = nullptr;

    const G4String noName = " ";
    G4String selectedName = "undefined";

    // Control flag for output message
    //  0: Silent
    //  1: Warning message
    //  2: More
    G4int verboseLevel = 1;

    G4bool readyToUse = false;
};

#include "G4ParticleTable.icc"

#endif
