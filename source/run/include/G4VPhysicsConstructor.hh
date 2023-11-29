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
// G4VPhysicsConstructor
//
// Class description:
//
// This class is a virtual class for constructing particles and processes.
// This class objects is being registered to G4VPhysicsList.
//
// User must implement following four virtual methods in the concrete class
// derived from this class:
//
// - virtual void ConstructParticle();
//     All necessary particle type will be instantiated.
// - virtual void ConstructProcess();
//     All physics processes will be instantiated and
//     registered to the process manager of each particle type.
//
// Only one physics constructor can be registered to a Modular Physics List
// for each "physics_type". Physics constructors with same "physics_type"
// can be replaced by using the method:
//   G4VModularPhysicsList::ReplacePhysics().

// Original author: H.Kurashige (Kobe University), 12 November 2000
// --------------------------------------------------------------------
#ifndef G4VPhysicsConstructor_hh
#define G4VPhysicsConstructor_hh 1

#include "G4ParticleTable.hh"
#include "G4PhysicsListHelper.hh"
#include "G4VUPLSplitter.hh"
#include "G4ios.hh"
#include "globals.hh"

#include "rundefs.hh"

#include <vector>

class G4PhysicsBuilderInterface;

class G4VPCData
{
    // Encapsulate the fields of class G4VPhysicsConstructor
    // that are per-thread.

  public:
    using PhysicsBuilders_V = std::vector<G4PhysicsBuilderInterface*>;
    void initialize();
    G4ParticleTable::G4PTblDicIterator* _aParticleIterator;

    PhysicsBuilders_V* _builders = nullptr;
};

// The type G4VPCManager is introduced to encapsulate the methods used by
// both the master thread and worker threads to allocate memory space for
// the fields encapsulated by the class G4VPCData. When each thread
// changes the value for these fields, it refers to them using a macro
// definition defined below. For every G4VPhysicsConstructor instance,
// there is a corresponding G4VPCData instance. All G4VPCData instances
// are organized by the class G4VUPLManager as an array.
// The field "int g4vuplInstanceID" is added to the class G4VUserPhysicsList.
// The value of this field in each G4VUserPhysicsList instance is the
// subscript of the corresponding G44VUPLData instance.
// In order to use the class G44VUPLManager, we add a static member in the class
// G4VUserPhysicsList as follows: "static G4VUPLManager subInstanceManager".
// Both the master thread and worker threads change the length of the array
// for G44VUPLData instances mutually along with G4VUserPhysicsList
// instances are created. For each worker thread, it dynamically creates ions.
// Consider any thread A, if there is any other thread which creates an ion.
// This ion is shared by the thread A. So the thread A leaves an empty space
// in the array of G4PDefData instances for the ion.
//
// Important Note: you may wonder why we are introducing this mechanism
//                 since there is only one PL for each application.
//                 This is true, in the sense that only one PL is allowed
//                 to be associated to a G4RunManager, however a user can
//                 instantiate as many PLs are needed and at run-time select
//                 one of the PLs to be used we thus need this mechanism to
//                 guarantee that the system works without problems in case of
//                 this (unusual) case. This may be reviewed in the future
//
using G4VPCManager = G4VUPLSplitter<G4VPCData>;
using G4VPhyscicsConstructorManager = G4VPCManager;

class G4VPhysicsConstructor
{
  public:
    G4VPhysicsConstructor(const G4String& = "");
    G4VPhysicsConstructor(const G4String& name, G4int physics_type);
    virtual ~G4VPhysicsConstructor();

    // This method will be invoked in the Construct() method.
    // Each particle type will be instantiated.
    virtual void ConstructParticle() = 0;

    // This method will be invoked in the Construct() method.
    // Each physics process will be instantiated and
    // registered to the process manager of each particle type.
    virtual void ConstructProcess() = 0;

    inline void SetPhysicsName(const G4String& = "");
    inline const G4String& GetPhysicsName() const;

    inline void SetPhysicsType(G4int);
    inline G4int GetPhysicsType() const;

    inline G4int GetInstanceID() const;
    static const G4VPCManager& GetSubInstanceManager();

    // Method called by kernel to destroy thread-local data, equivalent to
    // destructor in sequential mode. Derived classes implementing this
    // method, must also call this base class method.
    virtual void TerminateWorker();

    // Set/get control flag for output message
    //  0: Silent
    //  1: Warning message
    //  2: More
    // verbose level is set equal to physics list when registered.
    inline void SetVerboseLevel(G4int value);
    inline G4int GetVerboseLevel() const;

  protected:
    using PhysicsBuilder_V = G4VPCData::PhysicsBuilders_V;

    // Register a process to the particle type according to the ordering
    // parameter table. 'true' is returned if the process is registered
    // successfully.
    inline G4bool RegisterProcess(G4VProcess* process, G4ParticleDefinition* particle);

    G4ParticleTable::G4PTblDicIterator* GetParticleIterator() const;

    // This returns a copy of the vector of pointers.
    PhysicsBuilder_V GetBuilders() const;

    void AddBuilder(G4PhysicsBuilderInterface* bld);

  protected:
    G4int verboseLevel = 0;
    G4String namePhysics = "";
    G4int typePhysics = 0;

    G4ParticleTable* theParticleTable = nullptr;
    G4int g4vpcInstanceID = 0;
    G4RUN_DLL static G4VPCManager subInstanceManager;
};

// Inline methods implementations

inline void G4VPhysicsConstructor::SetVerboseLevel(G4int value)
{
  verboseLevel = value;
}

inline G4int G4VPhysicsConstructor::GetVerboseLevel() const
{
  return verboseLevel;
}

inline void G4VPhysicsConstructor::SetPhysicsName(const G4String& name)
{
  namePhysics = name;
}

inline const G4String& G4VPhysicsConstructor::GetPhysicsName() const
{
  return namePhysics;
}

inline void G4VPhysicsConstructor::SetPhysicsType(G4int val)
{
  if (val > 0) {
    typePhysics = val;
  }
}

inline G4int G4VPhysicsConstructor::GetPhysicsType() const
{
  return typePhysics;
}

inline G4bool G4VPhysicsConstructor::RegisterProcess(G4VProcess* process,
                                                     G4ParticleDefinition* particle)
{
  return G4PhysicsListHelper::GetPhysicsListHelper()->RegisterProcess(process, particle);
}

inline const G4VPCManager& G4VPhysicsConstructor::GetSubInstanceManager()
{
  return subInstanceManager;
}

#endif
