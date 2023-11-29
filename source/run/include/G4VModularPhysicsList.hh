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
// G4VModularPhysicsList
//
// Class description:
//
// This class is a subclass of G4VUserPhysicsList.
// The user should register his/her physics constructors by using:
//   G4VModularPhysicsList::RegsiterPhysics()
// to construt particles and processes.
//
// Only one physics constructor can be registered for each "physics_type".
// Physics constructors with same "physics_type" can be replaced using the
// G4VModularPhysicsList::ReplacePhysics() method.

// Original author: H.Kurashige (Kobe University), 12 November 2000
// --------------------------------------------------------------------
#ifndef G4VModularPhysicsList_hh
#define G4VModularPhysicsList_hh 1

#include "G4VPhysicsConstructor.hh"
#include "G4VUPLSplitter.hh"
#include "G4VUserPhysicsList.hh"
#include "G4ios.hh"
#include "globals.hh"

#include "rundefs.hh"

#include <vector>

class G4VMPLData
{
    // Encapsulate the fields of class G4VModularPhysicsList
    // that are per-thread.

  public:
    void initialize();
    using G4PhysConstVectorData = std::vector<G4VPhysicsConstructor*>;
    // See: https://jira-geant4.kek.jp/browse/DEV-284
    G4PhysConstVectorData* physicsVector = nullptr;
};

// The type G4VMPLManager is introduced to encapsulate the methods used by
// both the master thread and worker threads to allocate memory space for
// the fields encapsulated by the class G4VMPLData. When each thread
// changes the value for these fields, it refers to them using a macro
// definition defined below. For every G4VUserPhysicsList instance,
// there is a corresponding G4VMPLData instance. All G4VMPLData instances
// are organized by the class G4VMPLManager as an array.
// The field "int G4VMPLInstanceID" is added to the class G4VUserPhysicsList.
// The value of this field in each G4VUserPhysicsList instance is the
// subscript of the corresponding G44VUPLData instance.
// In order to use the class G44VUPLManager, we add a static member in the class
// G4VUserPhysicsList as follows: "static G4VMPLManager subInstanceManager".
// Both the master thread and worker threads change the length of the array
// for G44VUPLData instances mutually along with G4VUserPhysicsList
// instances are created.
//
using G4VMPLManager = G4VUPLSplitter<G4VMPLData>;
using G4VModularPhysicsListSubInstanceManager = G4VMPLManager;

class G4VModularPhysicsList : public virtual G4VUserPhysicsList
{
  public:
    G4VModularPhysicsList();
    ~G4VModularPhysicsList() override;

    // This method will be invoked in the Construct() method.
    // Each particle type will be instantiated.
    void ConstructParticle() override;

    // This method will be invoked in the Construct() method.
    // Each physics process will be instantiated and
    // registered to the process manager of each particle type.
    void ConstructProcess() override;

    // Register Physics Constructor.
    void RegisterPhysics(G4VPhysicsConstructor*);

    const G4VPhysicsConstructor* GetPhysics(G4int index) const;
    const G4VPhysicsConstructor* GetPhysics(const G4String& name) const;
    const G4VPhysicsConstructor* GetPhysicsWithType(G4int physics_type) const;

    // Replace the Physics Constructor.
    // The existing physics constructor with same physics_type as one of
    // the given physics constructor is replaced (existing physics will be
    // deleted). If a corresponding physics constructor is NOT found,
    // the given physics constructor is just added.
    void ReplacePhysics(G4VPhysicsConstructor*);

    // Remove the Physics Constructor from the list.
    void RemovePhysics(G4VPhysicsConstructor*);
    void RemovePhysics(G4int type);
    void RemovePhysics(const G4String& name);

    inline G4int GetInstanceID() const;
    static const G4VMPLManager& GetSubInstanceManager();
    void TerminateWorker() override;

    // Set/get control flag for output message
    //  0: Silent
    //  1: Warning message
    //  2: More
    // given verbose level is set to all physics constructors.
    void SetVerboseLevel(G4int value);
    G4int GetVerboseLevel() const;

  protected:
    // Protected copy constructor and assignment operator.
    G4VModularPhysicsList(const G4VModularPhysicsList&);
    G4VModularPhysicsList& operator=(const G4VModularPhysicsList&);

    using G4PhysConstVector = G4VMPLData::G4PhysConstVectorData;

    G4int verboseLevel = 0;
    G4int g4vmplInstanceID = 0;
    G4RUN_DLL static G4VMPLManager G4VMPLsubInstanceManager;
};

// Inline methods implementations

inline G4int G4VModularPhysicsList::GetVerboseLevel() const
{
  return verboseLevel;
}

inline G4int G4VModularPhysicsList::GetInstanceID() const
{
  return g4vmplInstanceID;
}

inline const G4VMPLManager& G4VModularPhysicsList::GetSubInstanceManager()
{
  return G4VMPLsubInstanceManager;
}

#endif
