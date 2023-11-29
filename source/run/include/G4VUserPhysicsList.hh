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
// G4VUserPhysicsList
//
// Class description:
//
// This class is an abstract class for constructing particles and processes.
// User must implement the following two pure virtual methods in the concrete
// class derived from this class:
// - G4VUserPhysicsList::ConstructParticle()
//     Construct particles
// - G4VUserPhysicsList::ConstructProcess()
//     Construct procesess and register them to particles.

// Original author: H.Kurashige (Kobe University), 9 January 1998
// --------------------------------------------------------------------
#ifndef G4VUserPhysicsList_hh
#define G4VUserPhysicsList_hh 1

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicsModelCatalog.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Threading.hh"
#include "G4VUPLSplitter.hh"
#include "G4ios.hh"
#include "globals.hh"

#include "rundefs.hh"
#include "tls.hh"

class G4UserPhysicsListMessenger;
class G4PhysicsListHelper;
class G4VProcess;

// Encapsulate the fields of class G4VUserPhysicsList that are per-thread.
class G4VUPLData
{
  public:
    void initialize();

    G4ParticleTable::G4PTblDicIterator* _theParticleIterator = nullptr;
    G4UserPhysicsListMessenger* _theMessenger = nullptr;
    G4PhysicsListHelper* _thePLHelper = nullptr;
    G4bool _fIsPhysicsTableBuilt = false;
    G4int _fDisplayThreshold = 0;
};

// The type G4VUPLManager is introduced to encapsulate the methods used by
// both the master thread and worker threads to allocate memory space for
// the fields encapsulated by the class G4VUPLData. When each thread
// changes the value for these fields, it refers to them using a macro
// definition defined below. For every G4VUserPhysicsList instance,
// there is a corresponding G4VUPLData instance. All G4VUPLData instances
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
//                 to be associated to a G4RunManager, however user can
//                 instantiate as many PLs are needed and at run-time select one
//                 of the PLs to be used we thus need this mechanism to
//                 guarantee that the system works without problems in case of
//                 this (unusual) case. This may be reviewed in the future
//
using G4VUPLManager = G4VUPLSplitter<G4VUPLData>;
using G4VUserPhysicsListSubInstanceManager = G4VUPLManager;

class G4VUserPhysicsList
{
  public:
    G4VUserPhysicsList();
    virtual ~G4VUserPhysicsList();
    G4VUserPhysicsList(const G4VUserPhysicsList&);
    G4VUserPhysicsList& operator=(const G4VUserPhysicsList&);

    // Each particle type will be instantiated.
    // This method is invoked by the RunManger.
    virtual void ConstructParticle() = 0;

    // By calling the "Construct" method,
    // process manager and processes are created.
    void Construct();

    // Each physics process will be instantiated and
    // registered to the process manager of each particle type.
    // Invoked in the Construct() method.
    virtual void ConstructProcess() = 0;

    // Sets a cut value for all particle types in the particle table.
    virtual void SetCuts();

    // Set/get the default cut value. Calling SetDefaultCutValue() causes
    // re-calcuration of cut values and physics tables just before the
    // next event loop.
    void SetDefaultCutValue(G4double newCutValue);
    G4double GetDefaultCutValue() const;

    // Invoke BuildPhysicsTable for all processes for all particles.
    // In case of "Retrieve" flag is ON, PhysicsTable will be
    // retrieved from files.
    void BuildPhysicsTable();

    // Prepare the PhysicsTable for specified particle type.
    void PreparePhysicsTable(G4ParticleDefinition*);

    // Build the PhysicsTable for specified particle type.
    void BuildPhysicsTable(G4ParticleDefinition*);

    // Store PhysicsTable together with both material and cut value
    // information in files under the specified directory.
    // Returns "true" if files are successfully created.
    G4bool StorePhysicsTable(const G4String& directory = ".");

    // Return true if "Retrieve" flag is ON.
    // (i.e. PhysicsTable will be retrieved from files).
    G4bool IsPhysicsTableRetrieved() const;
    G4bool IsStoredInAscii() const;

    // Get directory path for physics table files.
    const G4String& GetPhysicsTableDirectory() const;

    // Set "Retrieve" flag. Directory path can be set together.
    // Null string (default) means directory is not changed
    // from the current value.
    void SetPhysicsTableRetrieved(const G4String& directory = "");
    void SetStoredInAscii();

    // Reset "Retrieve" flag.
    void ResetPhysicsTableRetrieved();
    void ResetStoredInAscii();

    // Print out the List of registered particles types.
    void DumpList() const;

    // Request to print out information of cut values.
    // Printing will be performed when all tables are made.
    void DumpCutValuesTable(G4int flag = 1);

    // Triggers the print-out requested by the above method.
    // This method must be invoked by RunManager at the proper moment.
    void DumpCutValuesTableIfRequested();

    // Set/get control flag for output message
    //  0: Silent
    //  1: Warning message
    //  2: More
    void SetVerboseLevel(G4int value);
    G4int GetVerboseLevel() const;

    void UseCoupledTransportation(G4bool vl = true);

    // Invokes default SetCuts() method.
    // Note: cut values will not be overwritten.
    // Use of default SetCuts() method is recommended.
    void SetCutsWithDefault();

    // Sets a cut value for a particle type for the default region.
    void SetCutValue(G4double aCut, const G4String& pname);

    // Gets a cut value for a particle type for the default region.
    G4double GetCutValue(const G4String& pname) const;

    // Sets a cut value for a particle type for a region.
    void SetCutValue(G4double aCut, const G4String& pname, const G4String& rname);

    // Invoke SetCuts for specified particle for a region.
    // If the pointer to the region is NULL, the default region is used
    // In case of "Retrieve" flag is ON, cut values will be retrieved
    // from files.
    void SetParticleCuts(G4double cut, G4ParticleDefinition* particle, G4Region* region = nullptr);
    void SetParticleCuts(G4double cut, const G4String& particleName, G4Region* region = nullptr);

    // Invoke SetCuts() for all particles in a region.
    void SetCutsForRegion(G4double aCut, const G4String& rname);

    // Gets/sets the flag for ApplyCuts().
    void SetApplyCuts(G4bool value, const G4String& name);
    G4bool GetApplyCuts(const G4String& name) const;

    // Remove and delete ProcessManagers for all particles in the
    // Particle Table.
    void RemoveProcessManager();

    // Remove and delete TrackingManagers for all particles in the
    // Particle Table.
    void RemoveTrackingManager();

    // Add process manager for particles created on-the-fly.
    void AddProcessManager(G4ParticleDefinition* newParticle,
                           G4ProcessManager* newManager = nullptr);

    // Check consistencies of list of particles.
    void CheckParticleList();

    void DisableCheckParticleList();

    inline G4int GetInstanceID() const;
    static const G4VUPLManager& GetSubInstanceManager();

    // Used by Worker threads on the shared instance of physics-list
    // to initialise workers. Derived class re-implementing this method
    // must also call this base class method.
    virtual void InitializeWorker();

    // Destroy thread-local data. Note that derived classes
    // implementing this method should still call this base class one.
    virtual void TerminateWorker();

  protected:
    // User must invoke this method in his ConstructProcess()
    // implementation in order to enable particle transportation.
    void AddTransportation();

    // Register a process to the particle type
    // according to the ordering parameter table.
    // 'true' is returned if the process is registerd successfully.
    G4bool RegisterProcess(G4VProcess* process, G4ParticleDefinition* particle);

    // Build PhysicsTable for making the integral schema.
    void BuildIntegralPhysicsTable(G4VProcess*, G4ParticleDefinition*);

    // Retrieve PhysicsTable from files for process belonging to the particle.
    // Normal BuildPhysics procedure of processes will be invoked, if it
    // fails (in case of process's RetrievePhysicsTable() returns false).
    virtual void RetrievePhysicsTable(G4ParticleDefinition*, const G4String& directory,
                                      G4bool ascii = false);

    // Adds new ProcessManager to all particles in the Particle Table.
    // This function is used in Construct().
    void InitializeProcessManager();

    G4ParticleTable::G4PTblDicIterator* GetParticleIterator() const;

  protected:
    // The particle table has the complete List of existing particle types.
    G4ParticleTable* theParticleTable = nullptr;

    G4int verboseLevel = 1;

    // Default cut value for all particles
    G4double defaultCutValue = 1.0;
    G4bool isSetDefaultCutValue = false;

    // Pointer to ProductionCutsTable.
    G4ProductionCutsTable* fCutsTable = nullptr;

    // Flag to determine if physics table will be build from file or not.
    G4bool fRetrievePhysicsTable = false;
    G4bool fStoredInAscii = true;

    G4bool fIsCheckedForRetrievePhysicsTable = false;
    G4bool fIsRestoredCutValues = false;

    // Directory name for physics table files.
    G4String directoryPhysicsTable = ".";

    // Flag for CheckParticleList().
    G4bool fDisableCheckParticleList = false;

    // MT data
    G4int g4vuplInstanceID = 0;
    G4RUN_DLL static G4VUPLManager subInstanceManager;
};

// Inline methods implementations

inline void G4VUserPhysicsList::Construct()
{
#ifdef G4VERBOSE
  if (verboseLevel > 1) G4cout << "G4VUserPhysicsList::Construct()" << G4endl;
#endif

  if (G4Threading::IsMasterThread()) G4PhysicsModelCatalog::Initialize();

  InitializeProcessManager();

#ifdef G4VERBOSE
  if (verboseLevel > 1) G4cout << "Construct processes " << G4endl;
#endif
  ConstructProcess();
}

inline G4double G4VUserPhysicsList::GetDefaultCutValue() const
{
  return defaultCutValue;
}

inline G4int G4VUserPhysicsList::GetVerboseLevel() const
{
  return verboseLevel;
}

inline G4bool G4VUserPhysicsList::IsPhysicsTableRetrieved() const
{
  return fRetrievePhysicsTable;
}

inline G4bool G4VUserPhysicsList::IsStoredInAscii() const
{
  return fStoredInAscii;
}

inline const G4String& G4VUserPhysicsList::GetPhysicsTableDirectory() const
{
  return directoryPhysicsTable;
}

inline void G4VUserPhysicsList::SetStoredInAscii()
{
  fStoredInAscii = true;
}

inline void G4VUserPhysicsList::ResetPhysicsTableRetrieved()
{
  fRetrievePhysicsTable = false;
  fIsRestoredCutValues = false;
  fIsCheckedForRetrievePhysicsTable = false;
}

inline void G4VUserPhysicsList::ResetStoredInAscii()
{
  fStoredInAscii = false;
}

inline void G4VUserPhysicsList::DisableCheckParticleList()
{
  fDisableCheckParticleList = true;
}

inline G4int G4VUserPhysicsList::GetInstanceID() const
{
  return g4vuplInstanceID;
}

inline const G4VUPLManager& G4VUserPhysicsList::GetSubInstanceManager()
{
  return subInstanceManager;
}

#endif
