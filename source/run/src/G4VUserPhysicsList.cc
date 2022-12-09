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
// G4VUserPhysicsList implementation
//
// Original author: H.Kurashige (Kobe University), 9 January 1998
// --------------------------------------------------------------------

#include "G4PhysicsListHelper.hh"
#include "G4VUserPhysicsList.hh"

#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ProductionCuts.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4UserPhysicsListMessenger.hh"
#include "G4VTrackingManager.hh"
#include "G4ios.hh"
#include "globals.hh"

#include <fstream>
#include <iomanip>
#include <unordered_set>

// This static member is thread local. For each thread, it holds the array
// size of G4VUPLData instances.
//
#define G4MT_theMessenger                                                      \
  ((this->subInstanceManager.offset[this->g4vuplInstanceID])._theMessenger)
#define G4MT_thePLHelper                                                       \
  ((this->subInstanceManager.offset[this->g4vuplInstanceID])._thePLHelper)
#define fIsPhysicsTableBuilt                                                   \
  ((this->subInstanceManager.offset[this->g4vuplInstanceID])                   \
     ._fIsPhysicsTableBuilt)
#define fDisplayThreshold                                                      \
  ((this->subInstanceManager.offset[this->g4vuplInstanceID])._fDisplayThreshold)
#define theParticleIterator                                                    \
  ((this->subInstanceManager.offset[this->g4vuplInstanceID])                   \
     ._theParticleIterator)

// This field helps to use the class G4VUPLManager
//
G4VUPLManager G4VUserPhysicsList::subInstanceManager;

// --------------------------------------------------------------------
void G4VUPLData::initialize()
{
  _theParticleIterator  = G4ParticleTable::GetParticleTable()->GetIterator();
  _theMessenger         = nullptr;
  _thePLHelper          = G4PhysicsListHelper::GetPhysicsListHelper();
  _fIsPhysicsTableBuilt = false;
  _fDisplayThreshold    = 0;
}

// --------------------------------------------------------------------
G4VUserPhysicsList::G4VUserPhysicsList()
{
  g4vuplInstanceID = subInstanceManager.CreateSubInstance();  // AND
  // default cut value  (1.0mm)
  defaultCutValue = 1.0 * mm;

  // pointer to the particle table
  theParticleTable = G4ParticleTable::GetParticleTable();
  // theParticleIterator = theParticleTable->GetIterator();

  // pointer to the cuts table
  fCutsTable = G4ProductionCutsTable::GetProductionCutsTable();

  // set energy range for SetCut calcuration
  fCutsTable->SetEnergyRange(0.99 * keV, 100 * TeV);

  // UI Messenger
  // theMessenger = new G4UserPhysicsListMessenger(this);
  G4MT_theMessenger = new G4UserPhysicsListMessenger(this);  // AND

  // PhysicsListHelper
  // thePLHelper = G4PhysicsListHelper::GetPhysicsListHelper();
  // thePLHelper->SetVerboseLevel(verboseLevel);
  // G4MT_thePLHelper = G4PhysicsListHelper::GetPhysicsListHelper(); //AND
  G4MT_thePLHelper->SetVerboseLevel(verboseLevel);  // AND

  fIsPhysicsTableBuilt = false;
  fDisplayThreshold    = 0;
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::InitializeWorker()
{
  // Remember messengers are per-thread, so this needs to be done by each
  // worker and due to the presence of "this" cannot be done in
  // G4VUPLData::initialize()
  G4MT_theMessenger = new G4UserPhysicsListMessenger(this);
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::TerminateWorker()
{
  RemoveProcessManager();
  RemoveTrackingManager();
  delete G4MT_theMessenger;
  G4MT_theMessenger = nullptr;
}

// --------------------------------------------------------------------
G4VUserPhysicsList::~G4VUserPhysicsList()
{
  delete G4MT_theMessenger;
  G4MT_theMessenger = nullptr;

  RemoveProcessManager();
  RemoveTrackingManager();

  // invoke DeleteAllParticle
  theParticleTable->DeleteAllParticles();
}

// --------------------------------------------------------------------
G4VUserPhysicsList::G4VUserPhysicsList(const G4VUserPhysicsList& right)
  : verboseLevel(right.verboseLevel)
  , defaultCutValue(right.defaultCutValue)
  , isSetDefaultCutValue(right.isSetDefaultCutValue)
  , fRetrievePhysicsTable(right.fRetrievePhysicsTable)
  , fStoredInAscii(right.fStoredInAscii)
  , fIsCheckedForRetrievePhysicsTable(right.fIsCheckedForRetrievePhysicsTable)
  , fIsRestoredCutValues(right.fIsRestoredCutValues)
  , directoryPhysicsTable(right.directoryPhysicsTable)
  , fDisableCheckParticleList(right.fDisableCheckParticleList)
{
  g4vuplInstanceID = subInstanceManager.CreateSubInstance();
  // pointer to the particle table
  theParticleTable    = G4ParticleTable::GetParticleTable();
  theParticleIterator = theParticleTable->GetIterator();
  // pointer to the cuts table
  fCutsTable = G4ProductionCutsTable::GetProductionCutsTable();

  // UI Messenger
  G4MT_theMessenger = new G4UserPhysicsListMessenger(this);

  // PhysicsListHelper
  G4MT_thePLHelper = G4PhysicsListHelper::GetPhysicsListHelper();
  G4MT_thePLHelper->SetVerboseLevel(verboseLevel);

  fIsPhysicsTableBuilt = right.GetSubInstanceManager()
                           .offset[right.GetInstanceID()]
                           ._fIsPhysicsTableBuilt;
  fDisplayThreshold = right.GetSubInstanceManager()
                        .offset[right.GetInstanceID()]
                        ._fDisplayThreshold;
}

// --------------------------------------------------------------------
G4VUserPhysicsList&
G4VUserPhysicsList::operator=(const G4VUserPhysicsList& right)
{
  if(this != &right)
  {
    verboseLevel                      = right.verboseLevel;
    defaultCutValue                   = right.defaultCutValue;
    isSetDefaultCutValue              = right.isSetDefaultCutValue;
    fRetrievePhysicsTable             = right.fRetrievePhysicsTable;
    fStoredInAscii                    = right.fStoredInAscii;
    fIsCheckedForRetrievePhysicsTable = right.fIsCheckedForRetrievePhysicsTable;
    fIsRestoredCutValues              = right.fIsRestoredCutValues;
    directoryPhysicsTable             = right.directoryPhysicsTable;
    fIsPhysicsTableBuilt = right.GetSubInstanceManager()
                             .offset[right.GetInstanceID()]
                             ._fIsPhysicsTableBuilt;
    fDisplayThreshold = right.GetSubInstanceManager()
                          .offset[right.GetInstanceID()]
                          ._fDisplayThreshold;
    fDisableCheckParticleList = right.fDisableCheckParticleList;
  }
  return *this;
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::AddProcessManager(G4ParticleDefinition* newParticle,
                                           G4ProcessManager*)
{
  if(newParticle == nullptr)
    return;
  G4Exception("G4VUserPhysicsList::AddProcessManager", "Run0252",
              JustWarning, "This method is obsolete");
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::InitializeProcessManager()
{
  // Request lock for particle table access. Some changes are inside
  // this critical region.
#ifdef G4MULTITHREADED
  G4MUTEXLOCK(&G4ParticleTable::particleTableMutex());
  G4ParticleTable::lockCount()++;
#endif
  G4ParticleDefinition* gion =
    G4ParticleTable::GetParticleTable()->GetGenericIon();

  // loop over all particles in G4ParticleTable
  theParticleIterator->reset();
  while((*theParticleIterator)())
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager     = particle->GetProcessManager();

    if(pmanager == nullptr)
    {
      // create process manager if the particle does not have its own.
      pmanager = new G4ProcessManager(particle);
      particle->SetProcessManager(pmanager);
      if(particle->GetMasterProcessManager() == nullptr)
        particle->SetMasterProcessManager(pmanager);
#ifdef G4VERBOSE
      if(verboseLevel > 2)
      {
        G4cout << "G4VUserPhysicsList::InitializeProcessManager: creating "
                  "ProcessManager to "
               << particle->GetParticleName() << G4endl;
      }
#endif
    }
  }

  if(gion != nullptr)
  {
    G4ProcessManager* gionPM = gion->GetProcessManager();
    // loop over all particles once again (this time, with all general ions)
    theParticleIterator->reset(false);
    while((*theParticleIterator)())
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      if(particle->IsGeneralIon())
      {
        particle->SetProcessManager(gionPM);
#ifdef G4VERBOSE
        if(verboseLevel > 2)
        {
          G4cout << "G4VUserPhysicsList::InitializeProcessManager: copying "
                    "ProcessManager to "
                 << particle->GetParticleName() << G4endl;
        }
#endif
      }
    }
  }

  // release lock for particle table access
#ifdef G4MULTITHREADED
  G4MUTEXUNLOCK(&G4ParticleTable::particleTableMutex());
#endif
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::RemoveProcessManager()
{
  // Request lock for particle table access. Some changes are inside
  // this critical region.
#ifdef G4MULTITHREADED
  G4MUTEXLOCK(&G4ParticleTable::particleTableMutex());
  G4ParticleTable::lockCount()++;
#endif

  // loop over all particles in G4ParticleTable
  theParticleIterator->reset();
  while((*theParticleIterator)())
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    if(particle->GetInstanceID() <
       G4ParticleDefinitionSubInstanceManager::slavetotalspace())
    {
      if(particle->GetParticleSubType() != "generic" ||
         particle->GetParticleName() == "GenericIon")
      {
        G4ProcessManager* pmanager = particle->GetProcessManager();
        if(pmanager != nullptr)
          delete pmanager;
#ifdef G4VERBOSE
        if(verboseLevel > 2)
        {
          G4cout << "G4VUserPhysicsList::RemoveProcessManager: ";
          G4cout << "remove ProcessManager from ";
          G4cout << particle->GetParticleName() << G4endl;
        }
#endif
      }
      particle->SetProcessManager(nullptr);
    }
  }

  // release lock for particle table access
#ifdef G4MULTITHREADED
  G4MUTEXUNLOCK(&G4ParticleTable::particleTableMutex());
#endif
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::RemoveTrackingManager()
{
  // One tracking manager may be registered for multiple particles, make sure
  // to delete every object only once.
  std::unordered_set<G4VTrackingManager *> trackingManagers;

  // loop over all particles in G4ParticleTable
  theParticleIterator->reset();
  while((*theParticleIterator)())
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    if (auto *trackingManager = particle->GetTrackingManager())
    {
#ifdef G4VERBOSE
      if(verboseLevel > 2)
      {
        G4cout << "G4VUserPhysicsList::RemoveTrackingManager: ";
        G4cout << "remove TrackingManager from ";
        G4cout << particle->GetParticleName() << G4endl;
      }
#endif
      trackingManagers.insert(trackingManager);
      particle->SetTrackingManager(nullptr);
    }
  }

  for (G4VTrackingManager *tm : trackingManagers)
  {
    delete tm;
  }
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::SetCuts()
{
  if(!isSetDefaultCutValue)
  {
    SetDefaultCutValue(defaultCutValue);
  }

#ifdef G4VERBOSE
  if(verboseLevel > 1)
  {
    G4cout << "G4VUserPhysicsList::SetCuts:   " << G4endl;
    G4cout << "Cut for gamma: " << GetCutValue("gamma") / mm << "[mm]"
           << G4endl;
    G4cout << "Cut  for e-: " << GetCutValue("e-") / mm << "[mm]" << G4endl;
    G4cout << "Cut  for e+: " << GetCutValue("e+") / mm << "[mm]" << G4endl;
    G4cout << "Cut  for proton: " << GetCutValue("proton") / mm << "[mm]"
           << G4endl;
  }
#endif

  // dump Cut values if verboseLevel==3
  if(verboseLevel > 2)
  {
    DumpCutValuesTable();
  }
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::SetDefaultCutValue(G4double value)
{
  if(value < 0.0)
  {
#ifdef G4VERBOSE
    if(verboseLevel > 0)
    {
      G4cout << "G4VUserPhysicsList::SetDefaultCutValue: negative cut values"
             << "  :" << value / mm << "[mm]" << G4endl;
    }
#endif
    return;
  }

  defaultCutValue      = value;
  isSetDefaultCutValue = true;

  // set cut values for gamma at first and for e- and e+
  SetCutValue(defaultCutValue, "gamma");
  SetCutValue(defaultCutValue, "e-");
  SetCutValue(defaultCutValue, "e+");
  SetCutValue(defaultCutValue, "proton");

#ifdef G4VERBOSE
  if(verboseLevel > 1)
  {
    G4cout << "G4VUserPhysicsList::SetDefaultCutValue:"
           << "default cut value is changed to   :" << defaultCutValue / mm
           << "[mm]" << G4endl;
  }
#endif
}

// --------------------------------------------------------------------
G4double G4VUserPhysicsList::GetCutValue(const G4String& name) const
{
  std::size_t nReg = (G4RegionStore::GetInstance())->size();
  if(nReg == 0)
  {
#ifdef G4VERBOSE
    if(verboseLevel > 0)
    {
      G4cout << "G4VUserPhysicsList::GetCutValue "
             << " : No Default Region " << G4endl;
    }
#endif
    G4Exception("G4VUserPhysicsList::GetCutValue", "Run0253", FatalException,
                "No Default Region");
    return -1. * mm;
  }
  G4Region* region =
    G4RegionStore::GetInstance()->GetRegion("DefaultRegionForTheWorld", false);
  return region->GetProductionCuts()->GetProductionCut(name);
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::SetCutValue(G4double aCut, const G4String& name)
{
  SetParticleCuts(aCut, name);
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::SetCutValue(G4double aCut, const G4String& pname,
                                     const G4String& rname)
{
  G4Region* region = G4RegionStore::GetInstance()->GetRegion(rname);
  if(region != nullptr)
  {
    // set cut value
    SetParticleCuts(aCut, pname, region);
  }
  else
  {
#ifdef G4VERBOSE
    if(verboseLevel > 0)
    {
      G4cout << "G4VUserPhysicsList::SetCutValue "
             << " : No Region of " << rname << G4endl;
    }
#endif
  }
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::SetCutsWithDefault()
{
  SetDefaultCutValue(defaultCutValue);
  G4VUserPhysicsList::SetCuts();
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::SetCutsForRegion(G4double aCut, const G4String& rname)
{
  // set cut values for gamma at first and for e- and e+
  SetCutValue(aCut, "gamma", rname);
  SetCutValue(aCut, "e-", rname);
  SetCutValue(aCut, "e+", rname);
  SetCutValue(aCut, "proton", rname);
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::SetParticleCuts(G4double cut,
                                         G4ParticleDefinition* particle,
                                         G4Region* region)
{
  SetParticleCuts(cut, particle->GetParticleName(), region);
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::SetParticleCuts(G4double cut,
                                         const G4String& particleName,
                                         G4Region* region)
{
  if(cut < 0.0)
  {
#ifdef G4VERBOSE
    if(verboseLevel > 0)
    {
      G4cout << "G4VUserPhysicsList::SetParticleCuts: negative cut values"
             << "  :" << cut / mm << "[mm]"
             << " for " << particleName << G4endl;
    }
#endif
    return;
  }

  G4Region* world_region =
    G4RegionStore::GetInstance()->GetRegion("DefaultRegionForTheWorld", false);
  if(region == nullptr)
  {
    std::size_t nReg = G4RegionStore::GetInstance()->size();
    if(nReg == 0)
    {
#ifdef G4VERBOSE
      if(verboseLevel > 0)
      {
        G4cout << "G4VUserPhysicsList::SetParticleCuts "
               << " : No Default Region " << G4endl;
      }
#endif
      G4Exception("G4VUserPhysicsList::SetParticleCuts ", "Run0254",
                  FatalException, "No Default Region");
      return;
    }
    region = world_region;
  }

  if(!isSetDefaultCutValue)
  {
    SetDefaultCutValue(defaultCutValue);
  }

  G4ProductionCuts* pcuts = region->GetProductionCuts();
  if(region != world_region &&
     pcuts == G4ProductionCutsTable::GetProductionCutsTable()
                ->GetDefaultProductionCuts())
  {  // This region had no unique cuts yet but shares the default cuts.
     // Need to create a new object before setting the value.
    pcuts =
      new G4ProductionCuts(*(G4ProductionCutsTable::GetProductionCutsTable()
                               ->GetDefaultProductionCuts()));
    region->SetProductionCuts(pcuts);
  }
  pcuts->SetProductionCut(cut, particleName);
#ifdef G4VERBOSE
  if(verboseLevel > 2)
  {
    G4cout << "G4VUserPhysicsList::SetParticleCuts: "
           << "  :" << cut / mm << "[mm]"
           << " for " << particleName << G4endl;
  }
#endif
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::BuildPhysicsTable()
{
  // Prepare Physics table for all particles
  theParticleIterator->reset();
  while((*theParticleIterator)())
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    PreparePhysicsTable(particle);
  }

  // ask processes to prepare physics table
  if(fRetrievePhysicsTable)
  {
    fIsRestoredCutValues =
      fCutsTable->RetrieveCutsTable(directoryPhysicsTable, fStoredInAscii);
    // check if retrieve Cut Table successfully
    if(!fIsRestoredCutValues)
    {
#ifdef G4VERBOSE
      if(verboseLevel > 0)
      {
        G4cout << "G4VUserPhysicsList::BuildPhysicsTable"
               << " Retrieve Cut Table failed !!" << G4endl;
      }
#endif
      G4Exception("G4VUserPhysicsList::BuildPhysicsTable", "Run0255",
                  RunMustBeAborted, "Fail to retrieve Production Cut Table");
    }
    else
    {
#ifdef G4VERBOSE
      if(verboseLevel > 2)
      {
        G4cout << "G4VUserPhysicsList::BuildPhysicsTable"
               << "  Retrieve Cut Table successfully " << G4endl;
      }
#endif
    }
  }
  else
  {
#ifdef G4VERBOSE
    if(verboseLevel > 2)
    {
      G4cout << "G4VUserPhysicsList::BuildPhysicsTable"
             << " does not retrieve Cut Table but calculate " << G4endl;
    }
#endif
  }

  // Sets a value to particle
  // set cut values for gamma at first and for e- and e+
  G4String particleName;
  G4ParticleDefinition* GammaP = theParticleTable->FindParticle("gamma");
  if(GammaP)
    BuildPhysicsTable(GammaP);
  G4ParticleDefinition* EMinusP = theParticleTable->FindParticle("e-");
  if(EMinusP)
    BuildPhysicsTable(EMinusP);
  G4ParticleDefinition* EPlusP = theParticleTable->FindParticle("e+");
  if(EPlusP)
    BuildPhysicsTable(EPlusP);
  G4ParticleDefinition* ProtonP = theParticleTable->FindParticle("proton");
  if(ProtonP)
    BuildPhysicsTable(ProtonP);

  theParticleIterator->reset();
  while((*theParticleIterator)())
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    if(particle != GammaP && particle != EMinusP && particle != EPlusP &&
       particle != ProtonP)
    {
      BuildPhysicsTable(particle);
    }
  }

  // Set flag
  fIsPhysicsTableBuilt = true;
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::BuildPhysicsTable(G4ParticleDefinition* particle)
{
  if (auto *trackingManager = particle->GetTrackingManager())
  {
#ifdef G4VERBOSE
    if(verboseLevel > 2)
    {
      G4cout << "G4VUserPhysicsList::BuildPhysicsTable  "
             << "Calculate Physics Table for " << particle->GetParticleName()
             << " via custom TrackingManager" << G4endl;
    }
#endif
    trackingManager->BuildPhysicsTable(*particle);
    return;
  }

  // Change in order to share physics tables for two kind of process.

  if(particle->GetMasterProcessManager() == nullptr)
  {
#ifdef G4VERBOSE
    if(verboseLevel > 0)
    { G4cout
      << "#### G4VUserPhysicsList::BuildPhysicsTable() - BuildPhysicsTable("
      << particle->GetParticleName() << ") skipped..." << G4endl; }
#endif
    return;
  }
  if(fRetrievePhysicsTable)
  {
    if(!fIsRestoredCutValues)
    {
      // fail to retrieve cut tables
#ifdef G4VERBOSE
      if(verboseLevel > 0)
      {
        G4cout << "G4VUserPhysicsList::BuildPhysicsTable  "
               << "Physics table can not be retrieved and will be calculated "
               << G4endl;
      }
#endif
      fRetrievePhysicsTable = false;
    }
    else
    {
#ifdef G4VERBOSE
      if(verboseLevel > 2)
      {
        G4cout << "G4VUserPhysicsList::BuildPhysicsTable  "
               << " Retrieve Physics Table for " << particle->GetParticleName()
               << G4endl;
      }
#endif
      //  Retrieve PhysicsTable from files for proccesses
      RetrievePhysicsTable(particle, directoryPhysicsTable, fStoredInAscii);
    }
  }

#ifdef G4VERBOSE
  if(verboseLevel > 2)
  {
    G4cout << "G4VUserPhysicsList::BuildPhysicsTable  "
           << "Calculate Physics Table for " << particle->GetParticleName()
           << G4endl;
  }
#endif
  // Rebuild the physics tables for every process for this particle type
  // if particle is not ShortLived
  if(!particle->IsShortLived())
  {
    G4ProcessManager* pManager = particle->GetProcessManager();
    if(pManager == nullptr)
    {
#ifdef G4VERBOSE
      if(verboseLevel > 0)
      {
        G4cout << "G4VUserPhysicsList::BuildPhysicsTable "
               << " : No Process Manager for " << particle->GetParticleName()
               << G4endl;
        G4cout << particle->GetParticleName()
               << " should be created in your PhysicsList" << G4endl;
      }
#endif
      G4Exception("G4VUserPhysicsList::BuildPhysicsTable", "Run0271",
                  FatalException, "No process manager");
      return;
    }

    // Get processes from master thread;
    G4ProcessManager* pManagerShadow = particle->GetMasterProcessManager();

    G4ProcessVector* pVector = pManager->GetProcessList();
    if(pVector == nullptr)
    {
#ifdef G4VERBOSE
      if(verboseLevel > 0)
      {
        G4cout << "G4VUserPhysicsList::BuildPhysicsTable  "
               << " : No Process Vector for " << particle->GetParticleName()
               << G4endl;
      }
#endif
      G4Exception("G4VUserPhysicsList::BuildPhysicsTable", "Run0272",
                  FatalException, "No process Vector");
      return;
    }
#ifdef G4VERBOSE
    if(verboseLevel > 2)
    {
      G4cout << "G4VUserPhysicsList::BuildPhysicsTable %%%%%% "
             << particle->GetParticleName() << G4endl;
      G4cout << " ProcessManager : " << pManager
             << " ProcessManagerShadow : " << pManagerShadow << G4endl;
      for(G4int iv1 = 0; iv1 < (G4int)pVector->size(); ++iv1)
      {
        G4cout << "  " << iv1 << " - " << (*pVector)[iv1]->GetProcessName()
               << G4endl;
      }
      G4cout << "--------------------------------------------------------------"
             << G4endl;
      G4ProcessVector* pVectorShadow = pManagerShadow->GetProcessList();

      for(G4int iv2 = 0; iv2 < (G4int)pVectorShadow->size(); ++iv2)
      {
        G4cout << "  " << iv2 << " - "
               << (*pVectorShadow)[iv2]->GetProcessName() << G4endl;
      }
    }
#endif
    for(G4int j = 0; j < (G4int)pVector->size(); ++j)
    {
      // Andrea July 16th 2013 : migration to new interface...
      // Infer if we are in a worker thread or master thread
      // Master thread is the one in which the process manager
      // and process manager shadow pointers are the same
      if(pManagerShadow == pManager)
      {
        (*pVector)[j]->BuildPhysicsTable(*particle);
      }
      else
      {
        (*pVector)[j]->BuildWorkerPhysicsTable(*particle);
      }

    }  // End loop on processes vector
  }    // End if short-lived
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::PreparePhysicsTable(G4ParticleDefinition* particle)
{
  if (auto *trackingManager = particle->GetTrackingManager())
  {
    trackingManager->PreparePhysicsTable(*particle);
    return;
  }

  if(!(particle->GetMasterProcessManager()))
  {
    return;
  }
  // Prepare the physics tables for every process for this particle type
  // if particle is not ShortLived
  if(!particle->IsShortLived())
  {
    G4ProcessManager* pManager = particle->GetProcessManager();
    if(pManager == nullptr)
    {
#ifdef G4VERBOSE
      if(verboseLevel > 0)
      {
        G4cout << "G4VUserPhysicsList::PreparePhysicsTable  "
               << ": No Process Manager for " << particle->GetParticleName()
               << G4endl;
        G4cout << particle->GetParticleName()
               << " should be created in your PhysicsList" << G4endl;
      }
#endif
      G4Exception("G4VUserPhysicsList::PreparePhysicsTable", "Run0273",
                  FatalException, "No process manager");
      return;
    }

    // Get processes from master thread
    G4ProcessManager* pManagerShadow = particle->GetMasterProcessManager();

    G4ProcessVector* pVector = pManager->GetProcessList();
    if(pVector == nullptr)
    {
#ifdef G4VERBOSE
      if(verboseLevel > 0)
      {
        G4cout << "G4VUserPhysicsList::PreparePhysicsTable  "
               << ": No Process Vector for " << particle->GetParticleName()
               << G4endl;
      }
#endif
      G4Exception("G4VUserPhysicsList::PreparePhysicsTable", "Run0274",
                  FatalException, "No process Vector");
      return;
    }
    for(G4int j = 0; j < (G4int)pVector->size(); ++j)
    {
      // Andrea July 16th 2013 : migration to new interface...
      // Infer if we are in a worker thread or master thread
      // Master thread is the one in which the process manager
      // and process manager shadow pointers are the same
      if(pManagerShadow == pManager)
      {
        (*pVector)[j]->PreparePhysicsTable(*particle);
      }
      else
      {
        (*pVector)[j]->PrepareWorkerPhysicsTable(*particle);
      }
    }  // End loop on processes vector
  }    // End if pn ShortLived
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::BuildIntegralPhysicsTable(
  G4VProcess* process, G4ParticleDefinition* particle)
{
  // TODO Should we change this function?
  //*******************************************************************
  // Temporary addition to make the integral schema of electromagnetic
  // processes work.

  if((process->GetProcessName() == "Imsc") ||
     (process->GetProcessName() == "IeIoni") ||
     (process->GetProcessName() == "IeBrems") ||
     (process->GetProcessName() == "Iannihil") ||
     (process->GetProcessName() == "IhIoni") ||
     (process->GetProcessName() == "IMuIoni") ||
     (process->GetProcessName() == "IMuBrems") ||
     (process->GetProcessName() == "IMuPairProd"))
  {
#ifdef G4VERBOSE
    if(verboseLevel > 2)
    {
      G4cout << "G4VUserPhysicsList::BuildIntegralPhysicsTable  "
             << " BuildPhysicsTable is invoked for "
             << process->GetProcessName() << "(" << particle->GetParticleName()
             << ")" << G4endl;
    }
#endif
    process->BuildPhysicsTable(*particle);
  }
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::DumpList() const
{
  theParticleIterator->reset();
  G4int idx = 0;
  while((*theParticleIterator)())
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4cout << particle->GetParticleName();
    if((idx++ % 4) == 3)
    {
      G4cout << G4endl;
    }
    else
    {
      G4cout << ", ";
    }
  }
  G4cout << G4endl;
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::DumpCutValuesTable(G4int flag)
{
  fDisplayThreshold = flag;
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::DumpCutValuesTableIfRequested()
{
  if(fDisplayThreshold == 0)
    return;
  G4ProductionCutsTable::GetProductionCutsTable()->DumpCouples();
  fDisplayThreshold = 0;
}

// --------------------------------------------------------------------
G4bool G4VUserPhysicsList::StorePhysicsTable(const G4String& directory)
{
  G4bool ascii = fStoredInAscii;
  G4String dir = directory;
  if(dir.empty())
    dir = directoryPhysicsTable;
  else
    directoryPhysicsTable = dir;

  // store CutsTable info
  if(!fCutsTable->StoreCutsTable(dir, ascii))
  {
    G4Exception("G4VUserPhysicsList::StorePhysicsTable", "Run0281", JustWarning,
                "Fail to store Cut Table");
    return false;
  }
#ifdef G4VERBOSE
  if(verboseLevel > 2)
  {
    G4cout << "G4VUserPhysicsList::StorePhysicsTable   "
           << " Store material and cut values successfully" << G4endl;
  }
#endif

  G4bool success = true;

  // loop over all particles in G4ParticleTable
  theParticleIterator->reset();
  while((*theParticleIterator)())
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    // Store physics tables for every process for this particle type
    G4ProcessVector* pVector =
      (particle->GetProcessManager())->GetProcessList();
    for(G4int j = 0; j < (G4int)pVector->size(); ++j)
    {
      if(!(*pVector)[j]->StorePhysicsTable(particle, dir, ascii))
      {
        G4String comment = "Fail to store physics table for ";
        comment += (*pVector)[j]->GetProcessName();
        comment += "(" + particle->GetParticleName() + ")";
        G4Exception("G4VUserPhysicsList::StorePhysicsTable", "Run0282",
                    JustWarning, comment);
        success = false;
      }
    }
    // end loop over processes
  }
  // end loop over particles
  return success;
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::SetPhysicsTableRetrieved(const G4String& directory)
{
  fRetrievePhysicsTable = true;
  if(!directory.empty())
  {
    directoryPhysicsTable = directory;
  }
  fIsCheckedForRetrievePhysicsTable = false;
  fIsRestoredCutValues              = false;
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::RetrievePhysicsTable(G4ParticleDefinition* particle,
                                              const G4String& directory,
                                              G4bool ascii)
{
  G4bool success[100];
  // Retrieve physics tables for every process for this particle type
  G4ProcessVector* pVector = (particle->GetProcessManager())->GetProcessList();
  for(G4int j = 0; j < (G4int)pVector->size(); ++j)
  {
    success[j] =
      (*pVector)[j]->RetrievePhysicsTable(particle, directory, ascii);

    if(!success[j])
    {
#ifdef G4VERBOSE
      if(verboseLevel > 2)
      {
        G4cout << "G4VUserPhysicsList::RetrievePhysicsTable   "
               << " Fail to retrieve Physics Table for "
               << (*pVector)[j]->GetProcessName() << G4endl;
        G4cout << "Calculate Physics Table for " << particle->GetParticleName()
               << G4endl;
      }
#endif
      (*pVector)[j]->BuildPhysicsTable(*particle);
    }
  }
  for(G4int j = 0; j < (G4int)pVector->size(); ++j)
  {
    // temporary addition to make the integral schema
    if(!success[j])
      BuildIntegralPhysicsTable((*pVector)[j], particle);
  }
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::SetApplyCuts(G4bool value, const G4String& name)
{
#ifdef G4VERBOSE
  if(verboseLevel > 2)
  {
    G4cout << "G4VUserPhysicsList::SetApplyCuts for " << name << G4endl;
  }
#endif
  if(name == "all")
  {
    theParticleTable->FindParticle("gamma")->SetApplyCutsFlag(value);
    theParticleTable->FindParticle("e-")->SetApplyCutsFlag(value);
    theParticleTable->FindParticle("e+")->SetApplyCutsFlag(value);
    theParticleTable->FindParticle("proton")->SetApplyCutsFlag(value);
  }
  else
  {
    theParticleTable->FindParticle(name)->SetApplyCutsFlag(value);
  }
}

// --------------------------------------------------------------------
G4bool G4VUserPhysicsList::GetApplyCuts(const G4String& name) const
{
  return theParticleTable->FindParticle(name)->GetApplyCutsFlag();
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::CheckParticleList()
{
  if(!fDisableCheckParticleList)
  {
    G4MT_thePLHelper->CheckParticleList();
  }
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::AddTransportation()
{
  G4MT_thePLHelper->AddTransportation();
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::UseCoupledTransportation(G4bool vl)
{
  G4MT_thePLHelper->UseCoupledTransportation(vl);
}

// --------------------------------------------------------------------
G4bool G4VUserPhysicsList::RegisterProcess(G4VProcess* process,
                                           G4ParticleDefinition* particle)
{
  return G4MT_thePLHelper->RegisterProcess(process, particle);
}

// --------------------------------------------------------------------
G4ParticleTable::G4PTblDicIterator*
G4VUserPhysicsList::GetParticleIterator() const
{
  return (subInstanceManager.offset[g4vuplInstanceID])._theParticleIterator;
}

// --------------------------------------------------------------------
void G4VUserPhysicsList::SetVerboseLevel(G4int value)
{
  verboseLevel = value;
  // set verboseLevel for G4ProductionCutsTable same as one for
  // G4VUserPhysicsList:
  fCutsTable->SetVerboseLevel(verboseLevel);

  G4MT_thePLHelper->SetVerboseLevel(verboseLevel);

#ifdef G4VERBOSE
  if(verboseLevel > 1)
  {
    G4cout << "G4VUserPhysicsList::SetVerboseLevel  :"
           << " Verbose level is set to " << verboseLevel << G4endl;
  }
#endif
}
