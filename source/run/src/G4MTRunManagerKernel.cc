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
// G4MTRunManagerKernel implementation
//
// Author: M.Asai - July 2013
// --------------------------------------------------------------------
#include "G4MTRunManagerKernel.hh"

#include "G4AutoLock.hh"
#include "G4DecayTable.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTableIterator.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PhysicsVector.hh"
#include "G4PolyconeSide.hh"
#include "G4PolyhedraSide.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4StateManager.hh"
#include "G4UImanager.hh"
#include "G4UserWorkerInitialization.hh"
#include "G4UserWorkerThreadInitialization.hh"
#include "G4VDecayChannel.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserPhysicsList.hh"
#include "G4WorkerRunManager.hh"
#include "G4WorkerThread.hh"

std::vector<G4WorkerRunManager*>* G4MTRunManagerKernel::workerRMvector = nullptr;

G4ThreadLocal G4WorkerThread* G4MTRunManagerKernel::wThreadContext = nullptr;

// --------------------------------------------------------------------
namespace
{
G4Mutex workerRMMutex = G4MUTEX_INITIALIZER;
}

// --------------------------------------------------------------------
G4MTRunManagerKernel::G4MTRunManagerKernel() : G4RunManagerKernel(masterRMK)
{
  // This version of the constructor should never be called in sequential mode!
#ifndef G4MULTITHREADED
  G4ExceptionDescription msg;
  msg << "Geant4 code is compiled without multi-threading support "
         "(-DG4MULTITHREADED "
         "is set to off).";
  msg << " This type of RunManager can only be used in mult-threaded "
         "applications.";
  G4Exception("G4RunManagerKernel::G4RunManagerKernel()", "Run0109", FatalException, msg);
#endif
  G4AutoLock l(&workerRMMutex);
  if (workerRMvector == nullptr) workerRMvector = new std::vector<G4WorkerRunManager*>;
  l.unlock();
  // Set flag that a MT-type kernel has been instantiated
  G4Threading::SetMultithreadedApplication(true);
}

// --------------------------------------------------------------------
G4MTRunManagerKernel::~G4MTRunManagerKernel()
{
  G4AutoLock l(&workerRMMutex);
  if (workerRMvector != nullptr) {
    if (!workerRMvector->empty()) {
      G4ExceptionDescription msg;
      msg << "G4MTRunManagerKernel is to be deleted while " << workerRMvector->size()
          << " G4WorkerRunManager are still alive.";
      G4Exception("G4RunManagerKernel::~G4RunManagerKernel()", "Run10035", FatalException, msg);
    }
    delete workerRMvector;
    workerRMvector = nullptr;
  }
}

// --------------------------------------------------------------------
void G4MTRunManagerKernel::SetupShadowProcess() const
{
  // Behavior is the same as base class (sequential mode)
  // ShadowProcess pointer == process poitner
  G4RunManagerKernel::SetupShadowProcess();
}

// --------------------------------------------------------------------
G4WorkerThread* G4MTRunManagerKernel::GetWorkerThread()
{
  return wThreadContext;
}

// --------------------------------------------------------------------
void G4MTRunManagerKernel::StartThread(G4WorkerThread* context)
{
  //!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!!!!!! IMPORTANT !!!!!!!!!
  //!!!!!!!!!!!!!!!!!!!!!!!!!!
  // Here is not sequential anymore and G4UserWorkerThreadInitialization is
  // a shared user initialization class.
  // This means this method cannot use data members of G4RunManagerKernel
  // unless they are invariant ("read-only") and can be safely shared.
  // All the rest that is not invariant should be incapsulated into
  // the context (or, as for wThreadContext be G4ThreadLocal)
  //!!!!!!!!!!!!!!!!!!!!!!!!!!
  // #ifdef G4MULTITHREADED
  //     turnontpmalloc();
  // #endif
  G4Threading::WorkerThreadJoinsPool();
  wThreadContext = context;
  G4MTRunManager* masterRM = G4MTRunManager::GetMasterRunManager();

  //============================
  // Step-0: Thread ID
  //============================
  // Initialise per-thread stream-output
  // The following line is needed before we actually do IO initialisation
  // because the constructor of UI manager resets the IO destination.
  G4int thisID = wThreadContext->GetThreadId();
  G4Threading::G4SetThreadId(thisID);
  G4UImanager::GetUIpointer()->SetUpForAThread(thisID);

  //============================
  // Optimization: optional
  //============================
  // Enforce thread affinity if requested
  wThreadContext->SetPinAffinity(masterRM->GetPinAffinity());

  //============================
  // Step-1: Random number engine
  //============================
  // RNG Engine needs to be initialised by "cloning" the master one.
  const CLHEP::HepRandomEngine* masterEngine = masterRM->getMasterRandomEngine();
  masterRM->GetUserWorkerThreadInitialization()->SetupRNGEngine(masterEngine);

  //============================
  // Step-2: Initialise worker thread
  //============================
  if (masterRM->GetUserWorkerInitialization() != nullptr) {
    masterRM->GetUserWorkerInitialization()->WorkerInitialize();
  }
  if (masterRM->GetUserActionInitialization() != nullptr) {
    G4VSteppingVerbose* sv = masterRM->GetUserActionInitialization()->InitializeSteppingVerbose();
    if (sv != nullptr) {
      G4VSteppingVerbose::SetInstance(sv);
    }
  }
  // Now initialise worker part of shared objects (geometry/physics)
  wThreadContext->BuildGeometryAndPhysicsVector();
  G4WorkerRunManager* wrm = masterRM->GetUserWorkerThreadInitialization()->CreateWorkerRunManager();
  wrm->SetWorkerThread(wThreadContext);
  G4AutoLock wrmm(&workerRMMutex);
  workerRMvector->push_back(wrm);
  wrmm.unlock();

  //================================
  // Step-3: Setup worker run manager
  //================================
  // Set the detector and physics list to the worker thread. Share with master
  const G4VUserDetectorConstruction* detector = masterRM->GetUserDetectorConstruction();
  wrm->G4RunManager::SetUserInitialization(const_cast<G4VUserDetectorConstruction*>(detector));
  const G4VUserPhysicsList* physicslist = masterRM->GetUserPhysicsList();
  wrm->SetUserInitialization(const_cast<G4VUserPhysicsList*>(physicslist));

  //================================
  // Step-4: Initialise worker run manager
  //================================
  if (masterRM->GetUserActionInitialization() != nullptr) {
    masterRM->GetNonConstUserActionInitialization()->Build();
  }
  if (masterRM->GetUserWorkerInitialization() != nullptr) {
    masterRM->GetUserWorkerInitialization()->WorkerStart();
  }
  wrm->Initialize();

  //================================
  // Step5: Loop over requests from the master thread
  //================================
  // This function should enter a loop processing new runs and actions
  // requests from master. It should block until thread is ready
  // to terminate
  wrm->DoWork();

  //===============================
  // Step-6: Terminate worker thread
  //===============================
  if (masterRM->GetUserWorkerInitialization() != nullptr) {
    masterRM->GetUserWorkerInitialization()->WorkerStop();
  }

  wrmm.lock();
  for (auto itrWrm = workerRMvector->cbegin(); itrWrm != workerRMvector->cend(); ++itrWrm) {
    if ((*itrWrm) == wrm) {
      workerRMvector->erase(itrWrm);
      break;
    }
  }
  wrmm.unlock();
  delete wrm;

  //===============================
  // Step-7: Cleanup split classes
  //===============================
  wThreadContext->DestroyGeometryAndPhysicsVector();
  wThreadContext = nullptr;

  G4Threading::WorkerThreadLeavesPool();
}

// --------------------------------------------------------------------
void G4MTRunManagerKernel::SetUpDecayChannels()
{
  auto pItr = G4ParticleTable::GetParticleTable()->GetIterator();
  pItr->reset();
  while ((*pItr)()) {
    G4DecayTable* dt = pItr->value()->GetDecayTable();
    if (dt != nullptr) {
      G4int nCh = dt->entries();
      for (G4int i = 0; i < nCh; ++i) {
        dt->GetDecayChannel(i)->GetDaughter(0);
      }
    }
  }
}

// --------------------------------------------------------------------
void G4MTRunManagerKernel::BroadcastAbortRun(G4bool softAbort)
{
  G4AutoLock wrmm(&workerRMMutex);

  for (const auto& itr : *workerRMvector) {
    itr->AbortRun(softAbort);
  }
}
