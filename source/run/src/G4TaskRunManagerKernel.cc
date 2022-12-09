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

#include "G4TaskRunManagerKernel.hh"
#include "G4AutoLock.hh"
#include "G4RegionStore.hh"
#include "G4StateManager.hh"

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
#include "G4Run.hh"
#include "G4TaskManager.hh"
#include "G4TiMemory.hh"
#include "G4UImanager.hh"
#include "G4UserWorkerInitialization.hh"
#include "G4UserWorkerThreadInitialization.hh"
#include "G4VDecayChannel.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserPhysicsList.hh"
#include "G4WorkerTaskRunManager.hh"
#include "G4WorkerThread.hh"

#include <atomic>
#include <memory>

//============================================================================//

std::vector<G4String> G4TaskRunManagerKernel::initCmdStack = {};

//============================================================================//

G4TaskRunManagerKernel::G4TaskRunManagerKernel()
  : G4RunManagerKernel(masterRMK)
{
  // This version of the constructor should never be called in sequential mode!
#ifndef G4MULTITHREADED
  G4ExceptionDescription msg;
  msg << "Geant4 code is compiled without multi-threading support "
         "(-DG4MULTITHREADED "
         "is set to off).";
  msg << " This type of RunManager can only be used in mult-threaded "
         "applications.";
  G4Exception("G4RunManagerKernel::G4RunManagerKernel()", "Run0109",
              FatalException, msg);
#endif
  // Set flag that a MT-type kernel has been instantiated
  G4Threading::SetMultithreadedApplication(true);
}

//============================================================================//

G4TaskRunManagerKernel::~G4TaskRunManagerKernel() {}

//============================================================================//

void G4TaskRunManagerKernel::SetupShadowProcess() const
{
  // Behavior is the same as base class (sequential mode)
  // ShadowProcess pointer == process poitner
  G4RunManagerKernel::SetupShadowProcess();
}

//============================================================================//

namespace
{
  using WorkerRunManPtr_t = std::unique_ptr<G4WorkerTaskRunManager>;
  using WorkerThreadPtr_t = std::unique_ptr<G4WorkerThread>;

  WorkerRunManPtr_t& workerRM()
  {
    G4ThreadLocalStatic WorkerRunManPtr_t _instance{ nullptr };
    return _instance;
  }

  WorkerThreadPtr_t& context()
  {
    G4ThreadLocalStatic WorkerThreadPtr_t _instance{ nullptr };
    return _instance;
  }

}  // namespace

//============================================================================//

G4WorkerThread* G4TaskRunManagerKernel::GetWorkerThread()
{
  return context().get();
}

//============================================================================//

void G4TaskRunManagerKernel::InitializeWorker()
{
  if(context() && workerRM())
    return;

  G4TaskRunManager* mrm = G4TaskRunManager::GetMasterRunManager();
  if(G4MTRunManager::GetMasterThreadId() == G4ThisThread::get_id())
  {
    G4TaskManager* taskManager = mrm->GetTaskManager();
    auto _fut                  = taskManager->async(InitializeWorker);
    _fut->wait();
    return;
  }

  //!!!!!!!!!!!!!!!!!!!!!!!!!!
  //!!!!!! IMPORTANT !!!!!!!!!
  //!!!!!!!!!!!!!!!!!!!!!!!!!!
  // Here is not sequential anymore and G4UserWorkerThreadInitialization is
  // a shared user initialization class
  // This means this method cannot use data memebers of G4RunManagerKernel
  // unless they are invariant ("read-only") and can be safely shared.
  //  All the rest that is not invariant should be incapsualted into
  //  the context (or, as for wThreadContext be G4ThreadLocal)
  //!!!!!!!!!!!!!!!!!!!!!!!!!!

  G4Threading::WorkerThreadJoinsPool();
  context().reset(new G4WorkerThread);

  //============================
  // Step-0: Thread ID
  //============================
  // Initliazie per-thread stream-output
  // The following line is needed before we actually do IO initialization
  // becasue the constructor of UI manager resets the IO destination.
  context()->SetNumberThreads((G4int)mrm->GetThreadPool()->size());
  context()->SetThreadId(G4int(G4ThreadPool::get_this_thread_id() - 1));
  G4int thisID = context()->GetThreadId();
  G4Threading::G4SetThreadId(thisID);
  G4UImanager::GetUIpointer()->SetUpForAThread(thisID);

  //============================
  // Optimization: optional
  //============================
  // Enforce thread affinity if requested
  context()->SetPinAffinity(mrm->GetPinAffinity());

  //============================
  // Step-1: Random number engine
  //============================
  // RNG Engine needs to be initialized by "cloning" the master one.
  const CLHEP::HepRandomEngine* masterEngine = mrm->getMasterRandomEngine();
  mrm->GetUserWorkerThreadInitialization()->SetupRNGEngine(masterEngine);

  //============================
  // Step-2: Initialize worker thread
  //============================
  if(mrm->GetUserWorkerInitialization())
    mrm->GetUserWorkerInitialization()->WorkerInitialize();

  if(mrm->GetUserActionInitialization())
  {
    G4VSteppingVerbose* sv =
      mrm->GetUserActionInitialization()->InitializeSteppingVerbose();
    if(sv)
      G4VSteppingVerbose::SetInstance(sv);
  }
  // Now initialize worker part of shared objects (geometry/physics)
  context()->BuildGeometryAndPhysicsVector();
  workerRM().reset(static_cast<G4WorkerTaskRunManager*>(
    mrm->GetUserWorkerThreadInitialization()->CreateWorkerRunManager()));
  auto& wrm = workerRM();
  wrm->SetWorkerThread(context().get());

  //================================
  // Step-3: Setup worker run manager
  //================================
  // Set the detector and physics list to the worker thread. Share with master
  const G4VUserDetectorConstruction* detector =
    mrm->GetUserDetectorConstruction();
  wrm->G4RunManager::SetUserInitialization(
    const_cast<G4VUserDetectorConstruction*>(detector));
  const G4VUserPhysicsList* physicslist = mrm->GetUserPhysicsList();
  wrm->SetUserInitialization(const_cast<G4VUserPhysicsList*>(physicslist));

  //================================
  // Step-4: Initialize worker run manager
  //================================
  if(mrm->GetUserActionInitialization())
    mrm->GetNonConstUserActionInitialization()->Build();
  if(mrm->GetUserWorkerInitialization())
    mrm->GetUserWorkerInitialization()->WorkerStart();

  workerRM()->Initialize();

  for(auto& itr : initCmdStack)
    G4UImanager::GetUIpointer()->ApplyCommand(itr);

  wrm->ProcessUI();
}

//============================================================================//

void G4TaskRunManagerKernel::ExecuteWorkerInit()
{
  // because of TBB
  if(G4MTRunManager::GetMasterThreadId() == G4ThisThread::get_id())
  {
    G4TaskManager* taskManager =
      G4TaskRunManager::GetMasterRunManager()->GetTaskManager();
    auto _fut = taskManager->async(ExecuteWorkerInit);
    return _fut->get();
  }

  // this check is for TBB as there is not a way to run an initialization
  // routine on each thread
  if(!workerRM())
    InitializeWorker();

  auto& wrm = workerRM();
  assert(wrm.get() != nullptr);
  wrm->DoCleanup();
}

//============================================================================//

void G4TaskRunManagerKernel::ExecuteWorkerTask()
{
  // because of TBB
  if(G4MTRunManager::GetMasterThreadId() == G4ThisThread::get_id())
  {
    G4TaskManager* taskManager =
      G4TaskRunManager::GetMasterRunManager()->GetTaskManager();
    auto _fut = taskManager->async(ExecuteWorkerTask);
    return _fut->get();
  }

  // this check is for TBB as there is not a way to run an initialization
  // routine on each thread
  if(!workerRM())
    InitializeWorker();

  auto& wrm = workerRM();
  assert(wrm.get() != nullptr);
  wrm->DoWork();
}

//============================================================================//

void G4TaskRunManagerKernel::TerminateWorkerRunEventLoop()
{
  if(workerRM())
    TerminateWorkerRunEventLoop(workerRM().get());
}

//============================================================================//

void G4TaskRunManagerKernel::TerminateWorker()
{
  if(workerRM())
    TerminateWorker(workerRM().get());
  workerRM().reset();
  context().reset();
}

//============================================================================//

void G4TaskRunManagerKernel::TerminateWorkerRunEventLoop(
  G4WorkerTaskRunManager* wrm)
{
  if(!wrm)
    return;

  wrm->TerminateEventLoop();
  wrm->RunTermination();
}

//============================================================================//

void G4TaskRunManagerKernel::TerminateWorker(G4WorkerTaskRunManager* wrm)
{
  if(!wrm)
    return;

  //===============================
  // Step-6: Terminate worker thread
  //===============================
  G4TaskRunManager* mrm = G4TaskRunManager::GetMasterRunManager();
  if(mrm && mrm->GetUserWorkerInitialization())
    mrm->GetUserWorkerInitialization()->WorkerStop();

  G4WorkerThread* _context = wrm->GetWorkerThread();
  _context->DestroyGeometryAndPhysicsVector();

  G4Threading::WorkerThreadLeavesPool();
}

//============================================================================//

std::vector<G4String>& G4TaskRunManagerKernel::InitCommandStack()
{
  return initCmdStack;
}

//============================================================================//

void G4TaskRunManagerKernel::SetUpDecayChannels()
{
  G4ParticleTable::G4PTblDicIterator* pItr =
    G4ParticleTable::GetParticleTable()->GetIterator();
  pItr->reset();
  while((*pItr)())
  {
    G4DecayTable* dt = pItr->value()->GetDecayTable();
    if(dt)
    {
      G4int nCh = dt->entries();
      for(G4int i = 0; i < nCh; i++)
      {
        dt->GetDecayChannel(i)->GetDaughter(0);
      }
    }
  }
}

//============================================================================//

void G4TaskRunManagerKernel::BroadcastAbortRun(G4bool softAbort)
{
  G4ConsumeParameters(softAbort);
}

//============================================================================//
