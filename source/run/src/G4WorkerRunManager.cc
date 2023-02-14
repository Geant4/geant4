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
// G4WorkerRunManager implementation
//
// Original authors: X.Dong, A.Dotti - 2013
// --------------------------------------------------------------------

#include <fstream>
#include <sstream>

#include "G4WorkerRunManager.hh"
#include "G4MTRunManager.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4RNGHelper.hh"
#include "G4Run.hh"
#include "G4SDManager.hh"
#include "G4ScoringManager.hh"
#include "G4TiMemory.hh"
#include "G4Timer.hh"
#include "G4TransportationManager.hh"
#include "G4UImanager.hh"
#include "G4UserRunAction.hh"
#include "G4UserWorkerInitialization.hh"
#include "G4UserWorkerThreadInitialization.hh"
#include "G4VScoreNtupleWriter.hh"
#include "G4VScoringMesh.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VVisManager.hh"
#include "G4WorkerRunManagerKernel.hh"
#include "G4WorkerThread.hh"
#include "G4MTRunManager.hh"
#include "G4ParallelWorldProcessStore.hh"
#include "G4AutoLock.hh"

// --------------------------------------------------------------------
G4WorkerRunManager* G4WorkerRunManager::GetWorkerRunManager()
{
  return static_cast<G4WorkerRunManager*>(G4RunManager::GetRunManager());
}

// --------------------------------------------------------------------
G4WorkerRunManagerKernel* G4WorkerRunManager::GetWorkerRunManagerKernel()
{
  return static_cast<G4WorkerRunManagerKernel*>(GetWorkerRunManager()->kernel);
}

// --------------------------------------------------------------------
G4WorkerRunManager::G4WorkerRunManager()
  : G4RunManager(workerRM)
{
  // This constructor should never be called in non-multithreaded mode

#ifndef G4MULTITHREADED
  G4ExceptionDescription msg;
  msg << "Geant4 code is compiled without multi-threading support "
         "(-DG4MULTITHREADED "
         "is set to off).";
  msg << " This type of RunManager can only be used in mult-threaded "
         "applications.";
  G4Exception("G4WorkerRunManager::G4WorkerRunManager()", "Run0103",
              FatalException, msg);
#endif
  G4ParticleTable::GetParticleTable()->WorkerG4ParticleTable();
  G4ScoringManager* masterScM = G4MTRunManager::GetMasterScoringManager();
  if(masterScM != nullptr)
    G4ScoringManager::GetScoringManager();  // TLS instance for a worker

  // Properly initialise luxury level for Ranlux* engines...
  //
  if(dynamic_cast<const CLHEP::Ranlux64Engine*>(G4Random::getTheEngine()))
  {
    const CLHEP::Ranlux64Engine* theEngine =
      dynamic_cast<const CLHEP::Ranlux64Engine*>(G4Random::getTheEngine());
    luxury = theEngine->getLuxury();
  }
  else if(dynamic_cast<const CLHEP::RanluxEngine*>(G4Random::getTheEngine()))
  {
    const CLHEP::RanluxEngine* theEngine =
      dynamic_cast<const CLHEP::RanluxEngine*>(G4Random::getTheEngine());
    luxury = theEngine->getLuxury();
  }

  G4UImanager::GetUIpointer()->SetIgnoreCmdNotFound(true);

#ifdef G4MULTITHREADED
  G4VVisManager* pVVis = G4VVisManager::GetConcreteInstance();
  if(pVVis != nullptr)
  {
    pVVis->SetUpForAThread();
    visIsSetUp = true;
  }
#endif
}

// --------------------------------------------------------------------
G4WorkerRunManager::~G4WorkerRunManager()
{
  // Delete thread-local data process manager objects
  if(physicsList)
  {
    // physicsList->TerminateWorker();
    // physicsList->RemoveProcessManager();
  }

  // Put these pointers to zero: owned by master thread
  // If not to zero, the base class destructor will attempt to
  // delete them
  userDetector                   = nullptr;
  userWorkerInitialization       = nullptr;
  userWorkerThreadInitialization = nullptr;
  userActionInitialization       = nullptr;
  physicsList                    = nullptr;
  if(verboseLevel > 1)
    G4cout << "Destroying WorkerRunManager (" << this << ")" << G4endl;
}

// --------------------------------------------------------------------
void G4WorkerRunManager::InitializeGeometry()
{
  if(userDetector == nullptr)
  {
    G4Exception("G4RunManager::InitializeGeometry", "Run0033", FatalException,
                "G4VUserDetectorConstruction is not defined!");
    return;
  }
  if(fGeometryHasBeenDestroyed)
  {
    G4TransportationManager::GetTransportationManager()->ClearParallelWorlds();
  }

  // Step1: Get pointer to the physiWorld (note: needs to get the "super
  // pointer, i.e. the one shared by all threads"
  G4RunManagerKernel* masterKernel =
    G4MTRunManager::GetMasterRunManagerKernel();
  G4VPhysicalVolume* worldVol = masterKernel->GetCurrentWorld();
  // Step2:, Call a new "WorkerDefineWorldVolume( pointer from 2-, false);
  kernel->WorkerDefineWorldVolume(worldVol, false);
  kernel->SetNumberOfParallelWorld(masterKernel->GetNumberOfParallelWorld());
  // Step3: Call user's ConstructSDandField()
  userDetector->ConstructSDandField();
  userDetector->ConstructParallelSD();
  geometryInitialized = true;
}

// --------------------------------------------------------------------
void G4WorkerRunManager::RunInitialization()
{
#ifdef G4MULTITHREADED
  if(!visIsSetUp)
  {
    G4VVisManager* pVVis = G4VVisManager::GetConcreteInstance();
    if(pVVis != nullptr)
    {
      pVVis->SetUpForAThread();
      visIsSetUp = true;
    }
  }
#endif

  if(!(kernel->RunInitialization(fakeRun)))
    return;

  // Signal this thread can start event loop.
  // Note this will return only when all threads reach this point
  G4MTRunManager::GetMasterRunManager()->ThisWorkerReady();
  if(fakeRun)
    return;

  const G4UserWorkerInitialization* uwi =
    G4MTRunManager::GetMasterRunManager()->GetUserWorkerInitialization();
  CleanUpPreviousEvents();
  delete currentRun;
  currentRun = nullptr;

  if(fGeometryHasBeenDestroyed)
    G4ParallelWorldProcessStore::GetInstance()->UpdateWorlds();
  // Call a user hook: this is guaranteed all threads are "synchronized"
  if(uwi)
    uwi->WorkerRunStart();

  if(userRunAction)
    currentRun = userRunAction->GenerateRun();
  if(currentRun == nullptr)
    currentRun = new G4Run();

  currentRun->SetRunID(runIDCounter);
  currentRun->SetNumberOfEventToBeProcessed(numberOfEventToBeProcessed);

  currentRun->SetDCtable(DCtable);
  G4SDManager* fSDM = G4SDManager::GetSDMpointerIfExist();
  if(fSDM != nullptr)
  {
    currentRun->SetHCtable(fSDM->GetHCtable());
  }

  if(G4VScoreNtupleWriter::Instance())
  {
    auto hce            = (fSDM != nullptr) ? fSDM->PrepareNewEvent() : nullptr;
    isScoreNtupleWriter = G4VScoreNtupleWriter::Instance()->Book(hce);
    delete hce;
  }

  std::ostringstream oss;
  G4Random::saveFullState(oss);
  randomNumberStatusForThisRun = oss.str();
  currentRun->SetRandomNumberStatus(randomNumberStatusForThisRun);

  for(G4int i_prev = 0; i_prev < n_perviousEventsToBeStored; ++i_prev)
  {
    previousEvents->push_back(nullptr);
  }

  if(printModulo > 0 || verboseLevel > 0)
  {
    G4cout << "### Run " << currentRun->GetRunID()
           << " starts on worker thread "
           << G4Threading::G4GetThreadId() << "."
           << G4endl;
  }
  if(userRunAction != nullptr)
    userRunAction->BeginOfRunAction(currentRun);

#if defined(GEANT4_USE_TIMEMORY)
  workerRunProfiler.reset(new ProfilerConfig(currentRun));
#endif

  if(isScoreNtupleWriter)
  {
    G4VScoreNtupleWriter::Instance()->OpenFile();
  }

  if(storeRandomNumberStatus)
  {
    G4String fileN = "currentRun";
    if(rngStatusEventsFlag)
    {
      std::ostringstream os;
      os << "run" << currentRun->GetRunID();
      fileN = os.str();
    }
    StoreRNGStatus(fileN);
  }

  runAborted             = false;
  numberOfEventProcessed = 0;
}

// --------------------------------------------------------------------
void G4WorkerRunManager::DoEventLoop(G4int n_event, const char* macroFile,
                                     G4int n_select)
{
  if(userPrimaryGeneratorAction == nullptr)
  {
    G4Exception("G4RunManager::GenerateEvent()", "Run0032", FatalException,
                "G4VUserPrimaryGeneratorAction is not defined!");
  }

  // This is the same as in the sequential case, just the for-loop indexes are
  // different
  InitializeEventLoop(n_event, macroFile, n_select);

  // Reset random number seeds queue
  while(seedsQueue.size() > 0)
  {
    seedsQueue.pop();
  }
  // for each run, worker should receive at least one set of random number
  // seeds.
  runIsSeeded = false;

  // Event loop
  eventLoopOnGoing = true;
  ///////    G4int i_event = workerContext->GetThreadId();
  G4int i_event = -1;
  nevModulo     = -1;
  currEvID      = -1;

  while(eventLoopOnGoing)
  {
    ProcessOneEvent(i_event);
    if(eventLoopOnGoing)
    {
      TerminateOneEvent();
      if(runAborted)
      {
        eventLoopOnGoing = false;
      }
      //////        else
      //////        {
      //////          i_event += workerContext->GetNumberThreads();
      //////          eventLoopOnGoing = i_event<n_event;
      //////        }
    }
  }

  TerminateEventLoop();
}

// --------------------------------------------------------------------
void G4WorkerRunManager::ProcessOneEvent(G4int i_event)
{
  currentEvent = GenerateEvent(i_event);
  if(eventLoopOnGoing)
  {
    eventManager->ProcessOneEvent(currentEvent);
    AnalyzeEvent(currentEvent);
    UpdateScoring();
    if(currentEvent->GetEventID() < n_select_msg)
      G4UImanager::GetUIpointer()->ApplyCommand(msgText);
  }
}

// --------------------------------------------------------------------
G4Event* G4WorkerRunManager::GenerateEvent(G4int i_event)
{
  G4Event* anEvent          = new G4Event(i_event);
  G4long s1                   = 0;
  G4long s2                   = 0;
  G4long s3                   = 0;
  G4bool eventHasToBeSeeded = true;
  if(G4MTRunManager::SeedOncePerCommunication() == 1 && runIsSeeded)
  {
    eventHasToBeSeeded = false;
  }

  if(i_event < 0)
  {
    G4int nevM = G4MTRunManager::GetMasterRunManager()->GetEventModulo();
    if(nevM == 1)
    {
      eventLoopOnGoing = G4MTRunManager::GetMasterRunManager()->SetUpAnEvent(
        anEvent, s1, s2, s3, eventHasToBeSeeded);
      runIsSeeded = true;
    }
    else
    {
      if(nevModulo <= 0)
      {
        G4int nevToDo = G4MTRunManager::GetMasterRunManager()->SetUpNEvents(
          anEvent, &seedsQueue, eventHasToBeSeeded);
        if(nevToDo == 0)
        {
          eventLoopOnGoing = false;
        }
        else
        {
          currEvID  = anEvent->GetEventID();
          nevModulo = nevToDo - 1;
        }
      }
      else
      {
        if(G4MTRunManager::SeedOncePerCommunication() > 0)
          eventHasToBeSeeded = false;
        anEvent->SetEventID(++currEvID);
        --nevModulo;
      }
      if(eventLoopOnGoing && eventHasToBeSeeded)
      {
        s1 = seedsQueue.front();
        seedsQueue.pop();
        s2 = seedsQueue.front();
        seedsQueue.pop();
      }
    }

    if(!eventLoopOnGoing)
    {
      delete anEvent;
      return nullptr;
    }
  }
  else if(eventHasToBeSeeded)
  {
    // Need to reseed random number generator
    G4RNGHelper* helper = G4RNGHelper::GetInstance();
    s1                  = helper->GetSeed(i_event * 2);
    s2                  = helper->GetSeed(i_event * 2 + 1);
  }

  if(eventHasToBeSeeded)
  {
    G4long seeds[3] = { s1, s2, 0 };
    G4Random::setTheSeeds(seeds, luxury);
    runIsSeeded = true;
  }

  // Read from file seed.
  // Andrea Dotti 4 November 2015
  // This is required for strong-reproducibility, in MT mode we have that each
  // thread produces, for each event a status file, we want to do that.
  // Search a random file with the format run{%d}evt{%d}.rndm

  // This is the filename base constructed from run and event
  const auto filename = [&] {
    std::ostringstream os;
    os << "run" << currentRun->GetRunID() << "evt" << anEvent->GetEventID();
    return os.str();
  };

  G4bool RNGstatusReadFromFile = false;
  if(readStatusFromFile)
  {
    // Build full path of RNG status file for this event
    std::ostringstream os;
    os << filename() << ".rndm";
    const G4String& randomStatusFile = os.str();
    std::ifstream ifile(randomStatusFile.c_str());
    if(ifile)
    {  // File valid and readable
      RNGstatusReadFromFile = true;
      G4Random::restoreEngineStatus(randomStatusFile.c_str());
    }
  }

  if(storeRandomNumberStatusToG4Event == 1 ||
     storeRandomNumberStatusToG4Event == 3)
  {
    std::ostringstream oss;
    G4Random::saveFullState(oss);
    randomNumberStatusForThisEvent = oss.str();
    anEvent->SetRandomNumberStatus(randomNumberStatusForThisEvent);
  }

  if(storeRandomNumberStatus && !RNGstatusReadFromFile)
  {  // If reading from file, avoid to rewrite the same
    G4String fileN = "currentEvent";
    if(rngStatusEventsFlag)
    {
      fileN = filename();
    }
    StoreRNGStatus(fileN);
  }

  if(printModulo > 0 && anEvent->GetEventID() % printModulo == 0)
  {
    G4cout << "--> Event " << anEvent->GetEventID() << " starts";
    if(eventHasToBeSeeded)
    {
      G4cout << " with initial seeds (" << s1 << "," << s2 << ")";
    }
    G4cout << "." << G4endl;
  }
  userPrimaryGeneratorAction->GeneratePrimaries(anEvent);
  return anEvent;
}

// --------------------------------------------------------------------
void G4WorkerRunManager::MergePartialResults()
{
  // Merge partial results into global run
  G4MTRunManager* mtRM  = G4MTRunManager::GetMasterRunManager();
  G4ScoringManager* ScM = G4ScoringManager::GetScoringManagerIfExist();
  if(ScM != nullptr)
    mtRM->MergeScores(ScM);
  mtRM->MergeRun(currentRun);
}

// --------------------------------------------------------------------
void G4WorkerRunManager::RunTermination()
{
  if(!fakeRun)
  {
#if defined(GEANT4_USE_TIMEMORY)
    workerRunProfiler.reset();
#endif
    MergePartialResults();

    // Call a user hook: note this is before the next barrier
    // so threads execute this method asyncrhonouzly
    //(TerminateRun allows for synch via G4RunAction::EndOfRun)
    const G4UserWorkerInitialization* uwi =
      G4MTRunManager::GetMasterRunManager()->GetUserWorkerInitialization();
    if(uwi != nullptr)
      uwi->WorkerRunEnd();
  }

  G4RunManager::RunTermination();
  // Signal this thread has finished envent-loop.
  // Note this will return only whan all threads reach this point
  G4MTRunManager::GetMasterRunManager()->ThisWorkerEndEventLoop();
}

// --------------------------------------------------------------------
void G4WorkerRunManager::TerminateEventLoop()
{
  if(verboseLevel > 0 && !fakeRun)
  {
    timer->Stop();
    G4cout << "Thread-local run terminated." << G4endl;
    G4cout << "Run Summary" << G4endl;
    if(runAborted)
    {
      G4cout << "  Run Aborted after " << numberOfEventProcessed
             << " events processed." << G4endl;
    }
    else
    {
      G4cout << "  Number of events processed : " << numberOfEventProcessed
             << G4endl;
    }
    G4cout << "  " << *timer << G4endl;
  }
}

// --------------------------------------------------------------------
namespace
{
  G4Mutex ConstructScoringWorldsMutex = G4MUTEX_INITIALIZER;
}

// --------------------------------------------------------------------
void G4WorkerRunManager::ConstructScoringWorlds()
{
  using MeshShape = G4VScoringMesh::MeshShape;

  // Return if unnecessary
  G4ScoringManager* ScM = G4ScoringManager::GetScoringManagerIfExist();
  if(ScM == nullptr)
    return;
  G4int nPar = (G4int)ScM->GetNumberOfMesh();
  if(nPar < 1)
    return;

  // Update thread-local G4TransportationManager of all the world volumes
  kernel->WorkerUpdateWorldVolume();

  G4ScoringManager* masterScM = G4MTRunManager::GetMasterScoringManager();

  auto particleIterator = G4ParticleTable::GetParticleTable()->GetIterator();

  for(G4int iw = 0; iw < nPar; ++iw)
  {
    G4VScoringMesh* mesh = ScM->GetMesh(iw);
    if(fGeometryHasBeenDestroyed)
      mesh->GeometryHasBeenDestroyed();
    G4VPhysicalVolume* pWorld = nullptr;
    if(mesh->GetShape() != MeshShape::realWorldLogVol)
    {
      pWorld = G4TransportationManager::GetTransportationManager()
             ->IsWorldExisting(ScM->GetWorldName(iw));
      if(pWorld == nullptr)
      {
        G4ExceptionDescription ed;
        ed << "Mesh name <" << ScM->GetWorldName(iw)
           << "> is not found in the master thread.";
        G4Exception("G4WorkerRunManager::ConstructScoringWorlds()", "RUN79001",
                    FatalException, ed);
      }
    }
    if(!(mesh->GetMeshElementLogical()))
    {
      G4AutoLock l(&ConstructScoringWorldsMutex);
      G4VScoringMesh* masterMesh = masterScM->GetMesh(iw);
      mesh->SetMeshElementLogical(masterMesh->GetMeshElementLogical());
      l.unlock();

      if(mesh->GetShape() != MeshShape::realWorldLogVol)
      {
        G4ParallelWorldProcess* theParallelWorldProcess =
                                mesh->GetParallelWorldProcess();
        if(theParallelWorldProcess != nullptr)
        {
          theParallelWorldProcess->SetParallelWorld(ScM->GetWorldName(iw));
        }
        else
        {
          theParallelWorldProcess =
            new G4ParallelWorldProcess(ScM->GetWorldName(iw));
          mesh->SetParallelWorldProcess(theParallelWorldProcess);
          theParallelWorldProcess->SetParallelWorld(ScM->GetWorldName(iw));

          particleIterator->reset();
          while((*particleIterator)())
          {
            G4ParticleDefinition* particle = particleIterator->value();
            G4ProcessManager* pmanager     = particle->GetProcessManager();
            if(pmanager != nullptr)
            {
              pmanager->AddProcess(theParallelWorldProcess);
              if(theParallelWorldProcess->IsAtRestRequired(particle))
              {
                pmanager->SetProcessOrdering(theParallelWorldProcess, idxAtRest,
                                             9900);
              }
              pmanager->SetProcessOrderingToSecond(theParallelWorldProcess,
                                                   idxAlongStep);
              pmanager->SetProcessOrdering(theParallelWorldProcess, idxPostStep,
                                           9900);
            }  // if(pmanager)
          }    // while
        }
        theParallelWorldProcess->SetLayeredMaterialFlag(mesh->LayeredMassFlg());
      }
    }
    mesh->WorkerConstruct(pWorld);
  }
}

// --------------------------------------------------------------------
void G4WorkerRunManager::SetUserInitialization(G4UserWorkerInitialization*)
{
  G4Exception(
    "G4RunManager::SetUserInitialization(G4UserWorkerInitialization*)",
    "Run0118", FatalException,
    "This method should be used only with an instance of G4MTRunManager");
}

// --------------------------------------------------------------------
void G4WorkerRunManager::SetUserInitialization(
  G4UserWorkerThreadInitialization*)
{
  G4Exception(
    "G4RunManager::SetUserInitialization(G4UserWorkerThreadInitialization*)",
    "Run0119", FatalException,
    "This method should be used only with an instance of G4MTRunManager");
}

// --------------------------------------------------------------------
void G4WorkerRunManager::SetUserInitialization(G4VUserActionInitialization*)
{
  G4Exception(
    "G4RunManager::SetUserInitialization(G4VUserActionInitialization*)",
    "Run0120", FatalException,
    "This method should be used only with an instance of G4MTRunManager");
}

// --------------------------------------------------------------------
void G4WorkerRunManager::SetUserInitialization(G4VUserDetectorConstruction*)
{
  G4Exception(
    "G4RunManager::SetUserInitialization(G4VUserDetectorConstruction*)",
    "Run0121", FatalException,
    "This method should be used only with an instance of G4MTRunManager");
}

// --------------------------------------------------------------------
void G4WorkerRunManager::SetUserInitialization(G4VUserPhysicsList* pl)
{
  pl->InitializeWorker();
  G4RunManager::SetUserInitialization(pl);
}

// --------------------------------------------------------------------
void G4WorkerRunManager::SetUserAction(G4UserRunAction* userAction)
{
  G4RunManager::SetUserAction(userAction);
  if(userAction != nullptr)
    userAction->SetMaster(false);
}

// --------------------------------------------------------------------
void G4WorkerRunManager::SetupDefaultRNGEngine()
{
  const CLHEP::HepRandomEngine* mrnge =
    G4MTRunManager::GetMasterRunManager()->getMasterRandomEngine();
  const G4UserWorkerThreadInitialization* uwti =
    G4MTRunManager::GetMasterRunManager()->GetUserWorkerThreadInitialization();
  uwti->SetupRNGEngine(mrnge);
}

// Forward calls (avoid GCC compilation warnings)

// --------------------------------------------------------------------
void G4WorkerRunManager::SetUserAction(G4UserEventAction* ua)
{
  G4RunManager::SetUserAction(ua);
}

// --------------------------------------------------------------------
void G4WorkerRunManager::SetUserAction(G4VUserPrimaryGeneratorAction* ua)
{
  G4RunManager::SetUserAction(ua);
}

// --------------------------------------------------------------------
void G4WorkerRunManager::SetUserAction(G4UserStackingAction* ua)
{
  G4RunManager::SetUserAction(ua);
}

// --------------------------------------------------------------------
void G4WorkerRunManager::SetUserAction(G4UserTrackingAction* ua)
{
  G4RunManager::SetUserAction(ua);
}

// --------------------------------------------------------------------
void G4WorkerRunManager::SetUserAction(G4UserSteppingAction* ua)
{
  G4RunManager::SetUserAction(ua);
}

// --------------------------------------------------------------------
void G4WorkerRunManager::StoreRNGStatus(const G4String& fn)
{
  std::ostringstream os;
  os << randomNumberStatusDir << "G4Worker" << workerContext->GetThreadId()
     << "_" << fn << ".rndm";
  G4Random::saveEngineStatus(os.str().c_str());
}

// --------------------------------------------------------------------
void G4WorkerRunManager::rndmSaveThisRun()
{
  G4int runNumber = 0;
  if(currentRun != nullptr)
    runNumber = currentRun->GetRunID();
  if(!storeRandomNumberStatus)
  {
    G4cerr << "Warning from G4RunManager::rndmSaveThisRun():"
           << " Random number status was not stored prior to this run."
           << G4endl << "/random/setSavingFlag command must be issued. "
           << "Command ignored." << G4endl;
    return;
  }

  std::ostringstream oos;
  oos << "G4Worker" << workerContext->GetThreadId() << "_"
      << "currentRun.rndm"
      << "\0";
  G4String fileIn = randomNumberStatusDir + oos.str();

  std::ostringstream os;
  os << "run" << runNumber << ".rndm" << '\0';
  G4String fileOut = randomNumberStatusDir + os.str();

#ifdef WIN32
  G4String copCmd = "/control/shell copy " + fileIn + " " + fileOut;
#else
  G4String copCmd = "/control/shell cp " + fileIn + " " + fileOut;
#endif
  G4UImanager::GetUIpointer()->ApplyCommand(copCmd);
  if(verboseLevel > 0)
  {
    G4cout << fileIn << " is copied to " << fileOut << G4endl;
  }
}

// --------------------------------------------------------------------
void G4WorkerRunManager::rndmSaveThisEvent()
{
  if(currentEvent == nullptr)
  {
    G4cerr << "Warning from G4RunManager::rndmSaveThisEvent():"
           << " there is no currentEvent available." << G4endl
           << "Command ignored." << G4endl;
    return;
  }

  if(!storeRandomNumberStatus)
  {
    G4cerr << "Warning from G4RunManager::rndmSaveThisEvent():"
           << " Random number engine status is not available." << G4endl
           << "/random/setSavingFlag command must be issued "
           << "prior to the start of the run. Command ignored." << G4endl;
    return;
  }

  std::ostringstream oos;
  oos << "G4Worker" << workerContext->GetThreadId() << "_"
      << "currentEvent.rndm"
      << "\0";
  G4String fileIn = randomNumberStatusDir + oos.str();

  std::ostringstream os;
  os << "run" << currentRun->GetRunID() << "evt" << currentEvent->GetEventID()
     << ".rndm" << '\0';
  G4String fileOut = randomNumberStatusDir + os.str();

#ifdef WIN32
  G4String copCmd = "/control/shell copy " + fileIn + " " + fileOut;
#else
  G4String copCmd = "/control/shell cp " + fileIn + " " + fileOut;
#endif
  G4UImanager::GetUIpointer()->ApplyCommand(copCmd);
  if(verboseLevel > 0)
  {
    G4cout << fileIn << " is copied to " << fileOut << G4endl;
  }
}

// --------------------------------------------------------------------
void G4WorkerRunManager::DoWork()
{
  G4MTRunManager* mrm = G4MTRunManager::GetMasterRunManager();
  G4MTRunManager::WorkerActionRequest nextAction =
    mrm->ThisWorkerWaitForNextAction();
  while(nextAction != G4MTRunManager::WorkerActionRequest::ENDWORKER)
  {
    if(nextAction ==
       G4MTRunManager::WorkerActionRequest::NEXTITERATION)  // start the next
                                                            // run
    {
      // The following code deals with changing materials between runs
      static G4ThreadLocal G4bool skipInitialization = true;
      if(skipInitialization)
      {
        // re-initialization is not necessary for the first run
        skipInitialization = false;
      }
      else
      {
        //        ReinitializeGeometry();
        workerContext->UpdateGeometryAndPhysicsVectorFromMaster();
      }

      // Execute UI commands stored in the master UI manager
      std::vector<G4String> cmds = mrm->GetCommandStack();
      G4UImanager* uimgr         = G4UImanager::GetUIpointer(); // TLS instance
      for(auto it = cmds.cbegin(); it != cmds.cend(); ++it)
      {
        uimgr->ApplyCommand(*it);
      }
      // Start this run
      G4int numevents    = mrm->GetNumberOfEventsToBeProcessed();
      G4String macroFile = mrm->GetSelectMacro();
      G4int numSelect    = mrm->GetNumberOfSelectEvents();
      if(macroFile == "" || macroFile == " ")
      {
        this->BeamOn(numevents);
      }
      else
      {
        this->BeamOn(numevents, macroFile, numSelect);
      }
    }
    else if(nextAction == G4MTRunManager::WorkerActionRequest::PROCESSUI)
    {
      std::vector<G4String> cmds = mrm->GetCommandStack();
      G4UImanager* uimgr         = G4UImanager::GetUIpointer(); // TLS instance
      for(auto it = cmds.cbegin(); it != cmds.cend(); ++it)
      {
        uimgr->ApplyCommand(*it);
      }
      mrm->ThisWorkerProcessCommandsStackDone();
    }
    else
    {
      G4ExceptionDescription d;
      d << "Cannot continue, this worker has been requested an unknown action: "
        << static_cast<
             std::underlying_type<G4MTRunManager::WorkerActionRequest>::type>(
             nextAction);
      G4Exception("G4WorkerRunManager::DoWork", "Run0104", FatalException, d);
    }

    // Now wait for master thread to signal new action to be performed
    nextAction = mrm->ThisWorkerWaitForNextAction();
  }  // No more actions to perform

  return;
}
