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

#include "G4WorkerSubEvtRunManager.hh"

#include "G4AutoLock.hh"
#include "G4MTRunManager.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4ParallelWorldProcessStore.hh"
#include "G4RNGHelper.hh"
#include "G4Run.hh"
#include "G4SDManager.hh"
#include "G4ScoringManager.hh"
#include "G4SubEvtRunManager.hh"
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
#include "G4WorkerTaskRunManagerKernel.hh"
#include "G4WorkerThread.hh"

#include <fstream>
#include <sstream>

//============================================================================//

G4WorkerSubEvtRunManager* G4WorkerSubEvtRunManager::GetWorkerRunManager()
{
  return static_cast<G4WorkerSubEvtRunManager*>(G4RunManager::GetRunManager());
}

//============================================================================//

G4WorkerSubEvtRunManagerKernel* G4WorkerSubEvtRunManager::GetWorkerRunManagerKernel()
{
  return static_cast<G4WorkerSubEvtRunManagerKernel*>(GetWorkerRunManager()->kernel);
}

//============================================================================//

G4WorkerSubEvtRunManager::G4WorkerSubEvtRunManager(G4int seType)
{
  runManagerType = subEventWorkerRM; 
  SetSubEventType(seType);
}

void G4WorkerSubEvtRunManager::RunInitialization()
{
#ifdef G4MULTITHREADED
  if (!visIsSetUp) {
    G4VVisManager* pVVis = G4VVisManager::GetConcreteInstance();
    if (pVVis != nullptr) {
      pVVis->SetUpForAThread();
      visIsSetUp = true;
    }
  }
#endif
  runIsSeeded = false;

  if (!(kernel->RunInitialization(fakeRun))) return;

  // Signal this thread can start event loop.
  // Note this will return only when all threads reach this point
  G4MTRunManager::GetMasterRunManager()->ThisWorkerReady();
  if (fakeRun) return;

  const G4UserWorkerInitialization* uwi =
    G4MTRunManager::GetMasterRunManager()->GetUserWorkerInitialization();

  CleanUpPreviousEvents();

  delete currentRun;

  currentRun = nullptr;

  if (IfGeometryHasBeenDestroyed()) G4ParallelWorldProcessStore::GetInstance()->UpdateWorlds();

  // Call a user hook: this is guaranteed all threads are "synchronized"
  if (uwi != nullptr) uwi->WorkerRunStart();

  if (userRunAction != nullptr) currentRun = userRunAction->GenerateRun();
  if (currentRun == nullptr) currentRun = new G4Run();

  currentRun->SetRunID(runIDCounter);
  G4TaskRunManager* mrm = G4TaskRunManager::GetMasterRunManager();
  numberOfEventToBeProcessed = mrm->GetNumberOfEventsToBeProcessed();
  currentRun->SetNumberOfEventToBeProcessed(numberOfEventToBeProcessed);

  currentRun->SetDCtable(DCtable);
  G4SDManager* fSDM = G4SDManager::GetSDMpointerIfExist();
  if (fSDM != nullptr) {
    currentRun->SetHCtable(fSDM->GetHCtable());
  }

  if (G4VScoreNtupleWriter::Instance() != nullptr) {
    auto hce = (fSDM != nullptr) ? fSDM->PrepareNewEvent() : nullptr;
    isScoreNtupleWriter = G4VScoreNtupleWriter::Instance()->Book(hce);
    delete hce;
  }

  std::ostringstream oss;
  G4Random::saveFullState(oss);
  randomNumberStatusForThisRun = oss.str();
  currentRun->SetRandomNumberStatus(randomNumberStatusForThisRun);

  for (G4int i_prev = 0; i_prev < n_perviousEventsToBeStored; ++i_prev)
    previousEvents->push_back(nullptr);

  if (printModulo > 0 || verboseLevel > 0) {
    G4cout << "### Run " << currentRun->GetRunID() << " starts on worker thread "
           << G4Threading::G4GetThreadId() << "." << G4endl;
  }

  if (userRunAction != nullptr) userRunAction->BeginOfRunAction(currentRun);

  if (isScoreNtupleWriter) {
    G4VScoreNtupleWriter::Instance()->OpenFile();
  }

  if (storeRandomNumberStatus) {
    G4String fileN = "currentRun";
    if (rngStatusEventsFlag) {
      std::ostringstream os;
      os << "run" << currentRun->GetRunID();
      fileN = os.str();
    }
    StoreRNGStatus(fileN);
  }

  runAborted = false;
  numberOfEventProcessed = 0;
  if(verboseLevel > 0) timer->Start();
}

//============================================================================//

void G4WorkerSubEvtRunManager::DoEventLoop(G4int /*n_event*/, const char* /*macroFile*/,
                                           G4int /*n_select*/)
{
  // This method is not used in worker sub-event mode
  G4Exception("G4WorkerSubEvtRunManager::DoEventLoop()","SubEvtXXX001",FatalException,
    "This method is not used in the worker thread of sub-event parallel mode");
}

//============================================================================//

void G4WorkerSubEvtRunManager::ProcessOneEvent(G4int /*i_event*/)
{
  // This method is not used in worker sub-event mode
  G4Exception("G4WorkerSubEvtRunManager::ProcessOneEvent()","SubEvtXXX002",FatalException,
    "This method is not used in the worker thread of sub-event parallel mode");
}

//============================================================================//

G4Event* G4WorkerSubEvtRunManager::GenerateEvent(G4int /*i_event*/)
{
  // This method is not used in worker sub-event mode
  G4Exception("G4WorkerSubEvtRunManager::GenerateEvent()","SubEvtXXX003",FatalException,
    "This method is not used in the worker thread of sub-event parallel mode");
  return nullptr;
}

//============================================================================//

void G4WorkerSubEvtRunManager::RunTermination()
{
  if (!fakeRun && (currentRun != nullptr)) {
    MergePartialResults(true);

    // Call a user hook: note this is before the next barrier
    // so threads execute this method asyncrhonouzly
    //(TerminateRun allows for synch via G4RunAction::EndOfRun)
    const G4UserWorkerInitialization* uwi =
      G4MTRunManager::GetMasterRunManager()->GetUserWorkerInitialization();
    if (uwi != nullptr) uwi->WorkerRunEnd();
  }

  if (currentRun != nullptr) {
    G4RunManager::RunTermination();
  }
  // Signal this thread has finished envent-loop.
  // Note this will return only whan all threads reach this point
  G4MTRunManager::GetMasterRunManager()->ThisWorkerEndEventLoop();
}

//============================================================================//

void G4WorkerSubEvtRunManager::TerminateEventLoop()
{
  if (verboseLevel > 0 && !fakeRun) {
    timer->Stop();
    // prefix with thread # info due to how TBB calls this function
    G4String prefix = "[thread " + std::to_string(workerContext->GetThreadId()) + "] ";
    G4cout << prefix << "Thread-local run terminated." << G4endl;
    G4cout << prefix << "Run Summary" << G4endl;
    if (runAborted)
      G4cout << prefix << "  Run Aborted after " << numberOfEventProcessed << " sub-events processed."
             << G4endl;
    else
      G4cout << prefix << "  Number of sub-events processed : " << numberOfEventProcessed << G4endl;
    G4cout << prefix << "  " << *timer << G4endl;
  }
}

//============================================================================//

void G4WorkerSubEvtRunManager::SetupDefaultRNGEngine()
{
  const CLHEP::HepRandomEngine* mrnge =
    G4MTRunManager::GetMasterRunManager()->getMasterRandomEngine();
  assert(mrnge);  // Master has created RNG
  const G4UserWorkerThreadInitialization* uwti =
    G4MTRunManager::GetMasterRunManager()->GetUserWorkerThreadInitialization();
  uwti->SetupRNGEngine(mrnge);
}

//============================================================================//

void G4WorkerSubEvtRunManager::StoreRNGStatus(const G4String& fn)
{
  std::ostringstream os;
  os << randomNumberStatusDir << "G4Worker" << workerContext->GetThreadId() << "_" << fn << ".rndm";
  G4Random::saveEngineStatus(os.str().c_str());
}

//============================================================================//

void G4WorkerSubEvtRunManager::ProcessUI()
{
  G4TaskRunManager* mrm = G4TaskRunManager::GetMasterRunManager();
  if (mrm == nullptr) return;

  //------------------------------------------------------------------------//
  // Check UI commands not already processed
  auto command_stack = mrm->GetCommandStack();
  bool matching = (command_stack.size() == processedCommandStack.size());
  if (matching) {
    for (uintmax_t i = 0; i < command_stack.size(); ++i)
      if (processedCommandStack.at(i) != command_stack.at(i)) {
        matching = false;
        break;
      }
  }

  //------------------------------------------------------------------------//
  // Execute UI commands stored in the master UI manager
  if (!matching) {
    for (const auto& itr : command_stack)
      G4UImanager::GetUIpointer()->ApplyCommand(itr);
    processedCommandStack = std::move(command_stack);
  }
}

//============================================================================//

void G4WorkerSubEvtRunManager::DoCleanup()
{
  // Nothing to do for a run 

  //CleanUpPreviousEvents();
  //
  //delete currentRun;
  //currentRun = nullptr;
}

//============================================================================//

void G4WorkerSubEvtRunManager::DoWork()
{
  if(verboseLevel>1) {
    G4cout << "G4WorkerSubEvtRunManager::DoWork() starts.........." << G4endl;
  }

  G4SubEvtRunManager* mrm = G4SubEvtRunManager::GetMasterRunManager();
  G4bool newRun = false;
  const G4Run* run = mrm->GetCurrentRun();
  G4ThreadLocalStatic G4int runId = -1;
  if ((run != nullptr) && run->GetRunID() != runId) {
    runId = run->GetRunID();
    newRun = true;
    if (runId > 0) { ProcessUI(); }
  }

  G4bool reseedRequired = false;
  if (newRun) {
    G4bool cond = ConfirmBeamOnCondition();
    if (cond) {
      ConstructScoringWorlds();
      RunInitialization();
    }
    reseedRequired = true;
  }

  assert(workerContext != nullptr);
  workerContext->UpdateGeometryAndPhysicsVectorFromMaster();

  eventManager->UseSubEventParallelism(true);

  G4bool needMoreWork = true;
  while(needMoreWork)
  {
    G4bool notReady = false;
    G4long s1, s2, s3;
    auto subEv = mrm->GetSubEvent(fSubEventType, notReady, s1, s2, s3, reseedRequired);
    if(subEv==nullptr && notReady)
    {
      // Master is not yet ready for tasking a sub-event. 
      // Wait 1 second and retry.
      G4THREADSLEEP(1);
    }
    else if(subEv==nullptr)
    { 
      // No more sub-event to process
      // Report the results of previous sub-events if any
      G4Event* remainingE = nullptr;
      do
      {
        remainingE = eventManager->RetrieveCompletedSubEvent();
        if(remainingE)
        {
          mrm->SubEventFinished(remainingE->GetSubEvent(),remainingE);
          delete remainingE;
        }
      }
      while(remainingE);

      // Check if no more sub-event in event manager
      if(eventManager->GetNumberOfRemainingSubEvents()==0)
      { needMoreWork = false; }
      else
      { // Wait 1 second and revisit.
        G4THREADSLEEP(1);
        if(verboseLevel>1) {
          G4cout << "G4WorkerSubEvtRunManager::DoWork() - " 
                 << eventManager->GetNumberOfRemainingSubEvents()
                 << " sub-events are still incomplete in the event manager."<< G4endl;
        }
      }
    }
    else
    {
      // Let's work for this sub-event.
      if(reseedRequired)
      {
        G4long seeds[3] = {s1, s2, s3};
        G4Random::setTheSeeds(seeds, -1);
        reseedRequired = false;
      }

      // create a G4Event object for this sub-event. This G4Event object will contain output
      // to be merged into the master event.
      auto masterEvent = subEv->GetEvent();
      G4Event* ev = new G4Event(masterEvent->GetEventID());
      ev->FlagAsSubEvent(masterEvent,fSubEventType,subEv);
      ++numberOfEventProcessed;

      // Create a G4TrackVector as the input
      G4TrackVector* tv = new G4TrackVector();
      for(auto& stackedTrack : *subEv)
      {
        // tracks (and trajectories) stored in G4SubEvent object belong to the master thread
        // and thus they must not be deleted by the worker thread. They must be cloned.
        G4Track* tr = new G4Track();
        tr->CopyTrackInfo(*(stackedTrack.GetTrack()),false);
        tv->push_back(tr);
      }

      // Process this sub-event
      currentEvent = ev;
      eventManager->ProcessOneEvent(tv,ev);

      // We don't need following two lines, as they are taken care by the master
      //////AnalyzeEvent(ev);
      //////UpdateScoring();

      // Report the results to the master if sub-event is completed
      if(subEv->IsCompleted())
      {
        mrm->SubEventFinished(subEv,ev); 
        delete ev;
      }
    
      // Report the results of previous sub-events if any
      G4Event* remainingE = nullptr;
      do 
      {
        remainingE = eventManager->RetrieveCompletedSubEvent(); 
        if(remainingE)
        {
          mrm->SubEventFinished(remainingE->GetSubEvent(),remainingE);
          delete remainingE;
        }
      }
      while(remainingE);

      if(verboseLevel>2)
      {
        G4cout << "G4WorkerSubEvtRunManager::DoWork() - " 
               << eventManager->GetNumberOfRemainingSubEvents()
               << " sub-events are still incomplete in the event manager."<< G4endl;
      }

      // clean up
      delete tv;
    }
  }

  if(verboseLevel>1) {
    G4cout << "G4WorkerSubEvtRunManager::DoWork() completed.........." << G4endl;
  }
    
}

void G4WorkerSubEvtRunManager::SetSubEventType(G4int ty)
{
  auto* mrm = G4SubEvtRunManager::GetMasterRunManager();
  mrm->RegisterSubEvtWorker(this,ty);
  fSubEventType = ty;
}

//============================================================================//

void G4WorkerSubEvtRunManager::SetUserInitialization(G4UserWorkerInitialization*)
{
  G4Exception("G4WorkerSubEvtRunManager::SetUserInitialization(G4UserWorkerInitialization*)", "RunSE0118",
              FatalException, "This method should be used only with an instance of the master thread");
}

// --------------------------------------------------------------------
void G4WorkerSubEvtRunManager::SetUserInitialization(G4UserWorkerThreadInitialization*)
{
  G4Exception("G4WorkerSubEvtRunManager::SetUserInitialization(G4UserWorkerThreadInitialization*)", "RunSE0119",
              FatalException, "This method should be used only with an instance of the master thread");
}

// --------------------------------------------------------------------
void G4WorkerSubEvtRunManager::SetUserInitialization(G4VUserActionInitialization*)
{
  G4Exception("G4WorkerSubEvtRunManager::SetUserInitialization(G4VUserActionInitialization*)", "RunSE0120",
              FatalException, "This method should be used only with an instance of the master thread");
}

// --------------------------------------------------------------------
void G4WorkerSubEvtRunManager::SetUserInitialization(G4VUserDetectorConstruction*)
{
  G4Exception("G4WorkerSubEvtRunManager::SetUserInitialization(G4VUserDetectorConstruction*)", "RunSE0121",
              FatalException, "This method should be used only with an instance of the master thread");
}

// --------------------------------------------------------------------
void G4WorkerSubEvtRunManager::SetUserInitialization(G4VUserPhysicsList* pl)
{
  pl->InitializeWorker();
  G4RunManager::SetUserInitialization(pl);
}

// --------------------------------------------------------------------
void G4WorkerSubEvtRunManager::SetUserAction(G4UserRunAction*)
{
  G4Exception("G4WorkerSubEvtRunManager::SetUserAction(G4UserRunAction*)", "RunSE0221",
              FatalException, "This method should be used only with an instance of the master thread");
}

// Forward calls (avoid GCC compilation warnings)

// --------------------------------------------------------------------
void G4WorkerSubEvtRunManager::SetUserAction(G4UserEventAction* ua)
{
  G4RunManager::SetUserAction(ua);
}

// --------------------------------------------------------------------
void G4WorkerSubEvtRunManager::SetUserAction(G4VUserPrimaryGeneratorAction*)
{
  G4Exception("G4WorkerSubEvtRunManager::SetUserAction(G4VUserPrimaryGeneratorAction*)", "RunSE0223",
              FatalException, "This method should be used only with an instance of the master thread");
}

// --------------------------------------------------------------------
void G4WorkerSubEvtRunManager::SetUserAction(G4UserStackingAction* ua)
{
  G4RunManager::SetUserAction(ua);
}

// --------------------------------------------------------------------
void G4WorkerSubEvtRunManager::SetUserAction(G4UserTrackingAction* ua)
{
  G4RunManager::SetUserAction(ua);
}

// --------------------------------------------------------------------
void G4WorkerSubEvtRunManager::SetUserAction(G4UserSteppingAction* ua)
{
  G4RunManager::SetUserAction(ua);
}




