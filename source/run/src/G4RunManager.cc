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
// $Id: G4RunManager.cc 103291 2017-03-24 14:00:49Z gcosmo $
//
// 

// On Sun, to prevent conflict with ObjectSpace, G4Timer.hh has to be
// loaded *before* globals.hh...
#include "G4Timer.hh"

#include "G4RunManager.hh"
#include "G4RunManagerKernel.hh"
#include "G4MTRunManagerKernel.hh"
#include "G4WorkerRunManagerKernel.hh"

#include "G4StateManager.hh"
#include "G4ApplicationState.hh"
#include "Randomize.hh"
#include "G4Run.hh"
#include "G4RunMessenger.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4UserWorkerInitialization.hh"
#include "G4UserWorkerThreadInitialization.hh"
#include "G4UserRunAction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VPersistencyManager.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessTable.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Material.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ParallelWorldProcessStore.hh"
#include "G4ios.hh"
#include <sstream>


using namespace CLHEP;

G4ThreadLocal G4RunManager* G4RunManager::fRunManager = 0;

G4bool G4RunManager::fGeometryHasBeenDestroyed = false;
G4bool G4RunManager::IfGeometryHasBeenDestroyed() { return fGeometryHasBeenDestroyed; }

//The following lines are needed since G4VUserPhysicsList
//uses a #define theParticleIterator
#ifdef theParticleIterator
#undef theParticleIterator
#endif

G4RunManager* G4RunManager::GetRunManager()
{ return fRunManager; }

G4RunManager::G4RunManager()
:userDetector(0),physicsList(0),
 userActionInitialization(0),userWorkerInitialization(0),
 userWorkerThreadInitialization(0),
 userRunAction(0),userPrimaryGeneratorAction(0),userEventAction(0),
 userStackingAction(0),userTrackingAction(0),userSteppingAction(0),
 geometryInitialized(false),physicsInitialized(false),
 runAborted(false),initializedAtLeastOnce(false),
 geometryToBeOptimized(true),runIDCounter(0),
 verboseLevel(0),printModulo(-1),DCtable(0),
 currentRun(0),currentEvent(0),n_perviousEventsToBeStored(0),
 numberOfEventToBeProcessed(0),storeRandomNumberStatus(false),
 storeRandomNumberStatusToG4Event(0),rngStatusEventsFlag(false),
 currentWorld(0),nParallelWorlds(0),msgText(" "),n_select_msg(-1),
 numberOfEventProcessed(0),selectMacro(""),fakeRun(false)
{
  if(fRunManager)
  {
    G4Exception("G4RunManager::G4RunManager()", "Run0031",
                FatalException, "G4RunManager constructed twice.");
  }
  fRunManager = this;
      
  kernel = new G4RunManagerKernel();
  eventManager = kernel->GetEventManager();

  timer = new G4Timer();
  runMessenger = new G4RunMessenger(this);
  previousEvents = new std::list<G4Event*>;
  G4ParticleTable::GetParticleTable()->CreateMessenger();
  G4ProcessTable::GetProcessTable()->CreateMessenger();
  randomNumberStatusDir = "./";
  std::ostringstream oss;
  G4Random::saveFullState(oss);
  randomNumberStatusForThisRun = oss.str();
  randomNumberStatusForThisEvent = oss.str();
  runManagerType = sequentialRM;
}

G4RunManager::G4RunManager( RMType rmType )
:userDetector(0),physicsList(0),
 userActionInitialization(0),userWorkerInitialization(0),
 userWorkerThreadInitialization(0),
 userRunAction(0),userPrimaryGeneratorAction(0),userEventAction(0),
 userStackingAction(0),userTrackingAction(0),userSteppingAction(0),
 geometryInitialized(false),physicsInitialized(false),
 runAborted(false),initializedAtLeastOnce(false),
 geometryToBeOptimized(true),runIDCounter(0),
 verboseLevel(0),printModulo(-1),DCtable(0),
 currentRun(0),currentEvent(0),n_perviousEventsToBeStored(0),
 numberOfEventToBeProcessed(0),storeRandomNumberStatus(false),
 storeRandomNumberStatusToG4Event(0),rngStatusEventsFlag(false),
 currentWorld(0),nParallelWorlds(0),msgText(" "),n_select_msg(-1),
 numberOfEventProcessed(0),selectMacro(""),fakeRun(false)
{
  //This version of the constructor should never be called in sequential mode!
#ifndef G4MULTITHREADED
  G4ExceptionDescription msg;
  msg<<"Geant4 code is compiled without multi-threading support (-DG4MULTITHREADED is set to off).";
  msg<<" This type of RunManager can only be used in mult-threaded applications.";
  G4Exception("G4RunManager::G4RunManager(G4bool)","Run0107",FatalException,msg);
#endif
    
  if(fRunManager)
  {
    G4Exception("G4RunManager::G4RunManager()", "Run0031",
                    FatalException, "G4RunManager constructed twice.");
    return;
  }
  fRunManager = this;
    
  switch(rmType)
  {
   case masterRM:
    kernel = new G4MTRunManagerKernel();
    break;
   case workerRM:
    kernel = new G4WorkerRunManagerKernel();
    break;
   default:
    G4ExceptionDescription msgx;
    msgx<<" This type of RunManager can only be used in mult-threaded applications.";
    G4Exception("G4RunManager::G4RunManager(G4bool)","Run0108",FatalException,msgx);
    return;
  }
  runManagerType = rmType;

  eventManager = kernel->GetEventManager();
        
  timer = new G4Timer();
  runMessenger = new G4RunMessenger(this);
  previousEvents = new std::list<G4Event*>;
  G4ParticleTable::GetParticleTable()->CreateMessenger();
  G4ProcessTable::GetProcessTable()->CreateMessenger();
  randomNumberStatusDir = "./";
  std::ostringstream oss;
  G4Random::saveFullState(oss);
  randomNumberStatusForThisRun = oss.str();
  randomNumberStatusForThisEvent = oss.str();
}

G4RunManager::~G4RunManager()
{
  G4StateManager* pStateManager = G4StateManager::GetStateManager();
  // set the application state to the quite state
  if(pStateManager->GetCurrentState()!=G4State_Quit)
  {
     if(verboseLevel>0) G4cout << "G4 kernel has come to Quit state." << G4endl;
     pStateManager->SetNewState(G4State_Quit);
  }

  CleanUpPreviousEvents();
  if(currentRun) delete currentRun;
  delete timer;
  delete runMessenger;
  G4ParticleTable::GetParticleTable()->DeleteMessenger();
  G4ProcessTable::GetProcessTable()->DeleteMessenger();
  delete previousEvents;

  //The following will work for all RunManager types
  //if derived class does the correct thing in derived
  //destructor that is set to zero pointers of
  //user initialization objects for which does not have
  //ownership
  DeleteUserInitializations();
  if(userRunAction)
  {
    delete userRunAction;
    userRunAction = 0;
    if(verboseLevel>1) G4cout << "UserRunAction deleted." << G4endl;
  }
  if(userPrimaryGeneratorAction)
  {
    delete userPrimaryGeneratorAction;
    userPrimaryGeneratorAction = 0;
    if(verboseLevel>1) G4cout << "UserPrimaryGenerator deleted." << G4endl;
  }

  if(verboseLevel>1) G4cout << "RunManager is deleting RunManagerKernel." << G4endl;

  delete kernel;

  fRunManager = 0;
}

void G4RunManager::DeleteUserInitializations()
{
    if( userDetector )
    {
        delete userDetector;
        userDetector = 0;
        if(verboseLevel>1) G4cout << "UserDetectorConstruction deleted." << G4endl;
    }
    if(physicsList)
    {
        delete physicsList;
        physicsList = 0;
        if(verboseLevel>1) G4cout << "UserPhysicsList deleted." << G4endl;
    }
    if(userActionInitialization)
    {
        delete userActionInitialization;
        userActionInitialization = 0;
        if(verboseLevel>1) G4cout <<"UserActionInitialization deleted." << G4endl;
    }
    if(userWorkerInitialization)
    {
        delete userWorkerInitialization;
        userWorkerInitialization = 0;
        if(verboseLevel>1) G4cout <<"UserWorkerInitialization deleted." << G4endl;
    }
    if(userWorkerThreadInitialization)
    {
        delete userWorkerThreadInitialization;
        userWorkerThreadInitialization = 0;
        if(verboseLevel>1) G4cout <<"UserWorkerThreadInitialization deleted." << G4endl;
    }

}

void G4RunManager::BeamOn(G4int n_event,const char* macroFile,G4int n_select)
{
  if(n_event<=0) { fakeRun = true; }
  else { fakeRun = false; }
  G4bool cond = ConfirmBeamOnCondition();
  if(cond)
  {
    numberOfEventToBeProcessed = n_event;
    numberOfEventProcessed = 0;
    ConstructScoringWorlds();
    RunInitialization();
    DoEventLoop(n_event,macroFile,n_select);
    RunTermination();
  }
  fakeRun = false;
}

G4bool G4RunManager::ConfirmBeamOnCondition()
{
  G4StateManager* stateManager = G4StateManager::GetStateManager();

  G4ApplicationState currentState = stateManager->GetCurrentState();
  if(currentState!=G4State_PreInit && currentState!=G4State_Idle)
  {
    G4cerr << "Illegal application state - BeamOn() ignored." << G4endl;
    return false;
  }

  if(!initializedAtLeastOnce)
  {
    G4cerr << " Geant4 kernel should be initialized" << G4endl;
    G4cerr << "before the first BeamOn(). - BeamOn ignored." << G4endl;
    return false;
  }

  if(!geometryInitialized || !physicsInitialized)
  {
    if(verboseLevel>0)
    {
      G4cout << "Start re-initialization because " << G4endl;
      if(!geometryInitialized) G4cout << "  Geometry" << G4endl;
      if(!physicsInitialized)  G4cout << "  Physics processes" << G4endl;
      G4cout << "has been modified since last Run." << G4endl;
    }
    Initialize();
  }
  return true;
}

void G4RunManager::RunInitialization()
{
  if(!(kernel->RunInitialization(fakeRun))) return;

  runAborted = false;
  numberOfEventProcessed = 0;

  CleanUpPreviousEvents();
  if(currentRun) delete currentRun;
  currentRun = 0;

  if(fakeRun) return;

  if(fGeometryHasBeenDestroyed) G4ParallelWorldProcessStore::GetInstance()->UpdateWorlds();

  if(userRunAction) currentRun = userRunAction->GenerateRun();
  if(!currentRun) currentRun = new G4Run();

  currentRun->SetRunID(runIDCounter);
  currentRun->SetNumberOfEventToBeProcessed(numberOfEventToBeProcessed);

  currentRun->SetDCtable(DCtable);
  G4SDManager* fSDM = G4SDManager::GetSDMpointerIfExist();
  if(fSDM)
  { currentRun->SetHCtable(fSDM->GetHCtable()); }
  
  std::ostringstream oss;
    G4Random::saveFullState(oss);
  randomNumberStatusForThisRun = oss.str();
  currentRun->SetRandomNumberStatus(randomNumberStatusForThisRun);
  
  for(G4int i_prev=0;i_prev<n_perviousEventsToBeStored;i_prev++)
  { previousEvents->push_back((G4Event*)0); }

  if(printModulo>=0 || verboseLevel>0)
  { G4cout << "### Run " << currentRun->GetRunID() << " starts." << G4endl; }
  if(userRunAction) userRunAction->BeginOfRunAction(currentRun);

  if(storeRandomNumberStatus) {
      G4String fileN = "currentRun";
      if ( rngStatusEventsFlag ) {
          std::ostringstream os;
          os << "run" << currentRun->GetRunID();
          fileN = os.str();
      }
      StoreRNGStatus(fileN);
  }
}

void G4RunManager::DoEventLoop(G4int n_event,const char* macroFile,G4int n_select)
{
  InitializeEventLoop(n_event,macroFile,n_select);

// Event loop
  for(G4int i_event=0; i_event<n_event; i_event++ )
  {
    ProcessOneEvent(i_event);
    TerminateOneEvent();
    if(runAborted) break;
  }

  // For G4MTRunManager, TerminateEventLoop() is invoked after all threads are finished.
  if(runManagerType==sequentialRM) TerminateEventLoop();
}

void G4RunManager::InitializeEventLoop(G4int n_event,const char* macroFile,G4int n_select)
{
  if(verboseLevel>0) 
  { timer->Start(); }

  n_select_msg = n_select;
  if(macroFile!=0)
  { 
    if(n_select_msg<0) n_select_msg = n_event;
    msgText = "/control/execute ";
    msgText += macroFile;
    selectMacro = macroFile;
  }
  else
  {
      n_select_msg = -1;
      selectMacro = "";
  }
}

void G4RunManager::ProcessOneEvent(G4int i_event)
{
  currentEvent = GenerateEvent(i_event);
  eventManager->ProcessOneEvent(currentEvent);
  AnalyzeEvent(currentEvent);
  UpdateScoring();
  if(i_event<n_select_msg) G4UImanager::GetUIpointer()->ApplyCommand(msgText);
}

void G4RunManager::TerminateOneEvent()
{
  StackPreviousEvent(currentEvent);
  currentEvent = 0;
  numberOfEventProcessed++;
}

void G4RunManager::TerminateEventLoop()
{
  if(verboseLevel>0 && !fakeRun)
  {
    timer->Stop();
    G4cout << " Run terminated." << G4endl;
    G4cout << "Run Summary" << G4endl;
    if(runAborted)
    { G4cout << "  Run Aborted after " << numberOfEventProcessed << " events processed." << G4endl; }
    else
    { G4cout << "  Number of events processed : " << numberOfEventProcessed << G4endl; }
    G4cout << "  "  << *timer << G4endl;
  }
  ////////////////  G4ProductionCutsTable::GetProductionCutsTable()->PhysicsTableUpdated();
  fGeometryHasBeenDestroyed = false;
}

G4Event* G4RunManager::GenerateEvent(G4int i_event)
{
  if(!userPrimaryGeneratorAction)
  {
    G4Exception("G4RunManager::GenerateEvent()", "Run0032", FatalException,
                "G4VUserPrimaryGeneratorAction is not defined!");
    return 0;
  }

  G4Event* anEvent = new G4Event(i_event);

  if(storeRandomNumberStatusToG4Event==1 || storeRandomNumberStatusToG4Event==3)
  {
    std::ostringstream oss;
    G4Random::saveFullState(oss);
    randomNumberStatusForThisEvent = oss.str();
    anEvent->SetRandomNumberStatus(randomNumberStatusForThisEvent);
  }

  if(storeRandomNumberStatus) {
      G4String fileN = "currentEvent";
      if ( rngStatusEventsFlag ) {
          std::ostringstream os;
          os << "run" << currentRun->GetRunID() << "evt" << anEvent->GetEventID();
          fileN = os.str();
      }
      StoreRNGStatus(fileN);
  }  
    
  if(printModulo > 0 && anEvent->GetEventID()%printModulo == 0 )
  { G4cout << "--> Event " << anEvent->GetEventID() << " starts." << G4endl; }
  userPrimaryGeneratorAction->GeneratePrimaries(anEvent);
  return anEvent;
}

void G4RunManager::StoreRNGStatus(const G4String& fnpref)
{
    G4String fileN = randomNumberStatusDir + fnpref+".rndm";
    G4Random::saveEngineStatus(fileN);
}

void G4RunManager::AnalyzeEvent(G4Event* anEvent)
{
  G4VPersistencyManager* fPersM = G4VPersistencyManager::GetPersistencyManager();
  if(fPersM) fPersM->Store(anEvent);
  currentRun->RecordEvent(anEvent);
}

void G4RunManager::RunTermination()
{
  if(!fakeRun)
  {
    CleanUpUnnecessaryEvents(0);
    if(userRunAction) userRunAction->EndOfRunAction(currentRun);
    G4VPersistencyManager* fPersM = G4VPersistencyManager::GetPersistencyManager();
    if(fPersM) fPersM->Store(currentRun);
    runIDCounter++;
  }

  kernel->RunTermination();
}

void G4RunManager::CleanUpPreviousEvents()
{
  // Delete all events carried over from previous run.
  // This method is invoked at the beginning of the next run
  // or from the destructor of G4RunManager at the very end of
  // the program.
  // N.B. If ToBeKept() is true, the pointer of this event is
  // kept in G4Run of the previous run, and deleted along with
  // the deletion of G4Run.

  std::list<G4Event*>::iterator evItr = previousEvents->begin();
  while(evItr!=previousEvents->end())
  {
    G4Event* evt = *evItr;
    if(evt && !(evt->ToBeKept())) delete evt;
    evItr = previousEvents->erase(evItr); 
  }
}

void G4RunManager::CleanUpUnnecessaryEvents(G4int keepNEvents)
{
  // Delete events that are no longer necessary for post
  // processing such as visualization.
  // N.B. If ToBeKept() is true, the pointer of this event is
  // kept in G4Run of the previous run, and deleted along with
  // the deletion of G4Run.

  std::list<G4Event*>::iterator evItr = previousEvents->begin();
  while(evItr!=previousEvents->end())
  {
    if(G4int(previousEvents->size()) <= keepNEvents) return;

    G4Event* evt = *evItr;
    if(evt)
    {
      if(evt->GetNumberOfGrips()==0) 
      {
        if(!(evt->ToBeKept())) delete evt; 
        evItr = previousEvents->erase(evItr);
      }
      else
      { evItr++; }
    }
    else
    { evItr = previousEvents->erase(evItr); }
  }
}

void G4RunManager::StackPreviousEvent(G4Event* anEvent)
{
  if(anEvent->ToBeKept()) currentRun->StoreEvent(anEvent);
  
  if(n_perviousEventsToBeStored==0)
  {
    if(anEvent->GetNumberOfGrips()==0) 
    { if(!(anEvent->ToBeKept())) delete anEvent; }
    else
    { previousEvents->push_back(anEvent); }
  }

  CleanUpUnnecessaryEvents(n_perviousEventsToBeStored);
}

void G4RunManager::Initialize()
{
  G4StateManager* stateManager = G4StateManager::GetStateManager();
  G4ApplicationState currentState = stateManager->GetCurrentState();
  if(currentState!=G4State_PreInit && currentState!=G4State_Idle)
  {
    G4cerr << "Illegal application state - "
           << "G4RunManager::Initialize() ignored." << G4endl;
    return;
  }

  stateManager->SetNewState(G4State_Init);
  if(!geometryInitialized) InitializeGeometry();
  if(!physicsInitialized) InitializePhysics();
  initializedAtLeastOnce = true;
  if(stateManager->GetCurrentState()!=G4State_Idle)
  { stateManager->SetNewState(G4State_Idle); }
}

void G4RunManager::InitializeGeometry()
{
  if(!userDetector)
  {
    G4Exception("G4RunManager::InitializeGeometry", "Run0033",
                FatalException, "G4VUserDetectorConstruction is not defined!");
    return;
  }

  if(verboseLevel>1) G4cout << "userDetector->Construct() start." << G4endl;

  G4StateManager* stateManager = G4StateManager::GetStateManager();
  G4ApplicationState currentState = stateManager->GetCurrentState();
  if(currentState==G4State_PreInit || currentState==G4State_Idle)
  { stateManager->SetNewState(G4State_Init); }
  kernel->DefineWorldVolume(userDetector->Construct(),false);
  userDetector->ConstructSDandField();
  nParallelWorlds = userDetector->ConstructParallelGeometries();
  userDetector->ConstructParallelSD();
  kernel->SetNumberOfParallelWorld(nParallelWorlds);
  geometryInitialized = true;
  stateManager->SetNewState(currentState); 
}

void G4RunManager::InitializePhysics()
{
  G4StateManager* stateManager = G4StateManager::GetStateManager();
  G4ApplicationState currentState = stateManager->GetCurrentState();
  if(currentState==G4State_PreInit || currentState==G4State_Idle)
  { stateManager->SetNewState(G4State_Init); }
  if(physicsList)
  {
    kernel->InitializePhysics();
  }
  else
  {
    G4Exception("G4RunManager::InitializePhysics()", "Run0034",
                FatalException, "G4VUserPhysicsList is not defined!");
  }
  physicsInitialized = true;
  stateManager->SetNewState(currentState); 

}

void G4RunManager::AbortRun(G4bool softAbort)
{
  // This method is valid only for GeomClosed or EventProc state
  G4ApplicationState currentState = 
    G4StateManager::GetStateManager()->GetCurrentState();
  if(currentState==G4State_GeomClosed || currentState==G4State_EventProc)
  {
    runAborted = true;
    if(currentState==G4State_EventProc && !softAbort)
    {
      currentEvent->SetEventAborted();
      eventManager->AbortCurrentEvent();
    }
  }
  else
  {
    G4cerr << "Run is not in progress. AbortRun() ignored." << G4endl;
  }
}

void G4RunManager::AbortEvent()
{
  // This method is valid only for EventProc state
  G4ApplicationState currentState = 
    G4StateManager::GetStateManager()->GetCurrentState();
  if(currentState==G4State_EventProc)
  {
    currentEvent->SetEventAborted();
    eventManager->AbortCurrentEvent();
  }
  else
  {
    G4cerr << "Event is not in progress. AbortEevnt() ignored." << G4endl;
  }
}

void G4RunManager::DefineWorldVolume(G4VPhysicalVolume* worldVol,
                                     G4bool topologyIsChanged)
{
  kernel->DefineWorldVolume(worldVol,topologyIsChanged);
}

void G4RunManager::rndmSaveThisRun()
{
  G4int runNumber = 0;
  if(currentRun) runNumber = currentRun->GetRunID();
  if(!storeRandomNumberStatus) {
     G4cerr << "Warning from G4RunManager::rndmSaveThisRun():"
          << " Random number status was not stored prior to this run." 
	  << G4endl << "Command ignored." << G4endl;
     return;
  }
  
  G4String fileIn  = randomNumberStatusDir + "currentRun.rndm";
 
  std::ostringstream os;
  os << "run" << runNumber << ".rndm" << '\0';
  G4String fileOut = randomNumberStatusDir + os.str();  

  G4String copCmd = "/control/shell cp "+fileIn+" "+fileOut;
  G4UImanager::GetUIpointer()->ApplyCommand(copCmd);
  if(verboseLevel>0) G4cout << "currentRun.rndm is copied to file: " << fileOut << G4endl;    
}

void G4RunManager::rndmSaveThisEvent()
{
  if(!storeRandomNumberStatus || currentEvent == 0) {
     G4cerr << "Warning from G4RunManager::rndmSaveThisEvent():"
          << " there is no currentEvent or its RandomEngineStatus is not available."
	  << G4endl << "Command ignored." << G4endl;
     return;
  }
  
  G4String fileIn  = randomNumberStatusDir + "currentEvent.rndm";

  std::ostringstream os;
  os << "run" << currentRun->GetRunID() << "evt" << currentEvent->GetEventID()
     << ".rndm" << '\0';
  G4String fileOut = randomNumberStatusDir + os.str();       

  G4String copCmd = "/control/shell cp "+fileIn+" "+fileOut;
  G4UImanager::GetUIpointer()->ApplyCommand(copCmd);
  if(verboseLevel>0) G4cout << "currentEvent.rndm is copied to file: " << fileOut << G4endl;  
}
  
void G4RunManager::RestoreRandomNumberStatus(const G4String& fileN)
{
  G4String fileNameWithDirectory;
  if(fileN.index("/")==std::string::npos)
  { fileNameWithDirectory = randomNumberStatusDir+fileN; }
  else
  { fileNameWithDirectory = fileN; }
  
  G4Random::restoreEngineStatus(fileNameWithDirectory);
  if(verboseLevel>0) G4cout << "RandomNumberEngineStatus restored from file: "
         << fileNameWithDirectory << G4endl;
  G4Random::showEngineStatus();
}

void G4RunManager::DumpRegion(const G4String& rname) const
{
  kernel->DumpRegion(rname);
}

void G4RunManager::DumpRegion(G4Region* region) const
{
  kernel->DumpRegion(region);
}

#include "G4ScoringManager.hh"
#include "G4TransportationManager.hh"
#include "G4VScoringMesh.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"

#include "G4ScoringBox.hh"
#include "G4ScoringCylinder.hh"

void G4RunManager::ConstructScoringWorlds()
{
  G4ScoringManager* ScM = G4ScoringManager::GetScoringManagerIfExist();
  if(!ScM) return;

  G4int nPar = ScM->GetNumberOfMesh();
  if(nPar<1) return;

  G4ParticleTable::G4PTblDicIterator* theParticleIterator
   = G4ParticleTable::GetParticleTable()->GetIterator();
  for(G4int iw=0;iw<nPar;iw++)
  {
    G4VScoringMesh* mesh = ScM->GetMesh(iw);
    if(fGeometryHasBeenDestroyed) mesh->GeometryHasBeenDestroyed();

    G4VPhysicalVolume* pWorld
       = G4TransportationManager::GetTransportationManager()
         ->IsWorldExisting(ScM->GetWorldName(iw));
    if(!pWorld)
    {
      pWorld = G4TransportationManager::GetTransportationManager()
           ->GetParallelWorld(ScM->GetWorldName(iw));
      pWorld->SetName(ScM->GetWorldName(iw));

      G4ParallelWorldProcess* theParallelWorldProcess
        = mesh->GetParallelWorldProcess();
      if(theParallelWorldProcess)
      { theParallelWorldProcess->SetParallelWorld(ScM->GetWorldName(iw)); }
      else
      {
        theParallelWorldProcess = new G4ParallelWorldProcess(ScM->GetWorldName(iw));
        mesh->SetParallelWorldProcess(theParallelWorldProcess);
        theParallelWorldProcess->SetParallelWorld(ScM->GetWorldName(iw));

        theParticleIterator->reset();
        while( (*theParticleIterator)() ){
          G4ParticleDefinition* particle = theParticleIterator->value();
          G4ProcessManager* pmanager = particle->GetProcessManager();
          if(pmanager)
          {
            pmanager->AddProcess(theParallelWorldProcess);
            if(theParallelWorldProcess->IsAtRestRequired(particle))
            { pmanager->SetProcessOrdering(theParallelWorldProcess, idxAtRest, 9900); }
            pmanager->SetProcessOrderingToSecond(theParallelWorldProcess, idxAlongStep);
            pmanager->SetProcessOrdering(theParallelWorldProcess, idxPostStep, 9900);
          }
        }
      }
    }

    mesh->Construct(pWorld);
  }

  GeometryHasBeenModified();
}

void G4RunManager::UpdateScoring()
{
  G4ScoringManager* ScM = G4ScoringManager::GetScoringManagerIfExist();
  if(!ScM) return;
  G4int nPar = ScM->GetNumberOfMesh();
  if(nPar<1) return;

  G4HCofThisEvent* HCE = currentEvent->GetHCofThisEvent();
  if(!HCE) return;
  G4int nColl = HCE->GetCapacity();
  for(G4int i=0;i<nColl;i++)
  {
    G4VHitsCollection* HC = HCE->GetHC(i);
    if(HC) ScM->Accumulate(HC);
  }
}

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4SmartVoxelHeader.hh"
#include "G4SmartVoxelStat.hh"

void G4RunManager::ReOptimizeMotherOf(G4VPhysicalVolume* pPhys)
{
  G4LogicalVolume* pMotherL = pPhys->GetMotherLogical();
  if(pMotherL) ReOptimize(pMotherL);
}

void G4RunManager::ReOptimize(G4LogicalVolume* pLog)
{
  G4Timer localtimer;
  if(verboseLevel>1)
  { localtimer.Start(); }
  G4SmartVoxelHeader* header = pLog->GetVoxelHeader();
  delete header;
  header = new G4SmartVoxelHeader(pLog);
  pLog->SetVoxelHeader(header);
  if(verboseLevel>1)
  {
    localtimer.Stop();
    G4SmartVoxelStat stat(pLog,header,localtimer.GetSystemElapsed(),
                          localtimer.GetUserElapsed());
    G4cout << G4endl << "Voxelisation of logical volume <"
           << pLog->GetName() << ">" << G4endl;
    G4cout << " heads : " << stat.GetNumberHeads() << " - nodes : "
           << stat.GetNumberNodes() << " - pointers : " 
           << stat.GetNumberPointers() << G4endl;
    G4cout << " Memory used : " << (stat.GetMemoryUse()+512)/1024
           << "k - total time : " << stat.GetTotalTime()
           << " - system time : " << stat.GetSysTime() << G4endl;
  }
}

void G4RunManager::SetUserInitialization(G4VUserDetectorConstruction* userInit)
{ userDetector = userInit; }

void G4RunManager::SetUserInitialization(G4VUserPhysicsList* userInit)
{
  physicsList = userInit;
  kernel->SetPhysics(userInit);
}

void G4RunManager::SetUserInitialization(G4UserWorkerInitialization* /*userInit*/)
{
  G4Exception("G4RunManager::SetUserInitialization()", "Run3001", FatalException,
    "Base-class G4RunManager cannot take G4UserWorkerInitialization. Use G4MTRunManager.");
}

void G4RunManager::SetUserInitialization(G4UserWorkerThreadInitialization* /*userInit*/)
{
  G4Exception("G4RunManager::SetUserThreadInitialization()", "Run3001", FatalException,
    "Base-class G4RunManager cannot take G4UserWorkerThreadInitialization. Use G4MTRunManager.");
}

void G4RunManager::SetUserInitialization(G4VUserActionInitialization* userInit)
{
  userActionInitialization = userInit; 
  userActionInitialization->Build();
}

void G4RunManager::SetUserAction(G4UserRunAction* userAction)
{
  userRunAction = userAction;
}

void G4RunManager::SetUserAction(G4VUserPrimaryGeneratorAction* userAction)
{
  userPrimaryGeneratorAction = userAction;
}

void G4RunManager::SetUserAction(G4UserEventAction* userAction)
{
  eventManager->SetUserAction(userAction);
  userEventAction = userAction;
}

void G4RunManager::SetUserAction(G4UserStackingAction* userAction)
{
  eventManager->SetUserAction(userAction);
  userStackingAction = userAction;
}

void G4RunManager::SetUserAction(G4UserTrackingAction* userAction)
{
  eventManager->SetUserAction(userAction);
  userTrackingAction = userAction;
}

void G4RunManager::SetUserAction(G4UserSteppingAction* userAction)
{
  eventManager->SetUserAction(userAction);
  userSteppingAction = userAction;
}

void G4RunManager::GeometryHasBeenModified(G4bool prop)
{
  if(prop)
  { G4UImanager::GetUIpointer()->ApplyCommand("/run/geometryModified"); }
  else
  { kernel->GeometryHasBeenModified(); }
}

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RegionStore.hh"

void G4RunManager::ReinitializeGeometry(G4bool destroyFirst, G4bool prop)
{
  if(destroyFirst && G4Threading::IsMasterThread())
  {
    if(verboseLevel>0)
    { 
      G4cout<<"#### G4PhysicalVolumeStore, G4LogicalVolumeStore and G4SolidStore\n"
            <<"#### are wiped out."<<G4endl; 
    }
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();

    // remove all logical volume pointers from regions
    // exception: world logical volume pointer must be kept
    G4RegionStore* regionStore = G4RegionStore::GetInstance();
    std::vector<G4Region*>::iterator rItr;
    for(rItr = regionStore->begin();rItr != regionStore->end(); rItr++)
    {
      if((*rItr)->GetName()=="DefaultRegionForTheWorld") continue;
      //if((*rItr)->GetName()=="DefaultRegionForParallelWorld") continue;
      std::vector<G4LogicalVolume*>::iterator lvItr
        = (*rItr)->GetRootLogicalVolumeIterator();
      for(size_t iRLV = 0;iRLV < (*rItr)->GetNumberOfRootVolumes(); iRLV++)
      {
        (*rItr)->RemoveRootLogicalVolume(*lvItr,false);
        lvItr++;
      }
      if(verboseLevel>0)
      { G4cout<<"#### Region <"<<(*rItr)->GetName()<<"> is cleared."<<G4endl; }
    }

    // clear transportation manager
    fGeometryHasBeenDestroyed = true;
    G4TransportationManager::GetTransportationManager()->ClearParallelWorlds();
  }
  if(prop)
  { G4UImanager::GetUIpointer()->ApplyCommand("/run/reinitializeGeometry"); }
  else
  {
    kernel->GeometryHasBeenModified();
    geometryInitialized = false;
    // Notify the VisManager as well
    if(G4Threading::IsMasterThread())
    {
      G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
      if(pVVisManager) pVVisManager->GeometryHasChanged();
    }
  }
}

