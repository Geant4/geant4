// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RunManager.cc,v 1.18 2000-12-14 10:38:14 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

// On Sun, to prevent conflict with ObjectSpace, G4Timer.hh has to be
// loaded *before* globals.hh...
#include "G4Timer.hh"

#include "G4RunManager.hh"

#include "Randomize.hh"
#include "G4Run.hh"
#include "G4RunMessenger.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPhysicsList.hh"
#include "G4UserRunAction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeometryManager.hh"
#include "G4SDManager.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ApplicationState.hh"
#include "G4StateManager.hh"
#include "G4VPersistencyManager.hh"
#include "G4UImanager.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessTable.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"

#include "G4ios.hh"
#include "g4std/strstream"


G4RunManager* G4RunManager::fRunManager = NULL;

G4RunManager* G4RunManager::GetRunManager()
{ return fRunManager; }

G4RunManager::G4RunManager()
:userDetector(NULL),physicsList(NULL),
 userRunAction(NULL),userPrimaryGeneratorAction(NULL),userEventAction(NULL),
 userStackingAction(NULL),userTrackingAction(NULL),userSteppingAction(NULL),
 currentRun(NULL),currentEvent(NULL),n_perviousEventsToBeStored(0),
 geometryInitialized(false),physicsInitialized(false),cutoffInitialized(false),
 geometryNeedsToBeClosed(true),initializedAtLeastOnce(false),
 runAborted(false),
 geometryToBeOptimized(true),verboseLevel(0),DCtable(NULL),runIDCounter(0),
 storeRandomNumberStatus(0)
{
  if(fRunManager)
  { G4Exception("G4RunManager constructed twice."); }
  //G4UnitDefinition::BuildUnitsTable();
  fRunManager = this;
  eventManager = new G4EventManager();
  timer = new G4Timer();
  runMessenger = new G4RunMessenger(this);
  previousEvents = new G4RWTPtrOrderedVector<G4Event>;
  G4ParticleTable::GetParticleTable()->CreateMessenger();
  G4ProcessTable::GetProcessTable()->CreateMessenger();
  randomNumberStatusDir = "./";
  G4cout 
  << "**********************************************" << G4endl
  << " Geant4 version $Name: not supported by cvs2svn $" << G4endl
  << "                                (15-Dec-2000)" << G4endl
  << "             Copyright : Geant4 Collaboration" << G4endl
  << "**********************************************" << G4endl;
}

G4RunManager::~G4RunManager()
{
  if(verboseLevel>0) G4cout << "G4 kernel has come to Quit state." << G4endl;
  G4StateManager* pStateManager = G4StateManager::GetStateManager();
  pStateManager->SetNewState(Quit);

  if(verboseLevel>1) G4cout << "Deletion of G4 kernel class start." << G4endl;
  delete timer;
  delete runMessenger;
  physicsList->RemoveProcessManager();
  G4ParticleTable::GetParticleTable()->DeleteMessenger();
  G4ProcessTable::GetProcessTable()->DeleteMessenger();
  delete previousEvents;
  if(userDetector)
  {
    delete userDetector;
    if(verboseLevel>1) G4cout << "UserDetectorConstruction deleted." << G4endl;
  }
  if(physicsList)
  {
    delete physicsList;
    if(verboseLevel>1) G4cout << "UserPhysicsList deleted." << G4endl;
  }
  if(userRunAction)
  {
    delete userRunAction;
    if(verboseLevel>1) G4cout << "UserRunAction deleted." << G4endl;
  }
  if(userPrimaryGeneratorAction)
  {
    delete userPrimaryGeneratorAction;
    if(verboseLevel>1) G4cout << "UserPrimaryGenerator deleted." << G4endl;
  }
  G4SDManager* fSDM = G4SDManager::GetSDMpointerIfExist();
  if(fSDM)
  {
    delete fSDM;
    if(verboseLevel>1) G4cout << "G4SDManager deleted." << G4endl;
  }
  delete eventManager;
  if(verboseLevel>1) G4cout << "EventManager deleted." << G4endl;
  G4UImanager* pUImanager = G4UImanager::GetUIpointer();
  {
    if(pUImanager) delete pUImanager;
    if(verboseLevel>1) G4cout << "UImanager deleted." << G4endl;
  }
  if(pStateManager) 
  {
    delete pStateManager;
    if(verboseLevel>1) G4cout << "StateManager deleted." << G4endl;
  }
  if(verboseLevel>1) G4cout << "RunManager is deleting." << G4endl;
}

void G4RunManager::BeamOn(G4int n_event,const char* macroFile,G4int n_select)
{
  G4bool cond = ConfirmBeamOnCondition();
  if(cond)
  {
    RunInitialization();
    DoEventLoop(n_event,macroFile,n_select);
    RunTermination();
  }
}

G4bool G4RunManager::ConfirmBeamOnCondition()
{
  G4StateManager* stateManager = G4StateManager::GetStateManager();

  G4ApplicationState currentState = stateManager->GetCurrentState();
  if(currentState!=PreInit && currentState!=Idle)
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

  if(!(geometryInitialized && physicsInitialized && cutoffInitialized)) 
  {
    if(verboseLevel>0)
    {
      G4cout << "Start re-initialization because " << G4endl;
      if(!geometryInitialized) G4cout << "  Geometry" << G4endl;
      if(!physicsInitialized)  G4cout << "  Physics processes" << G4endl;
      if(!cutoffInitialized)   G4cout << "  Cutoff" << G4endl;
      G4cout << "has been modified since last Run." << G4endl;
    }
    Initialize();
  }

  return true;
}

void G4RunManager::RunInitialization()
{
  currentRun = new G4Run();
  currentRun->SetRunID(runIDCounter);

  currentRun->SetDCtable(DCtable);
  G4SDManager* fSDM = G4SDManager::GetSDMpointerIfExist();
  if(fSDM)
  { currentRun->SetHCtable(fSDM->GetHCtable()); }
  
  if(userRunAction) userRunAction->BeginOfRunAction(currentRun);
  if(geometryNeedsToBeClosed)
  {
    if(verboseLevel>1) G4cout << "Start closing geometry." << G4endl;
    G4GeometryManager* geomManager = G4GeometryManager::GetInstance();
    geomManager->OpenGeometry();
    geomManager->CloseGeometry(geometryToBeOptimized);
    geometryNeedsToBeClosed = false;
  }
  G4StateManager* stateManager = G4StateManager::GetStateManager();
  stateManager->SetNewState(GeomClosed);

  previousEvents->clearAndDestroy();
  for(G4int i_prev=0;i_prev<n_perviousEventsToBeStored;i_prev++)
  { previousEvents->insert((G4Event*)NULL); }

  runAborted = false;

  if(storeRandomNumberStatus==1 || storeRandomNumberStatus==-1) StoreRandomNumberStatus();
  
  if(verboseLevel>0) G4cout << "Start Run processing." << G4endl;
}

void G4RunManager::DoEventLoop(G4int n_event,const char* macroFile,G4int n_select)
{
  G4StateManager* stateManager = G4StateManager::GetStateManager();

  if(verboseLevel>0) 
  { timer->Start(); }

  G4String msg;
  if(macroFile!=NULL)
  { 
    if(n_select<0) n_select = n_event;
    msg = "/control/execute ";
    msg += macroFile;
  }
  else
  { n_select = -1; }

  G4int i_event;
  for( i_event=0; i_event<n_event; i_event++ )
  {
    stateManager->SetNewState(EventProc);

    currentEvent = GenerateEvent(i_event);

    eventManager->ProcessOneEvent(currentEvent);

    AnalyzeEvent(currentEvent);

    if(i_event<n_select) G4UImanager::GetUIpointer()->ApplyCommand(msg);
    stateManager->SetNewState(GeomClosed);
    StackPreviousEvent(currentEvent);
    currentEvent = NULL;
    if(runAborted) break;
  }

  if(verboseLevel>0)
  {
    timer->Stop();
    G4cout << "Run terminated." << G4endl;
    G4cout << "Run Summary" << G4endl;
    if(runAborted)
    { G4cout << "  Run Aborted after " << i_event << " events processed." << G4endl; }
    else
    { G4cout << "  Number of events processed : " << n_event << G4endl; }
    G4cout << "  "  << *timer << G4endl;
  }
}

G4Event* G4RunManager::GenerateEvent(G4int i_event)
{
  if(!userPrimaryGeneratorAction)
  {
    G4Exception
    ("G4RunManager::BeamOn - G4VUserPrimaryGeneratorAction is not defined.");
  }

  G4Event* anEvent = new G4Event(i_event);

  if(storeRandomNumberStatus==2 || storeRandomNumberStatus==-2) StoreRandomNumberStatus(anEvent->GetEventID());

  userPrimaryGeneratorAction->GeneratePrimaries(anEvent);
  return anEvent;
}

void G4RunManager::AnalyzeEvent(G4Event* anEvent)
{
  G4VPersistencyManager* fPersM = G4VPersistencyManager::GetPersistencyManager();
  if(fPersM) fPersM->Store(currentEvent);
  currentRun->RecordEvent(currentEvent);
}

void G4RunManager::RunTermination()
{
  G4StateManager* stateManager = G4StateManager::GetStateManager();

  previousEvents->clearAndDestroy();

  if(userRunAction) userRunAction->EndOfRunAction(currentRun);

  G4VPersistencyManager* fPersM = G4VPersistencyManager::GetPersistencyManager();
  if(fPersM) fPersM->Store(currentRun);
  delete currentRun;
  currentRun = NULL;
  runIDCounter++;

  stateManager->SetNewState(Idle);
}

void G4RunManager::StackPreviousEvent(G4Event* anEvent)
{
  G4Event* evt;
  if(n_perviousEventsToBeStored==0)
  { evt = anEvent; }
  else
  {
    previousEvents->prepend(anEvent);
    evt = previousEvents->removeLast();
  }
  delete evt;
}

void G4RunManager::Initialize()
{
  G4StateManager* stateManager = G4StateManager::GetStateManager();
  G4ApplicationState currentState = stateManager->GetCurrentState();
  if(currentState!=PreInit && currentState!=Idle)
  {
    G4cerr << "Illegal application state - "
         << "G4RunManager::Initialize() ignored." << G4endl;
    return;
  }

  stateManager->SetNewState(Init);
  if(!geometryInitialized) InitializeGeometry();
  if(!physicsInitialized) InitializePhysics();
  if(!cutoffInitialized) InitializeCutOff();
  stateManager->SetNewState(Idle);
  if(!initializedAtLeastOnce) initializedAtLeastOnce = true;
}

void G4RunManager::InitializeGeometry()
{
  if(!userDetector)
  {
    G4Exception
    ("G4RunManager::InitializeGeometry - G4VUserDetectorConstruction is not defined.");
  }

  if(verboseLevel>1) G4cout << "userDetector->Construct() start." << G4endl;
  DefineWorldVolume(userDetector->Construct());
  geometryInitialized = true;
}

void G4RunManager::InitializePhysics()
{
  if(physicsList)
  {
    if(verboseLevel>1) G4cout << "physicsList->Construct() start." << G4endl;
    physicsList->Construct();
  }
  else
  {
    G4Exception("G4VUserPhysicsList is not defined");
  }
  physicsInitialized = true;
}

void G4RunManager::InitializeCutOff()
{
  if(physicsList)
  {
    if(verboseLevel>1) G4cout << "physicsList->setCut() start." << G4endl;
    physicsList->SetCuts();
  }
  cutoffInitialized = true;
}
  
void G4RunManager::AbortRun()
{
  // This method is valid only for GeomClosed or EventProc state
  G4ApplicationState currentState = 
    G4StateManager::GetStateManager()->GetCurrentState();
  if(currentState==GeomClosed || currentState==EventProc)
  {
    runAborted = true;
    if(currentState==EventProc) eventManager->AbortCurrentEvent();
  }
  else
  {
    G4cerr << "Run is not in progress. AbortRun() ignored." << G4endl;
  }
}

void G4RunManager::DefineWorldVolume(G4VPhysicalVolume* worldVol)
{
  // set the world volume to the Navigator
  G4TransportationManager::GetTransportationManager()
    ->GetNavigatorForTracking()
    ->SetWorldVolume(worldVol);

  // Let VisManager know it
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager) pVVisManager->GeometryHasChanged();

  geometryNeedsToBeClosed = true;
}

void G4RunManager::StoreRandomNumberStatus(G4int eventID)
{
  G4String fileN = randomNumberStatusDir+"RandEngine";
  if(storeRandomNumberStatus>0 && currentRun != NULL)
  {
    char st[20];
    G4std::ostrstream os(st,20);
    os << currentRun->GetRunID() << '\0';
    fileN += "R";
    fileN += st;
  }
  if(storeRandomNumberStatus==2 && eventID>=0)
  {
    char st[20];
    G4std::ostrstream os(st,20);
    os << eventID << '\0';
    fileN += "E";
    fileN += st;
  }
  fileN += ".stat";
  HepRandom::saveEngineStatus(fileN);
}
  
void G4RunManager::RestoreRandomNumberStatus(G4String fileN)
{
  G4String fileNameWithDirectory;
  if(fileN.index("/")==G4std::string::npos)
  { fileNameWithDirectory = randomNumberStatusDir+fileN; }
  else
  { fileNameWithDirectory = fileN; }
  HepRandom::restoreEngineStatus(fileNameWithDirectory);
}








