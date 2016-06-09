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
// $Id$
//
// 

// On Sun, to prevent conflict with ObjectSpace, G4Timer.hh has to be
// loaded *before* globals.hh...
#include "G4Timer.hh"

#include "G4RunManager.hh"
#include "G4RunManagerKernel.hh"

#include "G4StateManager.hh"
#include "G4ApplicationState.hh"
#include "Randomize.hh"
#include "G4Run.hh"
#include "G4RunMessenger.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VUserDetectorConstruction.hh"
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
#include "G4ios.hh"
#include <sstream>

using namespace CLHEP;

G4RunManager* G4RunManager::fRunManager = 0;

G4RunManager* G4RunManager::GetRunManager()
{ return fRunManager; }

G4RunManager::G4RunManager()
:userDetector(0),physicsList(0),
 userRunAction(0),userPrimaryGeneratorAction(0),userEventAction(0),
 userStackingAction(0),userTrackingAction(0),userSteppingAction(0),
 geometryInitialized(false),physicsInitialized(false),
 runAborted(false),initializedAtLeastOnce(false),
 geometryToBeOptimized(true),runIDCounter(0),verboseLevel(0),DCtable(0),
 currentRun(0),currentEvent(0),n_perviousEventsToBeStored(0),
 numberOfEventToBeProcessed(0),storeRandomNumberStatus(false),
 storeRandomNumberStatusToG4Event(0),
 currentWorld(0),nParallelWorlds(0),msgText(" "),n_select_msg(-1),
 numberOfEventProcessed(0)
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
  previousEvents = new std::vector<G4Event*>;
  G4ParticleTable::GetParticleTable()->CreateMessenger();
  G4ProcessTable::GetProcessTable()->CreateMessenger();
  randomNumberStatusDir = "./";
  std::ostringstream oss;
  HepRandom::saveFullState(oss);
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

  if(currentRun) delete currentRun;
  delete timer;
  delete runMessenger;
  G4ParticleTable::GetParticleTable()->DeleteMessenger();
  G4ProcessTable::GetProcessTable()->DeleteMessenger();
  delete previousEvents;
  if(userDetector)
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

  delete kernel;

  if(verboseLevel>1) G4cout << "RunManager is deleted." << G4endl;
  fRunManager = 0;
}

void G4RunManager::BeamOn(G4int n_event,const char* macroFile,G4int n_select)
{
  G4bool cond = ConfirmBeamOnCondition();
  if(cond)
  {
    numberOfEventToBeProcessed = n_event;
    ConstructScoringWorlds();
    RunInitialization();
    if(n_event>0) DoEventLoop(n_event,macroFile,n_select);
    RunTermination();
  }
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
  if(!(kernel->RunInitialization())) return;
  if(currentRun) delete currentRun;
  currentRun = 0;

  if(userRunAction) currentRun = userRunAction->GenerateRun();
  if(!currentRun) currentRun = new G4Run();

  currentRun->SetRunID(runIDCounter);
  currentRun->SetNumberOfEventToBeProcessed(numberOfEventToBeProcessed);

  currentRun->SetDCtable(DCtable);
  G4SDManager* fSDM = G4SDManager::GetSDMpointerIfExist();
  if(fSDM)
  { currentRun->SetHCtable(fSDM->GetHCtable()); }
  
  std::ostringstream oss;
  HepRandom::saveFullState(oss);
  randomNumberStatusForThisRun = oss.str();
  currentRun->SetRandomNumberStatus(randomNumberStatusForThisRun);
  
  previousEvents->clear();
  for(G4int i_prev=0;i_prev<n_perviousEventsToBeStored;i_prev++)
  { previousEvents->push_back((G4Event*)0); }

  if(userRunAction) userRunAction->BeginOfRunAction(currentRun);

  if(storeRandomNumberStatus) {
    G4String fileN = randomNumberStatusDir + "currentRun.rndm"; 
    HepRandom::saveEngineStatus(fileN);
  }

  runAborted = false;
  numberOfEventProcessed = 0;
  if(verboseLevel>0) G4cout << "Start Run processing." << G4endl;
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

  TerminateEventLoop();
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
  }
  else
  { n_select_msg = -1; }
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
  if(verboseLevel>0)
  {
    timer->Stop();
    G4cout << "Run terminated." << G4endl;
    G4cout << "Run Summary" << G4endl;
    if(runAborted)
    { G4cout << "  Run Aborted after " << numberOfEventProcessed << " events processed." << G4endl; }
    else
    { G4cout << "  Number of events processed : " << numberOfEventProcessed << G4endl; }
    G4cout << "  "  << *timer << G4endl;
  }
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
    HepRandom::saveFullState(oss);
    randomNumberStatusForThisEvent = oss.str();
    anEvent->SetRandomNumberStatus(randomNumberStatusForThisEvent);
  }

  if(storeRandomNumberStatus) {
    G4String fileN = randomNumberStatusDir + "currentEvent.rndm"; 
    HepRandom::saveEngineStatus(fileN);
  }  
    
  userPrimaryGeneratorAction->GeneratePrimaries(anEvent);
  return anEvent;
}

void G4RunManager::AnalyzeEvent(G4Event* anEvent)
{
  G4VPersistencyManager* fPersM = G4VPersistencyManager::GetPersistencyManager();
  if(fPersM) fPersM->Store(anEvent);
  currentRun->RecordEvent(anEvent);
}

void G4RunManager::RunTermination()
{
  for(size_t itr=0;itr<previousEvents->size();itr++)
  {
    G4Event* prevEv =  (*previousEvents)[itr];
    if((prevEv) && !(prevEv->ToBeKept())) delete prevEv;
  }
  previousEvents->clear();
  for(G4int i_prev=0;i_prev<n_perviousEventsToBeStored;i_prev++)
  { previousEvents->push_back((G4Event*)0); }

  if(userRunAction) userRunAction->EndOfRunAction(currentRun);

  G4VPersistencyManager* fPersM = G4VPersistencyManager::GetPersistencyManager();
  if(fPersM) fPersM->Store(currentRun);
  runIDCounter++;

  kernel->RunTermination();
}

void G4RunManager::StackPreviousEvent(G4Event* anEvent)
{
  if(anEvent->ToBeKept()) currentRun->StoreEvent(anEvent);
  G4Event* evt;
  if(n_perviousEventsToBeStored==0)
  { evt = anEvent; }
  else
  {
    previousEvents->insert(previousEvents->begin(),anEvent);
    evt = previousEvents->back();
    previousEvents->pop_back();
  }
  if(evt && !(evt->ToBeKept())) delete evt;
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

  if(!geometryInitialized) InitializeGeometry();
  if(!physicsInitialized) InitializePhysics();
  initializedAtLeastOnce = true;
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
  kernel->DefineWorldVolume(userDetector->Construct(),false);
  nParallelWorlds = userDetector->ConstructParallelGeometries();
  kernel->SetNumberOfParallelWorld(nParallelWorlds);
  geometryInitialized = true;
}

void G4RunManager::InitializePhysics()
{
  if(physicsList)
  {
    if(verboseLevel>1) G4cout << "physicsList->Construct() start." << G4endl;
    kernel->InitializePhysics();
  }
  else
  {
    G4Exception("G4RunManager::InitializePhysics()", "Run0034",
                FatalException, "G4VUserPhysicsList is not defined!");
  }
  physicsInitialized = true;
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
  
  HepRandom::restoreEngineStatus(fileNameWithDirectory);
  if(verboseLevel>0) G4cout << "RandomNumberEngineStatus restored from file: "
         << fileNameWithDirectory << G4endl;
  HepRandom::showEngineStatus();	 
}

void G4RunManager::DumpRegion(const G4String& rname) const
{
//  kernel->UpdateRegion();
  kernel->DumpRegion(rname);
}

void G4RunManager::DumpRegion(G4Region* region) const
{
//  kernel->UpdateRegion();
  kernel->DumpRegion(region);
}

#include "G4ScoringManager.hh"
#include "G4TransportationManager.hh"
#include "G4VScoringMesh.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ParallelWorldScoringProcess.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"

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
    G4VPhysicalVolume* pWorld
       = G4TransportationManager::GetTransportationManager()
         ->IsWorldExisting(ScM->GetWorldName(iw));
    if(!pWorld)
    {
      pWorld = G4TransportationManager::GetTransportationManager()
         ->GetParallelWorld(ScM->GetWorldName(iw));
      pWorld->SetName(ScM->GetWorldName(iw));

      G4ParallelWorldScoringProcess* theParallelWorldScoringProcess
        = new G4ParallelWorldScoringProcess(ScM->GetWorldName(iw));
      theParallelWorldScoringProcess->SetParallelWorld(ScM->GetWorldName(iw));

      theParticleIterator->reset();
      while( (*theParticleIterator)() ){
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        if(pmanager)
        {
          pmanager->AddProcess(theParallelWorldScoringProcess);
          if(theParallelWorldScoringProcess->IsAtRestRequired(particle))
          { pmanager->SetProcessOrdering(theParallelWorldScoringProcess, idxAtRest, 9999); }
          pmanager->SetProcessOrderingToSecond(theParallelWorldScoringProcess, idxAlongStep);
          pmanager->SetProcessOrdering(theParallelWorldScoringProcess, idxPostStep, 9999);
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

void G4RunManager::ReOptimizeMotherOf(G4VPhysicalVolume* pPhys)
{
  G4LogicalVolume* pMotherL = pPhys->GetMotherLogical();
  if(pMotherL) ReOptimize(pMotherL);
}

void G4RunManager::ReOptimize(G4LogicalVolume* pLog)
{
  G4SmartVoxelHeader* header = pLog->GetVoxelHeader();
  delete header;
  header = new G4SmartVoxelHeader(pLog);
  pLog->SetVoxelHeader(header);
}

