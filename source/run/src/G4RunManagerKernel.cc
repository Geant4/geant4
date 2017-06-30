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
// $Id: G4RunManagerKernel.cc 103291 2017-03-24 14:00:49Z gcosmo $
//
//

#include "G4RunManagerKernel.hh"

#include <vector>

#include "G4StateManager.hh"
#include "G4ApplicationState.hh"
#include "G4ExceptionHandler.hh"
#include "G4PrimaryTransformer.hh"
#include "G4GeometryManager.hh"
#include "G4FieldManagerStore.hh"
#include "G4NavigationHistoryPool.hh"
#include "G4PathFinder.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VUserPhysicsList.hh"

#include "G4ParticleTable.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"
#include "G4ProductionCutsTable.hh"
#include "G4SDManager.hh"
#include "G4ParallelWorldProcessStore.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"
#include "G4Version.hh"
#include "G4ios.hh"

#include "G4MTRunManager.hh"
#include "G4AllocatorList.hh"

#include "G4AutoLock.hh"

#ifdef G4FPE_DEBUG
  #include "G4FPEDetection.hh"
#endif

//The following lines are needed since G4VUserPhysicsList uses a #define theParticleIterator
#ifdef theParticleIterator
#undef theParticleIterator
#endif

G4ThreadLocal G4RunManagerKernel* G4RunManagerKernel::fRunManagerKernel = 0;

G4RunManagerKernel* G4RunManagerKernel::GetRunManagerKernel()
{ return fRunManagerKernel; }

G4RunManagerKernel::G4RunManagerKernel()
: physicsList(0),currentWorld(0),
 geometryInitialized(false),physicsInitialized(false),
 geometryToBeOptimized(true),
 physicsNeedsToBeReBuilt(true),verboseLevel(0),
 numberOfParallelWorld(0),geometryNeedsToBeClosed(true),
 numberOfStaticAllocators(0)
{
#ifdef G4FPE_DEBUG
  InvalidOperationDetection();
#endif
  G4AllocatorList* allocList = G4AllocatorList::GetAllocatorListIfExist();
  if(allocList) numberOfStaticAllocators = allocList->Size();
  defaultExceptionHandler = new G4ExceptionHandler();
  if(fRunManagerKernel)
  {
    G4Exception("G4RunManagerKernel::G4RunManagerKernel()","Run0001",
                FatalException,"More than one G4RunManagerKernel is constructed.");
  }
  fRunManagerKernel = this;

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  if(particleTable->entries()>0)
  {
    // No particle should be registered beforehand
    G4ExceptionDescription ED;
    ED<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<G4endl;
    ED<<" G4RunManagerKernel fatal exception"<<G4endl;
    ED<<"  -- Following particles have already been registered"<<G4endl;
    ED<<"     before G4RunManagerKernel is instantiated."<<G4endl;
    for(int i=0;i<particleTable->entries();i++)
    { ED<<"     "<<particleTable->GetParticle(i)->GetParticleName()<<G4endl; }
    ED<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<G4endl;
    G4Exception("G4RunManagerKernel::G4RunManagerKernel()","Run0002",
       FatalException,ED);
  }
  
  // construction of Geant4 kernel classes
  eventManager = new G4EventManager();

  defaultRegion = new G4Region("DefaultRegionForTheWorld"); // deleted by store
  defaultRegionForParallelWorld = new G4Region("DefaultRegionForParallelWorld"); // deleted by store
  defaultRegion->SetProductionCuts(
       G4ProductionCutsTable::GetProductionCutsTable()->GetDefaultProductionCuts());
  defaultRegionForParallelWorld->SetProductionCuts(
       G4ProductionCutsTable::GetProductionCutsTable()->GetDefaultProductionCuts());

  runManagerKernelType = sequentialRMK;
  // set the initial application state
  G4StateManager::GetStateManager()->SetNewState(G4State_PreInit);

  // version banner
  G4String vs = G4Version;
  vs = vs.substr(1,vs.size()-2);
  versionString = " Geant4 version ";
  versionString += vs;
  versionString += "   ";
  versionString += G4Date;
  G4cout << G4endl
    << "*************************************************************" << G4endl
    << versionString << G4endl
    << "                      Copyright : Geant4 Collaboration" << G4endl
    << "                      Reference : NIM A 506 (2003), 250-303" << G4endl
    << "                            WWW : http://cern.ch/geant4" << G4endl
    << "*************************************************************" << G4endl
    << G4endl;
}

G4RunManagerKernel::G4RunManagerKernel(RMKType rmkType)
: physicsList(0),currentWorld(0),
geometryInitialized(false),physicsInitialized(false),
geometryToBeOptimized(true),
physicsNeedsToBeReBuilt(true),verboseLevel(0),
numberOfParallelWorld(0),geometryNeedsToBeClosed(true),
 numberOfStaticAllocators(0)
{
//This version of the constructor should never be called in sequential mode!
#ifndef G4MULTITHREADED
    G4ExceptionDescription msg;
    msg<<"Geant4 code is compiled without multi-threading support (-DG4MULTITHREADED is set to off).";
    msg<<" This type of RunManagerKernel can only be used in mult-threaded applications.";
    G4Exception("G4RunManagerKernel::G4RunManagerKernel(G4bool)","Run0105",FatalException,msg);
#endif

#ifdef G4FPE_DEBUG
   if ( G4Threading::IsMasterThread() ) {
	   InvalidOperationDetection();
   }
#endif
    
    defaultExceptionHandler = new G4ExceptionHandler();
    if(fRunManagerKernel)
    {
        G4Exception("G4RunManagerKernel::G4RunManagerKernel()","Run0001",
                    FatalException,"More than one G4RunManagerKernel is constructed.");
    }
    fRunManagerKernel = this;
    // construction of Geant4 kernel classes
    eventManager = new G4EventManager();
    
    switch(rmkType)
    {
     case masterRMK:
        //Master thread behvior
        defaultRegion = new G4Region("DefaultRegionForTheWorld"); // deleted by store
        defaultRegionForParallelWorld = new G4Region("DefaultRegionForParallelWorld"); // deleted by store
        defaultRegion->SetProductionCuts(
               G4ProductionCutsTable::GetProductionCutsTable()->GetDefaultProductionCuts());
        defaultRegionForParallelWorld->SetProductionCuts(
               G4ProductionCutsTable::GetProductionCutsTable()->GetDefaultProductionCuts());
        break;
     case workerRMK:
        //Worker thread behavior
        defaultRegion = G4RegionStore::GetInstance()->GetRegion("DefaultRegionForTheWorld", true);
        defaultRegionForParallelWorld
            = G4RegionStore::GetInstance()->GetRegion("DefaultRegionForParallelWorld", true);
        break;
     default:   
        defaultRegion = 0;
        defaultRegionForParallelWorld = 0;
        G4ExceptionDescription msgx;
        msgx<<" This type of RunManagerKernel can only be used in mult-threaded applications.";
        G4Exception("G4RunManagerKernel::G4RunManagerKernel(G4bool)","Run0106",FatalException,msgx);
    }
    runManagerKernelType = rmkType;
    
    // set the initial application state
    G4StateManager::GetStateManager()->SetNewState(G4State_PreInit);
    
    // version banner
    G4String vs = G4Version;
    vs = vs.substr(1,vs.size()-2);
    switch(rmkType)
    {
     case masterRMK:
      versionString = " Geant4 version ";
      versionString += vs;
      versionString += "   ";
      versionString += G4Date;
      G4cout << G4endl
       << "*************************************************************" << G4endl
       << versionString << G4endl
       << "  << in Multi-threaded mode >> " << G4endl
       << "                      Copyright : Geant4 Collaboration" << G4endl
       << "                      Reference : NIM A 506 (2003), 250-303" << G4endl
       << "                            WWW : http://cern.ch/geant4" << G4endl
       << "*************************************************************" << G4endl
       << G4endl;
       break;
     default:
      if(verboseLevel) {
       versionString = " Local thread RunManagerKernel version ";
       versionString += vs;
       G4cout << G4endl
       << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << G4endl
       << versionString << G4endl
       << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << G4endl
       << G4endl;
      }
    }

#ifdef G4MULTITHREADED
    G4UnitDefinition::GetUnitsTable().Synchronize();
#endif
}

void G4RunManagerKernel::SetupDefaultRegion()
{
    if(runManagerKernelType==workerRMK) return;

    // Remove old world logical volume from the default region, if exist
    if(defaultRegion->GetNumberOfRootVolumes())
    {
        if(defaultRegion->GetNumberOfRootVolumes()>size_t(1))
        {
            G4Exception("G4RunManager::SetupDefaultRegion",
                        "Run0005",
                        FatalException,
                        "Default world region should have a unique logical volume.");
        }
        std::vector<G4LogicalVolume*>::iterator lvItr
        = defaultRegion->GetRootLogicalVolumeIterator();
        defaultRegion->RemoveRootLogicalVolume(*lvItr,false);
        if(verboseLevel>1) G4cout
            << "Obsolete world logical volume is removed from the default region." << G4endl;
    }

}

G4RunManagerKernel::~G4RunManagerKernel()
{
  G4StateManager* pStateManager = G4StateManager::GetStateManager();
  // set the application state to the quite state
  if(pStateManager->GetCurrentState()!=G4State_Quit)
  {
    if(verboseLevel>0) G4cout << "G4 kernel has come to Quit state." << G4endl;
    pStateManager->SetNewState(G4State_Quit);
  }

  // open geometry for deletion
  G4GeometryManager::GetInstance()->OpenGeometry();

  // deletion of Geant4 kernel classes
  G4ParallelWorldProcessStore* pwps = G4ParallelWorldProcessStore::GetInstanceIfExist();
  if(pwps) delete pwps;
  G4SDManager* fSDM = G4SDManager::GetSDMpointerIfExist();
  if(fSDM)
  {
    delete fSDM;
    if(verboseLevel>1) G4cout << "G4SDManager deleted." << G4endl;
  }
  delete eventManager;
  if(verboseLevel>1) G4cout << "EventManager deleted." << G4endl;

  G4UnitDefinition::ClearUnitsTable();
  if(verboseLevel>1) G4cout << "Units table cleared." << G4endl;

  // deletion of path-finder field-manager store, geometry and transportation manager
  G4PathFinder* pFinder = G4PathFinder::GetInstanceIfExist();
  if (pFinder) delete pFinder;
  G4FieldManagerStore* fmStore = G4FieldManagerStore::GetInstanceIfExist();
  if (fmStore) delete fmStore;
  G4GeometryManager* gManager = G4GeometryManager::GetInstanceIfExist();
  if (gManager) delete gManager;
  G4TransportationManager* tManager = G4TransportationManager::GetInstanceIfExist();
  if (tManager)
  {
    delete tManager;
    if(verboseLevel>1) G4cout << "TransportationManager deleted." << G4endl;
  }

  // deletion of navigation levels
  if(verboseLevel>1) G4NavigationHistoryPool::GetInstance()->Print();
  delete G4NavigationHistoryPool::GetInstance();

  // deletion of allocators
  G4AllocatorList* allocList = G4AllocatorList::GetAllocatorListIfExist();
  if(allocList)
  {
    allocList->Destroy(numberOfStaticAllocators,verboseLevel);
    delete allocList; 
    if(verboseLevel>1) G4cout << "G4Allocator objects are deleted." << G4endl;
  }

  G4UImanager* pUImanager = G4UImanager::GetUIpointer();
  if((runManagerKernelType==workerRMK) && (verboseLevel>0))
  {
    G4cout << "Thread-local UImanager is to be deleted." << G4endl 
           << "There should not be any thread-local G4cout/G4cerr hereafter."
           << G4endl;
    verboseLevel = 0;
  }
  if(pUImanager) delete pUImanager;
  if(verboseLevel>1) G4cout << "UImanager deleted." << G4endl;
  
  delete pStateManager; 
  if(verboseLevel>1) G4cout << "StateManager deleted." << G4endl;
  delete defaultExceptionHandler;
  if(verboseLevel>0) G4cout << "RunManagerKernel is deleted. Good bye :)" << G4endl;
  fRunManagerKernel = 0;
}

void G4RunManagerKernel::WorkerUpdateWorldVolume()
{
  G4MTRunManager* masterRM = G4MTRunManager::GetMasterRunManager();
  G4TransportationManager* transM = G4TransportationManager::GetTransportationManager();
  G4MTRunManager::masterWorlds_t masterWorlds= masterRM->GetMasterWorlds();
  G4MTRunManager::masterWorlds_t::iterator itrMW = masterWorlds.begin();
  for(;itrMW!=masterWorlds.end();itrMW++)
  {
    G4VPhysicalVolume* wv = (*itrMW).second;
    G4VPhysicalVolume* pWorld
       = G4TransportationManager::GetTransportationManager()
         ->IsWorldExisting(wv->GetName());
    if(!pWorld)
    { transM->RegisterWorld(wv); }
  }
}

void G4RunManagerKernel::WorkerDefineWorldVolume(G4VPhysicalVolume* worldVol,
                                     G4bool topologyIsChanged)
{
  G4StateManager*    stateManager = G4StateManager::GetStateManager();
  G4ApplicationState currentState = stateManager->GetCurrentState();
  if(currentState!=G4State_Init)
  {
    if(!(currentState==G4State_Idle||currentState==G4State_PreInit))
    {
      G4cout << "Current application state is "
        << stateManager->GetStateString(currentState) << G4endl;
      G4Exception("G4RunManagerKernel::DefineWorldVolume",
		"DefineWorldVolumeAtIncorrectState",
		FatalException,
		"Geant4 kernel is not Init state : Method ignored.");
      return;
    } else {
      //G4Exception("G4RunManagerKernel::DefineWorldVolume",
      //        "DefineWorldVolumeAtIncorrectState",
      //        JustWarning,
      //        "Geant4 kernel is not Init state : Assuming Init state.");
      //G4cout<<"Warning : Geant4 kernel is not Init state : Assuming Init state."
      //      <<G4endl;
      stateManager->SetNewState(G4State_Init); 
    }
  }

  currentWorld = worldVol;
  G4MTRunManager* masterRM = G4MTRunManager::GetMasterRunManager();
  G4TransportationManager* transM = G4TransportationManager::GetTransportationManager();
  G4MTRunManager::masterWorlds_t masterWorlds= masterRM->GetMasterWorlds();
  G4MTRunManager::masterWorlds_t::iterator itrMW = masterWorlds.begin();
  for(;itrMW!=masterWorlds.end();itrMW++)
  {
    if((*itrMW).first == 0)
    {
      if((*itrMW).second != currentWorld)
      {
        G4Exception("G4RunManagerKernel::WorkerDefineWorldVolume","RUN3091",
                    FatalException,"Mass world is inconsistent");
      }
      transM->SetWorldForTracking((*itrMW).second);
    }
    else
    {
      transM->RegisterWorld((*itrMW).second);
    }
  }

  if(topologyIsChanged) geometryNeedsToBeClosed = true;

  // Notify the VisManager as well                                                                                
  if(G4Threading::IsMasterThread())
  {
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if(pVVisManager) pVVisManager->GeometryHasChanged();
  }

  geometryInitialized = true;
  stateManager->SetNewState(currentState); 
  if(physicsInitialized && currentState!=G4State_Idle)
  { stateManager->SetNewState(G4State_Idle); }
}

void G4RunManagerKernel::DefineWorldVolume(G4VPhysicalVolume* worldVol,
                                     G4bool topologyIsChanged)
{
  G4StateManager*    stateManager = G4StateManager::GetStateManager();
  G4ApplicationState currentState = stateManager->GetCurrentState();
  if(currentState!=G4State_Init)
  {
    if(!(currentState==G4State_Idle||currentState==G4State_PreInit))
    {
      G4cout << "Current application state is "
        << stateManager->GetStateString(currentState) << G4endl;
      G4Exception("G4RunManagerKernel::DefineWorldVolume",
		"DefineWorldVolumeAtIncorrectState",
		FatalException,
		"Geant4 kernel is not Init state : Method ignored.");
      return;
    } else {
      //G4Exception("G4RunManagerKernel::DefineWorldVolume",
      //        "DefineWorldVolumeAtIncorrectState",
      //        JustWarning,
      //        "Geant4 kernel is not Init state : Assuming Init state.");
      //G4cout<<"Warning : Geant4 kernel is not Init state : Assuming Init state."
      //      <<G4endl;
      stateManager->SetNewState(G4State_Init); 
    }
  }

  // The world volume MUST NOT have a region defined by the user
  if(worldVol->GetLogicalVolume()->GetRegion())
  {
    if(worldVol->GetLogicalVolume()->GetRegion()!=defaultRegion)
    {
      G4ExceptionDescription ED;
      ED << "The world volume has a user-defined region <"
           << worldVol->GetLogicalVolume()->GetRegion()->GetName()
           << ">." << G4endl;
      ED << "World would have a default region assigned by RunManagerKernel."
         << G4endl;
      G4Exception("G4RunManager::DefineWorldVolume",
                "Run0004", FatalException, ED);
    }
  }

  SetupDefaultRegion();
    
  // Accept the world volume
  currentWorld = worldVol; 

  // Set the default region to the world

  G4LogicalVolume* worldLog = currentWorld->GetLogicalVolume();
  worldLog->SetRegion(defaultRegion);
  defaultRegion->AddRootLogicalVolume(worldLog);
  if(verboseLevel>1) G4cout << worldLog->GetName()
   << " is registered to the default region." << G4endl;

  // Set the world volume, notify the Navigator and reset its state
  G4TransportationManager::GetTransportationManager()
      ->SetWorldForTracking(currentWorld);
  if(topologyIsChanged) geometryNeedsToBeClosed = true;
  
  // Notify the VisManager as well
  if(G4Threading::IsMasterThread())
  {
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if(pVVisManager) pVVisManager->GeometryHasChanged();
  }

  geometryInitialized = true;
  stateManager->SetNewState(currentState); 
  if(physicsInitialized && currentState!=G4State_Idle)
  { stateManager->SetNewState(G4State_Idle); }
} 
  
void G4RunManagerKernel::SetPhysics(G4VUserPhysicsList* uPhys)
{
  physicsList = uPhys;

  if(runManagerKernelType==workerRMK) return;

  SetupPhysics();
  if(verboseLevel>2) G4ParticleTable::GetParticleTable()->DumpTable();
  if(verboseLevel>1)
  {
    G4cout << "List of instantiated particles ============================================" << G4endl;
    G4int nPtcl = G4ParticleTable::GetParticleTable()->entries();
    for(G4int i=0;i<nPtcl;i++)
    {
      G4ParticleDefinition* pd = G4ParticleTable::GetParticleTable()->GetParticle(i);
      G4cout << pd->GetParticleName() << " ";
      if(i%10==9) G4cout << G4endl;
    }
    G4cout << G4endl;
  }
}

#include "G4IonTable.hh"
#include "G4ParticleTableIterator.hh"
#include "G4IonConstructor.hh"
#include "G4Geantino.hh"

void G4RunManagerKernel::SetupPhysics()
{
    G4ParticleTable::GetParticleTable()->SetReadiness();

    physicsList->ConstructParticle();

    // For sanity reason
    G4Geantino::GeantinoDefinition();
    G4ParticleDefinition* gion = G4ParticleTable::GetParticleTable()->GetGenericIon();
    if(gion)
    { G4IonConstructor::ConstructParticle(); }
    G4ParticleTable::GetParticleTable()->GetIonTable()->InitializeLightIons();

    G4ParticleTable::G4PTblDicIterator* pItr = G4ParticleTable::GetParticleTable()->GetIterator();
    pItr->reset();
    while( (*pItr)() )
    {
      G4ParticleDefinition* particle = pItr->value();
      if(!(particle->IsGeneralIon())) particle->SetParticleDefinitionID();
    }

    if(gion)
    {
      G4int gionId = gion->GetParticleDefinitionID();
      pItr->reset(false);
      while( (*pItr)() )
      {
        G4ParticleDefinition* particle = pItr->value();
        if(particle->IsGeneralIon()) particle->SetParticleDefinitionID(gionId);
      }
    }
#ifdef G4MULTITHREADED
    G4UnitDefinition::GetUnitsTable().Synchronize();
#endif
}

namespace {
    G4Mutex initphysicsmutex = G4MUTEX_INITIALIZER;
}

void G4RunManagerKernel::InitializePhysics()
{
  G4StateManager*    stateManager = G4StateManager::GetStateManager();
  G4ApplicationState currentState = stateManager->GetCurrentState();
  if(currentState!=G4State_Init)
  {
    G4cout << "Current application state is "
        << stateManager->GetStateString(currentState) << G4endl;
    if(!(currentState==G4State_Idle||currentState==G4State_PreInit))
    {
      G4Exception("G4RunManagerKernel::InitializePhysics",
		"InitializePhysicsIncorrectState",
		FatalException,
		"Geant4 kernel is not Init state : Method ignored.");
      return;
    } else {
      //G4Exception("G4RunManagerKernel::DefineWorldVolume",
	//"DefineWorldVolumeAtIncorrectState",
	//JustWarning,
	//"Geant4 kernel is not Init state : Assuming Init state.");
      G4cout<<"Warning : Geant4 kernel is not Init state : Assuming Init state."
            <<G4endl;
      stateManager->SetNewState(G4State_Init); 
    }
  }

  if(!physicsList)
  {
    G4Exception("G4RunManagerKernel::InitializePhysics",
                "Run0012", FatalException,
                "G4VUserPhysicsList is not defined");
    return;
  }

  if(verboseLevel>1) G4cout << "physicsList->Construct() start." << G4endl;
  if(numberOfParallelWorld>0) physicsList->UseCoupledTransportation();
  physicsList->Construct();
    
  if(verboseLevel>1) G4cout << "physicsList->CheckParticleList() start." << G4endl;
  physicsList->CheckParticleList();
    //Cannot assume that SetCuts and CheckRegions are thread safe. We need to mutex
    //Report from valgrind --tool=drd
  G4AutoLock l(&initphysicsmutex);
  if ( G4Threading::IsMasterThread() ) {
	  if(verboseLevel>1) G4cout << "physicsList->setCut() start." << G4endl;
	  physicsList->SetCuts();

  }
  CheckRegions();
  l.unlock();

/*******************
//  static G4bool createIsomerOnlyOnce = false;
//  if(G4Threading::IsMultithreadedApplication() && G4Threading::IsMasterThread())
//  {
//    if(!createIsomerOnlyOnce)
//    {
//      createIsomerOnlyOnce = true;
//      G4ParticleDefinition* gion = G4ParticleTable::GetParticleTable()->GetGenericIon();
//      if(gion)
//      {
//        G4ParticleTable::GetParticleTable()->GetIonTable()->CreateAllIsomer();
//        PropagateGenericIonID();
//      }
//    }
//  }
*********************/

  physicsInitialized = true;
#ifdef G4MULTITHREADED
  G4UnitDefinition::GetUnitsTable().Synchronize();
#endif
  stateManager->SetNewState(currentState); 
  if(geometryInitialized && currentState!=G4State_Idle)
  { stateManager->SetNewState(G4State_Idle); }
}

G4bool G4RunManagerKernel::RunInitialization(G4bool fakeRun)
{
  G4StateManager*    stateManager = G4StateManager::GetStateManager();
  G4ApplicationState currentState = stateManager->GetCurrentState();

  if(!geometryInitialized) 
  { 
    G4Exception("G4RunManagerKernel::RunInitialization",
                "Run0021",
                JustWarning,
                "Geometry has not yet initialized : method ignored.");
    return false;
  }
  
  if(!physicsInitialized) 
  { 
    G4Exception("G4RunManagerKernel::RunInitialization",
                "Run0022",
                JustWarning,
                "Physics has not yet initialized : method ignored.");
    return false;
  }

  if( currentState != G4State_Idle )
  { 
    G4Exception("G4RunManagerKernel::RunInitialization",
                "Run0023",
                JustWarning,
                "Geant4 kernel not in Idle state : method ignored.");
    return false;
  }

  if(geometryNeedsToBeClosed) CheckRegularGeometry();

  PropagateGenericIonID();
  SetupShadowProcess();
  UpdateRegion();
  BuildPhysicsTables(fakeRun);

  if(geometryNeedsToBeClosed)
  {
    ResetNavigator();
    // CheckRegularGeometry();
    // Notify the VisManager as well
    if(G4Threading::IsMasterThread())
    {
      G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
      if(pVVisManager) pVVisManager->GeometryHasChanged();
    }
  }
 
  GetPrimaryTransformer()->CheckUnknown();

#ifdef G4MULTITHREADED
  G4UnitDefinition::GetUnitsTable().Synchronize();
#endif
  stateManager->SetNewState(G4State_GeomClosed);
  return true;
}

void G4RunManagerKernel::PropagateGenericIonID()
{
  G4ParticleDefinition* gion = G4ParticleTable::GetParticleTable()->GetGenericIon();
  if(gion)
  {
    //G4ParticleTable::GetParticleTable()->GetIonTable()->CreateAllIsomer();
    G4int gionId = gion->GetParticleDefinitionID();
    G4ParticleTable::G4PTblDicIterator* pItr = G4ParticleTable::GetParticleTable()->GetIterator();
    pItr->reset(false);
    while( (*pItr)() )
    {
      G4ParticleDefinition* particle = pItr->value();
      if(particle->IsGeneralIon()) particle->SetParticleDefinitionID(gionId);
    }
  }
}

void G4RunManagerKernel::RunTermination()
{
  if ( runManagerKernelType != workerRMK )
      G4ProductionCutsTable::GetProductionCutsTable()->PhysicsTableUpdated();
  G4StateManager::GetStateManager()->SetNewState(G4State_Idle); 
}

void G4RunManagerKernel::ResetNavigator()
{
  if(runManagerKernelType==workerRMK)
  { geometryNeedsToBeClosed = false; return; }

  // We have to tweak the navigator's state in case a geometry has been
  // modified between runs. By the following calls we ensure that navigator's
  // state is reset properly. It is required the geometry to be closed
  // and previous optimisations to be cleared.
  
  G4GeometryManager* geomManager = G4GeometryManager::GetInstance();
  if(verboseLevel>1) G4cout << "Start closing geometry." << G4endl;

  geomManager->OpenGeometry();
  geomManager->CloseGeometry(geometryToBeOptimized, verboseLevel>1);
 
  geometryNeedsToBeClosed = false;
}

void G4RunManagerKernel::UpdateRegion()
{
  G4StateManager*    stateManager = G4StateManager::GetStateManager();
  G4ApplicationState currentState = stateManager->GetCurrentState();
  if( currentState != G4State_Idle )
  { 
    G4Exception("G4RunManagerKernel::UpdateRegion",
                "Run0024",
                JustWarning,
                "Geant4 kernel not in Idle state : method ignored.");
    return;
  }

  if(runManagerKernelType==workerRMK) return;

  CheckRegions();

  G4RegionStore::GetInstance()->UpdateMaterialList(currentWorld);

  G4ProductionCutsTable::GetProductionCutsTable()->UpdateCoupleTable(currentWorld);
}

void G4RunManagerKernel::BuildPhysicsTables(G4bool fakeRun)
{
  if( G4ProductionCutsTable::GetProductionCutsTable()->IsModified()
  || physicsNeedsToBeReBuilt)
  {
#ifdef G4MULTITHREADED
    if(runManagerKernelType==masterRMK)
    {
      // make sure workers also rebuild physics tables
      G4UImanager* pUImanager = G4UImanager::GetUIpointer();
      pUImanager->ApplyCommand("/run/physicsModified");
    }
#endif
    physicsList->BuildPhysicsTable();
    ////G4ProductionCutsTable::GetProductionCutsTable()->PhysicsTableUpdated();
    physicsNeedsToBeReBuilt = false;
  }

  if(!fakeRun && verboseLevel>1) DumpRegion();
  if(!fakeRun && verboseLevel>0) physicsList->DumpCutValuesTable();
  if(!fakeRun) physicsList->DumpCutValuesTableIfRequested();
}

void G4RunManagerKernel::CheckRegions()
{
  G4TransportationManager* transM = G4TransportationManager::GetTransportationManager();
  size_t nWorlds = transM->GetNoWorlds();
  std::vector<G4VPhysicalVolume*>::iterator wItr;
  for(size_t i=0;i<G4RegionStore::GetInstance()->size();i++)
  { 
    G4Region* region = (*(G4RegionStore::GetInstance()))[i];

    //Let each region have a pointer to the world volume where it belongs to.
    //G4Region::SetWorld() checks if the region belongs to the given world and set it
    //only if it does. Thus, here we go through all the registered world volumes.
    region->SetWorld(0); // reset
    region->UsedInMassGeometry(false);
    region->UsedInParallelGeometry(false);
    wItr = transM->GetWorldsIterator();
    for(size_t iw=0;iw<nWorlds;iw++)
    {
      if(region->BelongsTo(*wItr))
      {
        if(*wItr==currentWorld)
        { region->UsedInMassGeometry(true); }
        else
        { region->UsedInParallelGeometry(true); }
      }
      region->SetWorld(*wItr);
      wItr++;
    }

    G4ProductionCuts* cuts = region->GetProductionCuts();
    if(!cuts)
    {
      if(region->IsInMassGeometry())
      {
        G4cout << "Warning : Region <" << region->GetName()
             << "> does not have specific production cuts," << G4endl
             << "even though it appears in the current tracking world." << G4endl;
        G4cout << "Default cuts are used for this region." << G4endl;
      }

      if(region->IsInMassGeometry()||region->IsInParallelGeometry())
      {
        region->SetProductionCuts(
          G4ProductionCutsTable::GetProductionCutsTable()
             ->GetDefaultProductionCuts());
      }
    }
  }

  //
  // If a parallel world has no region, set default region for parallel world
  //

  wItr = transM->GetWorldsIterator();
  for(size_t iw=0;iw<nWorlds;iw++)
  {
    //G4cout << "+++ " << (*wItr)->GetName() << G4endl;
    if(*wItr!=currentWorld)
    {
      G4LogicalVolume* pwLogical = (*wItr)->GetLogicalVolume();
      if(!(pwLogical->GetRegion()))
      {
        pwLogical->SetRegion(defaultRegionForParallelWorld);
        defaultRegionForParallelWorld->AddRootLogicalVolume(pwLogical);
        //G4cout << "+++++ defaultRegionForParallelWorld is set to "
        //       << (*wItr)->GetName() << " +++++" << G4endl;
      }
    }
    wItr++;
  }

}

void G4RunManagerKernel::DumpRegion(const G4String& rname) const
{
  G4Region* region = G4RegionStore::GetInstance()->GetRegion(rname);
  if(region) DumpRegion(region);
}

void G4RunManagerKernel::DumpRegion(G4Region* region) const
{
  if(!region)
  {
    for(size_t i=0;i<G4RegionStore::GetInstance()->size();i++)
    { DumpRegion((*(G4RegionStore::GetInstance()))[i]); }
  }
  else
  {
    if(G4Threading::IsWorkerThread()) return;
    G4cout << G4endl;
    G4cout << "Region <" << region->GetName() << "> -- ";
    if(region->GetWorldPhysical())
    {
      G4cout << " -- appears in <" 
           << region->GetWorldPhysical()->GetName() << "> world volume";
    }
    else
    { G4cout << " -- is not associated to any world."; }
    G4cout << G4endl;
    if(region->IsInMassGeometry())
    { G4cout << " This region is in the mass world." << G4endl; }
    if(region->IsInParallelGeometry())
    { G4cout << " This region is in the parallel world." << G4endl; }

    G4cout << " Root logical volume(s) : ";
    size_t nRootLV = region->GetNumberOfRootVolumes();
    std::vector<G4LogicalVolume*>::iterator lvItr = region->GetRootLogicalVolumeIterator();
    for(size_t j=0;j<nRootLV;j++)
    { G4cout << (*lvItr)->GetName() << " "; lvItr++; }
    G4cout << G4endl;

    G4cout << " Pointers : G4VUserRegionInformation[" << region->GetUserInformation() 
           << "], G4UserLimits[" << region->GetUserLimits() 
           << "], G4FastSimulationManager[" << region->GetFastSimulationManager()
           << "], G4UserSteppingAction[" << region->GetRegionalSteppingAction() << "]" << G4endl;
    
    G4cout << " Materials : ";
    std::vector<G4Material*>::const_iterator mItr = region->GetMaterialIterator();
    size_t nMaterial = region->GetNumberOfMaterials();
    for(size_t iMate=0;iMate<nMaterial;iMate++)
    {
      G4cout << (*mItr)->GetName() << " ";
      mItr++;
    }
    G4cout << G4endl;
    G4ProductionCuts* cuts = region->GetProductionCuts();
    if(!cuts && region->IsInMassGeometry())
    {
      G4cerr << "Warning : Region <" << region->GetName()
             << "> does not have specific production cuts." << G4endl;
      G4cerr << "Default cuts are used for this region." << G4endl;
      region->SetProductionCuts(
          G4ProductionCutsTable::GetProductionCutsTable()->GetDefaultProductionCuts());
    }
    else if(cuts)
    {
      G4cout << " Production cuts : "
             << "  gamma "
             << G4BestUnit(cuts->GetProductionCut("gamma"),"Length")
             << "     e- "
             << G4BestUnit(cuts->GetProductionCut("e-"),"Length")
             << "     e+ "
             << G4BestUnit(cuts->GetProductionCut("e+"),"Length")
             << " proton "
             << G4BestUnit(cuts->GetProductionCut("proton"),"Length")
             << G4endl;
    }
  }
}

#include "G4LogicalVolumeStore.hh"
void G4RunManagerKernel::CheckRegularGeometry()
{
  G4LogicalVolumeStore* store = G4LogicalVolumeStore::GetInstance();
  for(G4LogicalVolumeStore::iterator pos=store->begin(); pos!=store->end(); pos++)
  {
    if((*pos)&&((*pos)->GetNoDaughters()==1))
    {
      if((*pos)->GetDaughter(0)->IsRegularStructure())
      {
        SetScoreSplitter();
        return;
      }
    }
  }
}
        
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4VProcess.hh"
G4bool G4RunManagerKernel::ConfirmCoupledTransportation()
{
  G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();
  G4ParticleTable::G4PTblDicIterator* theParticleIterator = theParticleTable->GetIterator();
  theParticleIterator->reset();
  while((*theParticleIterator)())
  {
    G4ParticleDefinition* pd = theParticleIterator->value();
    G4ProcessManager* pm = pd->GetProcessManager();
    if(pm)
    {
      G4ProcessVector* pv = pm->GetAlongStepProcessVector(typeDoIt);
      G4VProcess* p = (*pv)[0];
      return ( (p->GetProcessName()) == "CoupledTransportation" );
    }
  }
  return false;
}

#include "G4ScoreSplittingProcess.hh"
void G4RunManagerKernel::SetScoreSplitter()
{
  G4ScoreSplittingProcess* pSplitter = new G4ScoreSplittingProcess();
  G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();
  G4ParticleTable::G4PTblDicIterator* theParticleIterator = theParticleTable->GetIterator();

  // Ensure that Process is added only once to the particles' process managers
  static G4ThreadLocal bool InitSplitter=false; 
  if( ! InitSplitter ) {
    InitSplitter = true;

    theParticleIterator->reset();
    while( (*theParticleIterator)() )  
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      if(pmanager)
	{ pmanager->AddDiscreteProcess(pSplitter); }
    }

    if(verboseLevel>0) 
    {
      G4cout << "G4RunManagerKernel -- G4ScoreSplittingProcess is appended to all particles." << G4endl;
    }
  }
}

void G4RunManagerKernel::SetupShadowProcess() const
{
    G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();
    G4ParticleTable::G4PTblDicIterator* theParticleIterator = theParticleTable->GetIterator();
    theParticleIterator->reset();
    //loop on particles and get process manager from there list of processes
    while((*theParticleIterator)())
    {
        G4ParticleDefinition* pd = theParticleIterator->value();
        G4ProcessManager* pm = pd->GetProcessManager();
        if(pm)
        {
            G4ProcessVector& procs = *(pm->GetProcessList());
            for ( G4int idx = 0 ; idx<procs.size() ; ++idx)
            {
                const G4VProcess* masterP = procs[idx]->GetMasterProcess();
                if ( ! masterP )
                {
                    //Process does not have an associated shadow master process
                    //We are in master mode or sequential
                    procs[idx]->SetMasterProcess(const_cast<G4VProcess*>(procs[idx]));
                }
            }
        }
    }
}
