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

#include "G4MTRunManagerKernel.hh"
#include "G4RegionStore.hh"
#include "G4StateManager.hh"
#include "G4AutoLock.hh"

std::vector<G4WorkerRunManager*>* G4MTRunManagerKernel::workerRMvector = 0;

namespace {
 G4Mutex workerRMMutex = G4MUTEX_INITIALIZER;
}

G4MTRunManagerKernel::G4MTRunManagerKernel() : G4RunManagerKernel(masterRMK)
{
    //This version of the constructor should never be called in sequential mode!
#ifndef G4MULTITHREADED
    G4ExceptionDescription msg;
    msg<<"Geant4 code is compiled without multi-threading support (-DG4MULTITHREADED is set to off).";
    msg<<" This type of RunManager can only be used in mult-threaded applications.";
    G4Exception("G4RunManagerKernel::G4RunManagerKernel()","Run0035",FatalException,msg);
#endif
    if(!workerRMvector) workerRMvector = new std::vector<G4WorkerRunManager*>;
    //Set flag that a MT-type kernel has been instantiated
    G4Threading::SetMultithreadedApplication(true);
}

G4MTRunManagerKernel::~G4MTRunManagerKernel()
{
    if(!workerRMvector)
    {
      if(workerRMvector->size()>0)
      {
        G4ExceptionDescription msg;
        msg<<"G4MTRunManagerKernel is to be deleted while "
           <<workerRMvector->size()<<" G4WorkerRunManager are still alive.";
        G4Exception("G4RunManagerKernel::~G4RunManagerKernel()",
                    "Run10035",FatalException,msg);
      }
      delete workerRMvector;
      workerRMvector = 0;
    }  
}

void G4MTRunManagerKernel::SetupShadowProcess() const
{
    //Behavior is the same as base class (sequential mode)
    //ShadowProcess pointer == process poitner
    G4RunManagerKernel::SetupShadowProcess();
}

#include "G4WorkerRunManager.hh"
#include "G4UserWorkerInitialization.hh"
#include "G4UserWorkerThreadInitialization.hh"
#include "G4VUserActionInitialization.hh"
#include "G4WorkerThread.hh"
#include "G4UImanager.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVReplica.hh"
#include "G4Region.hh"
#include "G4Material.hh"
#include "G4PhysicsVector.hh"
#include "G4VDecayChannel.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4MaterialTable.hh"
#include "G4PolyconeSide.hh"
#include "G4PolyhedraSide.hh"
#include "G4PVParameterised.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4VModularPhysicsList.hh"

G4ThreadLocal G4WorkerThread* G4MTRunManagerKernel::wThreadContext = 0;
G4WorkerThread* G4MTRunManagerKernel::GetWorkerThread() 
{ return wThreadContext; }

void* G4MTRunManagerKernel::StartThread(void* context)
{
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
//#ifdef G4MULTITHREADED
//    turnontpmalloc();
//#endif

  wThreadContext = static_cast<G4WorkerThread*>(context);  
  G4MTRunManager* masterRM = G4MTRunManager::GetMasterRunManager();

    
  //============================
  //Step-0: Thread ID
  //============================
  //Initliazie per-thread stream-output
  //The following line is needed before we actually do IO initialization
  //becasue the constructor of UI manager resets the IO destination.
  G4int thisID = wThreadContext->GetThreadId();
  G4Threading::G4SetThreadId(thisID);
  G4UImanager::GetUIpointer()->SetUpForAThread(thisID);

  //============================
  //Optimization Step
  //============================
  //Enforce thread affinity if requested
#if !defined(WIN32)
  if ( masterRM->GetPinAffinity() != 0 ) {
	  G4cout<<"AFFINITY SET"<<G4endl;
      //Assign this thread to cpus in a round robin way
	  G4int offset = masterRM->GetPinAffinity();
	  G4int cpuindex = 0;
	  if ( std::abs(offset)>G4Threading::G4GetNumberOfCores() ) {
		  G4Exception("G4MTRunManagerKernel::StarThread","Run0035",JustWarning,"Cannot set thread affinity, affinity parameter larger than number of cores");
	  }
	  if ( offset == 0 ) {
		  offset = 1;
		  G4Exception("G4MTRunManagerKernel::StarThread","Run0035",JustWarning,"Affinity parameter==0, using 1 instead.");
	  }
	  if (offset>0) { //Start assigning affinity to given CPU
		  --offset;
		  cpuindex = (thisID+offset) % G4Threading::G4GetNumberOfCores(); //Round robin
	  } else {//Exclude the given CPU
		  offset *= -1;
		  --offset;
		  G4int myidx = thisID%(G4Threading::G4GetNumberOfCores()-1);
		  cpuindex = myidx + (myidx>=offset);
	  }
	  G4cout<<"AFFINITY:"<<cpuindex<<G4endl;
      //Avoid compilation warning in C90 standard w/o MT
#if defined(G4MULTITHREADED)
      G4Thread t = G4THREADSELF();
#else
      G4Thread t;
#endif
      G4bool success = G4Threading::G4SetPinAffinity(cpuindex,t);
      if ( ! success ) {
          G4Exception("G4MTRunManagerKernel::StarThread","Run0035",JustWarning,"Cannot set thread affinity.");
      }
  }
#endif
    
  //============================
  //Step-1: Random number engine
  //============================
  //RNG Engine needs to be initialized by "cloning" the master one.
  const CLHEP::HepRandomEngine* masterEngine = masterRM->getMasterRandomEngine();
  masterRM->GetUserWorkerThreadInitialization()->SetupRNGEngine(masterEngine);

  //============================
  //Step-2: Initialize worker thread
  //============================
  if(masterRM->GetUserWorkerInitialization())
  { masterRM->GetUserWorkerInitialization()->WorkerInitialize(); }
  if(masterRM->GetUserActionInitialization())
  {
      G4VSteppingVerbose* sv = masterRM->GetUserActionInitialization()->InitializeSteppingVerbose();
      if ( sv ) { G4VSteppingVerbose::SetInstance(sv); }
  }
  //Now initialize worker part of shared objects (geometry/physics)
  wThreadContext->BuildGeometryAndPhysicsVector();
  G4WorkerRunManager* wrm
        = masterRM->GetUserWorkerThreadInitialization()->CreateWorkerRunManager();
  wrm->SetWorkerThread(wThreadContext);
  G4AutoLock wrmm(&workerRMMutex);
  workerRMvector->push_back(wrm);
  wrmm.unlock();

  //================================
  //Step-3: Setup worker run manager
  //================================
  // Set the detector and physics list to the worker thread. Share with master
  const G4VUserDetectorConstruction* detector = masterRM->GetUserDetectorConstruction();
  wrm->G4RunManager::SetUserInitialization(const_cast<G4VUserDetectorConstruction*>(detector));
  const G4VUserPhysicsList* physicslist = masterRM->GetUserPhysicsList();
  wrm->SetUserInitialization(const_cast<G4VUserPhysicsList*>(physicslist));

  //================================
  //Step-4: Initialize worker run manager
  //================================
  if(masterRM->GetUserActionInitialization())
  { masterRM->GetNonConstUserActionInitialization()->Build(); }
  if(masterRM->GetUserWorkerInitialization())
  { masterRM->GetUserWorkerInitialization()->WorkerStart(); }
  wrm->Initialize();

  //================================
  //Step5: Loop over requests from the master thread 
  //================================
  G4MTRunManager::WorkerActionRequest nextAction = masterRM->ThisWorkerWaitForNextAction();
  while( nextAction != G4MTRunManager::ENDWORKER )
  {
    if( nextAction == G4MTRunManager::NEXTITERATION ) // start the next run
    {
      //The following code deals with changing materials between runs
      static G4ThreadLocal G4bool skipInitialization = true;
      if(skipInitialization)
      { 
        // re-initialization is not necessary for the first run
        skipInitialization = false;
      }
      else
      {
//        ReinitializeGeometry();
          wThreadContext->UpdateGeometryAndPhysicsVectorFromMaster();          
      }

      // Execute UI commands stored in the masther UI manager
      std::vector<G4String> cmds = masterRM->GetCommandStack();
      G4UImanager* uimgr = G4UImanager::GetUIpointer(); //TLS instance
      std::vector<G4String>::const_iterator it = cmds.begin();
      for(;it!=cmds.end();it++)
      { uimgr->ApplyCommand(*it); }
      //Start this run
      G4int numevents = masterRM->GetNumberOfEventsToBeProcessed();
      G4String macroFile = masterRM->GetSelectMacro();
      G4int numSelect = masterRM->GetNumberOfSelectEvents();
      if ( macroFile == "" || macroFile == " " )
      {
          wrm->BeamOn(numevents);
      }
      else
      {
          wrm->BeamOn(numevents,macroFile,numSelect);
      }
    }
    else
    {
      G4ExceptionDescription d;
      d<<"Cannot continue, this worker has been requested an unknwon action: "
       <<nextAction<<" expecting: ENDWORKER(=" <<G4MTRunManager::ENDWORKER
       <<") or NEXTITERATION(="<<G4MTRunManager::NEXTITERATION<<")";
      G4Exception("G4MTRunManagerKernel::StartThread","Run0035",FatalException,d);
    }

    //Now wait for master thread to signal new action to be performed
    nextAction = masterRM->ThisWorkerWaitForNextAction();
  } //No more actions to perform

  //===============================
  //Step-6: Terminate worker thread
  //===============================
  if(masterRM->GetUserWorkerInitialization())
  { masterRM->GetUserWorkerInitialization()->WorkerStop(); }

  wrmm.lock();
  std::vector<G4WorkerRunManager*>::iterator itrWrm = workerRMvector->begin();
  for(;itrWrm!=workerRMvector->end();itrWrm++)
  {
    if((*itrWrm)==wrm)
    {
      workerRMvector->erase(itrWrm);
      break;
    }
  }
  wrmm.unlock();
  delete wrm;
  wThreadContext->DestroyGeometryAndPhysicsVector();
  wThreadContext = 0;

  return static_cast<void*>(0);
}

//Now moved to G4WorkerThread
//void G4MTRunManagerKernel::ReinitializeGeometry()
//{
//  G4AutoLock wrmm(&workerRMMutex);
//  //=================================================
//  //Step-0: keep sensitive detector and field manager
//  //=================================================
//  typedef std::map<G4LogicalVolume*,std::pair<G4VSensitiveDetector*,G4FieldManager*> > LV2SDFM;
//  LV2SDFM lvmap;
//  G4PhysicalVolumeStore* mphysVolStore = G4PhysicalVolumeStore::GetInstance(); 
//  for(size_t ip=0;ip<mphysVolStore->size();ip++)
//  {
//    G4VPhysicalVolume* pv = (*mphysVolStore)[ip];
//    G4LogicalVolume *lv = pv->GetLogicalVolume();
//    G4VSensitiveDetector* sd = lv->GetSensitiveDetector();
//    G4FieldManager* fm = lv->GetFieldManager();
//    if(sd||fm) lvmap[lv] = std::make_pair(sd,fm);
//  }
//
//  //===========================
//  //Step-1: Clean the instances
//  //===========================
//  const_cast<G4LVManager&>(G4LogicalVolume::GetSubInstanceManager()).FreeSlave();
//  const_cast<G4PVManager&>(G4VPhysicalVolume::GetSubInstanceManager()).FreeSlave();
//  const_cast<G4PVRManager&>(G4PVReplica::GetSubInstanceManager()).FreeSlave();
//  const_cast<G4RegionManager&>(G4Region::GetSubInstanceManager()).FreeSlave();
//  const_cast<G4PlSideManager&>(G4PolyconeSide::GetSubInstanceManager()).FreeSlave();
//  const_cast<G4PhSideManager&>(G4PolyhedraSide::GetSubInstanceManager()).FreeSlave();
//
//  //===========================
//  //Step-2: Re-create instances
//  //===========================
//  const_cast<G4LVManager&>(G4LogicalVolume::GetSubInstanceManager()).SlaveCopySubInstanceArray();
//  const_cast<G4PVManager&>(G4VPhysicalVolume::GetSubInstanceManager()).SlaveCopySubInstanceArray();
//  const_cast<G4PVRManager&>(G4PVReplica::GetSubInstanceManager()).SlaveCopySubInstanceArray();
//  const_cast<G4RegionManager&>(G4Region::GetSubInstanceManager()).SlaveInitializeSubInstance();
//  const_cast<G4PlSideManager&>(G4PolyconeSide::GetSubInstanceManager()).SlaveInitializeSubInstance();
//  const_cast<G4PhSideManager&>(G4PolyhedraSide::GetSubInstanceManager()).SlaveInitializeSubInstance();
//
//  //===============================
//  //Step-3: Re-initialize instances
//  //===============================
//  for(size_t ip=0;ip<mphysVolStore->size();ip++)
//  {
//    G4VPhysicalVolume* physVol = (*mphysVolStore)[ip];
//    G4LogicalVolume* g4LogicalVolume = physVol->GetLogicalVolume();
//    G4VSolid* g4VSolid = g4LogicalVolume->GetMasterSolid(); // shadow pointer
//    G4PVReplica* g4PVReplica = 0;
//    g4PVReplica =  dynamic_cast<G4PVReplica*>(physVol);
//    if(g4PVReplica) // if the volume is a replica
//    {
//      G4VSolid *slaveg4VSolid = g4VSolid->Clone();
//      g4LogicalVolume->InitialiseWorker(g4LogicalVolume,slaveg4VSolid,0);
//    }
//    else
//    { g4LogicalVolume->InitialiseWorker(g4LogicalVolume,g4VSolid,0); }
//  }
//
//  //===================================================
//  //Step-4: Restore sensitive detector and field manaer
//  //===================================================
//  LV2SDFM::const_iterator it = lvmap.begin();
//  for(; it!=lvmap.end() ; ++it )
//  {
//    G4LogicalVolume* lv = it->first;
//    G4VSensitiveDetector* sd = (it->second).first;
//    G4FieldManager* fm = (it->second).second;
//    lv->SetFieldManager(fm, false);
//    lv->SetSensitiveDetector(sd);
//  }
//  wrmm.unlock();
//}

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTableIterator.hh"
#include "G4DecayTable.hh"
#include "G4VDecayChannel.hh"

void G4MTRunManagerKernel::SetUpDecayChannels()
{
  G4ParticleTable::G4PTblDicIterator* pItr
    = G4ParticleTable::GetParticleTable()->GetIterator();
  pItr->reset();
  while((*pItr)())
  {
    G4DecayTable* dt = pItr->value()->GetDecayTable();
    if(dt)
    {
      G4int nCh = dt->entries();
      for(G4int i=0;i<nCh;i++)
      { dt->GetDecayChannel(i)->GetDaughter(0); }
    }
  }
}

void G4MTRunManagerKernel::BroadcastAbortRun(G4bool softAbort)
{
  G4AutoLock wrmm(&workerRMMutex);
  std::vector<G4WorkerRunManager*>::iterator itr = workerRMvector->begin();
  for(;itr!=workerRMvector->end();itr++)
  { (*itr)->AbortRun(softAbort); }
}

