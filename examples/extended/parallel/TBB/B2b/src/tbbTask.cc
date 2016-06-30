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
#include "tbbTask.hh"
#include "tbbWorkerRunManager.hh"
#include "G4Threading.hh"
#include "G4String.hh"
#include "G4WorkerRunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "G4UserWorkerThreadInitialization.hh"
#include "G4WorkerThread.hh"
#include "G4MTRunManagerKernel.hh"
#include "G4AutoLock.hh"
#include "G4UserWorkerInitialization.hh"
#include "G4VUserActionInitialization.hh"

//Equivalent to G4MTRunManagerKernel::StartThread

tbbTask::tbbTask(G4int anId,
                     tbb::concurrent_queue<const G4Run*>* out,
                     G4int nevts) 
: m_nEvents(nevts),
  m_taskID(anId),
  m_output(out),
  m_beamOnCondition(false)
{
}

tbbTask::~tbbTask()
{
}

#include <tbb/atomic.h>
namespace {
    tbb::atomic<int> counter;
}

tbb::task* tbbTask::execute()
{
  // In tbb we do not have anymore the concept of thread:
  // tasks run "somewhere"; but there is no control over where.
  // This somewhere is a thread, created, controlled and 'switched'
  // to our task by the TBB runtime system.
  // 
  // The "pedantic" way to proceed is: recreate the "context"
  // from scratch every time (context=local run manager, 
  // geometry, etc)
  // This would be a clear waste of resources: there is no 
  // need to initialize multiple times, when the same thread 
  // is used several times, ie for several tasks.

  // In the current version to avoid this we initialize only
  //   i) the first time we run on each thread, and
  //  ii) at the start of a new run on each thread (reinitialization).
  // We are confident that this works because TBB co-works with TLS.
  // In addition this ensures that minimal changes needed to G4 code base.
  //
  // Note 1: that this "thread" is responsible for 1 or more TBB task, 
  //         e.g. at least one event.
  // Note 2: In this first example, we do not care about memory usage:
  //         The resources required depend only on the total number of 
  //         threads used for at least one event during the simulation.
  // Note 3: It is possible to do better - this is left for another example.
  //         Here is a sketch how:
  //         If there is a set maximum <N> threads which will be doing 
  //         simulation at any point, and that this is smaller than the 
  //         number of TBB tasks which can run simulataneously. 
  //         What would need to be done is to "acquire" a resource (workspace) 
  //         which holds all the memory required for a running task/thread,
  //         and release it at the end of the task to be re-used by other tasks
  //         possibly on different threads. This is demonstrated in another 
  //         example.

  static G4ThreadLocal tbbWorkerRunManager* localRM = 0;

  
  //Is this task running for the first time?
  //How to re-initialize between runs????
  if (! localRM ) {
    //It's a new thread, basically repeat what is being done in 
    //G4MTRunManagerKernel::StarThread with an 
    //important difference, do not process data or wait for master to 
    //communicate what to do, it will
    //not happen!
    G4MTRunManager* masterRM = G4MTRunManager::GetMasterRunManager();
    //============================
    //Step-0: Thread ID
    //============================
    //Initliazie per-thread stream-output
    //The following line is needed before we actually do I/O initialization
    //because the constructor of UI manager resets the I/O destination.
      G4int thisId = counter.fetch_and_increment();
      G4Threading::G4SetThreadId( thisId );
      G4UImanager::GetUIpointer()->SetUpForAThread( thisId );
  
    //============================
    //Step-1: Random number engine
    //============================
    //RNG Engine needs to be initialized by "cloning" the master one.
    const CLHEP::HepRandomEngine* masterEngine =
          masterRM->getMasterRandomEngine();
    masterRM->GetUserWorkerThreadInitialization()->SetupRNGEngine(masterEngine);
  
    //============================
    //Step-2: Initialize worker thread
    //============================
    if(masterRM->GetUserWorkerInitialization())
      masterRM->GetUserWorkerInitialization()->WorkerInitialize();
    if(masterRM->GetUserActionInitialization()) {
        G4VSteppingVerbose* sv =
        masterRM->GetUserActionInitialization()->InitializeSteppingVerbose();
        if (sv) G4VSteppingVerbose::SetInstance(sv);
    }
    //Now initialize worker part of shared objects (geometry/physics)
    G4WorkerThread::BuildGeometryAndPhysicsVector();
    localRM = static_cast<tbbWorkerRunManager*>(
                    masterRM->GetUserWorkerThreadInitialization()->CreateWorkerRunManager()
                                                );
      G4cout<<localRM<<G4endl;
    //localRM->SetWorkerThread(wThreadContext);
//    G4AutoLock wrmm(&workerRMMutex);
//    G4MTRunManagerKernel::workerRMvector->push_back(localRM); //<<<<?????? ANDREA TBB
//    wrmm.unlock();

    //================================
    //Step-3: Setup worker run manager
    //================================
    // Set the detector and physics list to the worker thread. Share with master
    const G4VUserDetectorConstruction* detector = 
        masterRM->GetUserDetectorConstruction();
    localRM->G4RunManager::SetUserInitialization(
        const_cast<G4VUserDetectorConstruction*>(detector));
    const G4VUserPhysicsList* physicslist = masterRM->GetUserPhysicsList();
    localRM->SetUserInitialization(const_cast<G4VUserPhysicsList*>(physicslist));
    
    //================================
    //Step-4: Initialize worker run manager
    //================================
    if(masterRM->GetUserActionInitialization())
      { masterRM->GetNonConstUserActionInitialization()->Build(); }
    if(masterRM->GetUserWorkerInitialization())
      { masterRM->GetUserWorkerInitialization()->WorkerStart(); }
    localRM->Initialize();
    
    //Problem at this point is if there is more than one run...
    // Execute UI commands stored in the masther UI manager
    std::vector<G4String> cmds = masterRM->GetCommandStack();
    G4UImanager* uimgr = G4UImanager::GetUIpointer(); //TLS instance
    std::vector<G4String>::const_iterator it = cmds.begin();
    for(;it!=cmds.end();it++)
    { uimgr->ApplyCommand(*it); }
    //Start this run
    G4String macroFile = masterRM->GetSelectMacro();
    G4int numSelect = masterRM->GetNumberOfSelectEvents();
    if ( macroFile == "" || macroFile == " " )
    {
        localRM->BeamOn(m_nEvents,0,numSelect);
    }
    else
    {
        localRM->BeamOn(m_nEvents,macroFile,numSelect);
    }
    //======= NEW TBB SPECIFIC ===========
    //Step-5: Initialize and start run
    //====================================
    // bla-bla-bla-bla
    // This is basically BeamOn 
    m_beamOnCondition = localRM->ConfirmBeamOnCondition();
    if (m_beamOnCondition) {
      localRM->SetNumberOfEventsToBeProcessed( m_nEvents );
      localRM->ConstructScoringWorlds(); 
      localRM->RunInitialization();
      //Register this G4Run in output
      //Note: the idea is that we are going to accumulate everything in 
      //G4Run or derived class. We let the kernel know this is the object
      //where the output is accumulated for the tbb::tasks that run on 
      //this thread.
      if ( m_output ) {}
    }
  }
  assert(localRM!=0);
  if ( m_beamOnCondition ) {
    localRM->DoEventLoop( m_nEvents );
    //localRM->RunTermination(); //<<<< How to call this??? ANDREA TBB
  }

  //Cannot be done here since thread can be re-used by other task
  //===============================
  //Step-6: Terminate worker thread
  //===============================
  // if(masterRM->GetUserWorkerInitialization())
  // { masterRM->GetUserWorkerInitialization()->WorkerStop(); }

  // wrmm.lock();
// std::vector<G4WorkerRunManager*>::iterator itrWrm = workerRMvector->begin();
  // for(;itrWrm!=workerRMvector->end();itrWrm++)
  // {
  //   if((*itrWrm)==wrm)
  //   {
  //     workerRMvector->erase(itrWrm);
  //     break;
  //   }
  // }
  // wrmm.unlock();
    
  // wThreadContext->DestroyGeometryAndPhysicsVector();
  // wThreadContext = 0;

  return static_cast<tbb::task*>(NULL);
}
