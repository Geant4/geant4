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

#include "G4MTRunManager.hh"
#include "G4MTRunManagerKernel.hh"
#include "G4Timer.hh"
#include "G4StateManager.hh"
#include "G4ScoringManager.hh"
#include "G4TransportationManager.hh"
#include "G4VUserActionInitialization.hh"
#include "G4UserWorkerInitialization.hh"
#include "G4UserWorkerThreadInitialization.hh"
#include "G4WorkerThread.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4AutoLock.hh"
#include "G4WorkerRunManager.hh"
#include "G4UserRunAction.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Timer.hh"

G4ScoringManager* G4MTRunManager::masterScM = 0;
G4MTRunManager::masterWorlds_t G4MTRunManager::masterWorlds = G4MTRunManager::masterWorlds_t();
G4MTRunManager* G4MTRunManager::fMasterRM = 0;
G4int G4MTRunManager::seedOncePerCommunication = 0;

namespace {
 G4Mutex cmdHandlingMutex = G4MUTEX_INITIALIZER;
 G4Mutex scorerMergerMutex = G4MUTEX_INITIALIZER;
 G4Mutex runMergerMutex = G4MUTEX_INITIALIZER;
 G4Mutex setUpEventMutex = G4MUTEX_INITIALIZER;
}

G4MTRunManager* G4MTRunManager::GetMasterRunManager()
{
////////#ifdef G4MULTITHREADED
    return fMasterRM;
////////#else
////////    return G4RunManager::GetRunManager();
////////#endif
}

G4RunManagerKernel* G4MTRunManager::GetMasterRunManagerKernel()
{
    return fMasterRM->kernel;
}

G4MTRunManagerKernel* G4MTRunManager::GetMTMasterRunManagerKernel()
{
    return fMasterRM->MTkernel;
}

G4MTRunManager::G4MTRunManager() : G4RunManager(masterRM),
    nworkers(2),forcedNwokers(-1),pinAffinity(0),
    masterRNGEngine(0),
    nextActionRequest(WorkerActionRequest::UNDEFINED),
    eventModuloDef(0),eventModulo(1),
    nSeedsUsed(0),nSeedsFilled(0),
    nSeedsMax(10000),nSeedsPerEvent(2)
{
    if ( fMasterRM )
    {
        G4Exception("G4MTRunManager::G4MTRunManager", "Run0110",FatalException,
                    "Another instance of a G4MTRunManager already exists.");
    }
    fMasterRM = this;
    MTkernel = static_cast<G4MTRunManagerKernel*>(kernel);
#ifndef G4MULTITHREADED
    G4ExceptionDescription msg;
    msg << "Geant4 code is compiled without multi-threading support"
        << "(-DG4MULTITHREADED is set to off).\n";
    msg << "G4MTRunManager can only be used in multi-threaded applications.";
    G4Exception("G4MTRunManager::G4MTRunManager","Run0111",FatalException,msg);
#endif

    G4int numberOfStaticAllocators = kernel->GetNumberOfStaticAllocators();
    if(numberOfStaticAllocators>0)
    {
      G4ExceptionDescription msg1;
      msg1 << "There are " << numberOfStaticAllocators
           << " static G4Allocator objects detected.\n"
           << "In multi-threaded mode, all G4Allocator objects must be dynamicly instantiated.";
      G4Exception("G4MTRunManager::G4MTRunManager","Run1035",FatalException,msg1);
    }
    G4UImanager::GetUIpointer()->SetMasterUIManager(true);
    masterScM = G4ScoringManager::GetScoringManagerIfExist();

    //Check if a default RandomNumberGenerator has been created by user,
    // if not create default one
    //Note this call forces creation of defaults if not already there
    //G4Random::getTheEngine(); //User did not specify RNG, create defaults
    //Now remember the master instance of the RNG Engine
    masterRNGEngine = G4Random::getTheEngine();

    numberOfEventToBeProcessed = 0;
    randDbl = new double[nSeedsPerEvent*nSeedsMax];

    char* env = getenv("G4FORCENUMBEROFTHREADS");
    if(env)
    {
      G4String envS = env;
      if(envS=="MAX"||envS=="max")
      { forcedNwokers = G4Threading::G4GetNumberOfCores(); }
      else
      {
        std::istringstream is(env);
        G4int val = -1;
        is >> val;
        if(val>0)
        { forcedNwokers = val; }
        else
        {
          G4ExceptionDescription msg2;
          msg2 << "Environment variable G4FORCENUMBEROFTHREADS has an invalid value <"
               << envS << ">. It has to be an integer or a word \"max\".\n"
               << "G4FORCENUMBEROFTHREADS is ignored.";
          G4Exception("G4MTRunManager::G4MTRunManager","Run1039",JustWarning,msg2);
        }
      }
      if(forcedNwokers>0)
      {
        nworkers = forcedNwokers;
        G4cout << "### Number of threads is forced to " << forcedNwokers 
               << " by Environment variable G4FORCENUMBEROFTHREADS." << G4endl;
      }
    }
}

G4MTRunManager::~G4MTRunManager()
{
    //TODO: Currently does not work due to concurrent deletion of something
    //      that is shared:
    //G4ProcessTable::DeleteMessenger from ~G4RunManager
    //G4cout<<"Destroy MTRunManager"<<G4endl;//ANDREA
    TerminateWorkers();
    delete [] randDbl;
}

void G4MTRunManager::StoreRNGStatus(const G4String& fn )
{
    std::ostringstream os;
    os << randomNumberStatusDir << "G4Master_"<<fn <<".rndm";
    G4Random::saveEngineStatus(os.str().c_str());
}

void G4MTRunManager::SetNumberOfThreads(G4int n )
{
    if ( threads.size() != 0 )
    {
      G4ExceptionDescription msg;
      msg << "Number of threads cannot be changed at this moment \n"
          << "(old threads are still alive). Method ignored.";
      G4Exception("G4MTRunManager::SetNumberOfThreads(G4int)",
                  "Run0112", JustWarning, msg);
    }
    else if ( forcedNwokers > 0 )
    {
      G4ExceptionDescription msg;
      msg << "Number of threads is forced to " << forcedNwokers
          << " by G4FORCENUMBEROFTHREADS shell variable.\n"
          << "Method ignored.";
      G4Exception("G4MTRunManager::SetNumberOfThreads(G4int)",
                  "Run0113", JustWarning, msg);
    }
    else  
    {
      nworkers = n;
    }
}

void G4MTRunManager::Initialize()
{
    G4RunManager::Initialize();

    // make sure all worker threads are set up.
    BeamOn(0);
    SetRunIDCounter(0);
    ///G4UImanager::GetUIpointer()->SetIgnoreCmdNotFound(true);
}

////void G4MTRunManager::TerminateEventLoop()
////{
////    //Nothing to do
////}
void G4MTRunManager::ProcessOneEvent(G4int)
{
    //Nothing to do
}
void G4MTRunManager::TerminateOneEvent()
{
    //Nothing to do
}

void G4MTRunManager::PrepareCommandsStack() {
    G4AutoLock l(&cmdHandlingMutex);
    uiCmdsForWorkers.clear();
    std::vector<G4String>* cmdCopy = G4UImanager::GetUIpointer()->GetCommandStack();
    for ( std::vector<G4String>::const_iterator it = cmdCopy->begin() ;
         it != cmdCopy->end(); ++it )
        uiCmdsForWorkers.push_back(*it);
    cmdCopy->clear();
    delete cmdCopy;
}

std::vector<G4String> G4MTRunManager::GetCommandStack()
{
    G4AutoLock l(&cmdHandlingMutex);
    return uiCmdsForWorkers;
}

void G4MTRunManager::CreateAndStartWorkers()
{
    //Now loop on requested number of workers
    //This will also start the workers
    //Currently we do not allow to change the
    //number of threads: threads area created once
    if ( threads.size() == 0 ) {
        for ( G4int nw = 0 ; nw<nworkers; ++nw) {
            //Create a new worker and remember it
            G4WorkerThread* context = new G4WorkerThread;
            context->SetNumberThreads(nworkers);
            context->SetThreadId(nw);
            G4Thread* thread = userWorkerThreadInitialization->CreateAndStartWorker(context);
            threads.push_back(thread);
        }
    }
    //Signal to threads they can start a new run
    NewActionRequest(WorkerActionRequest::NEXTITERATION);
}


void G4MTRunManager::InitializeEventLoop(G4int n_event, const char* macroFile, G4int n_select)
{
  MTkernel->SetUpDecayChannels();
  numberOfEventToBeProcessed = n_event;
  numberOfEventProcessed = 0;

  if(!fakeRun)
  {
    nSeedsUsed = 0;
    nSeedsFilled = 0;

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

    //initialize seeds
    //If user did not implement InitializeSeeds,
    // use default: nSeedsPerEvent seeds per event 
    if( eventModuloDef > 0 )
    {
      eventModulo = eventModuloDef;
      if(eventModulo > numberOfEventToBeProcessed/nworkers)
      {
        eventModulo = numberOfEventToBeProcessed/nworkers;
        if(eventModulo<1) eventModulo =1;
        G4ExceptionDescription msgd;
        msgd << "Event modulo is reduced to " << eventModulo
            << " to distribute events to all threads.";
        G4Exception("G4MTRunManager::InitializeEventLoop()",
                  "Run10035", JustWarning, msgd);
      }
    }
    else
    {
      eventModulo = int(std::sqrt(double(numberOfEventToBeProcessed/nworkers)));
      if(eventModulo<1) eventModulo =1;
    }
    if ( InitializeSeeds(n_event) == false && n_event>0 )
    {
        G4RNGHelper* helper = G4RNGHelper::GetInstance();
        switch(seedOncePerCommunication)
        {
         case 0:
          nSeedsFilled = n_event;
          break;
         case 1:
          nSeedsFilled = nworkers;
          break;
         case 2:
          nSeedsFilled = n_event/eventModulo + 1;
          break;
         default:
          G4ExceptionDescription msgd;
          msgd << "Parameter value <" << seedOncePerCommunication
               << "> of seedOncePerCommunication is invalid. It is reset to 0." ;
          G4Exception("G4MTRunManager::InitializeEventLoop()",
                  "Run10036", JustWarning, msgd);
          seedOncePerCommunication = 0;
          nSeedsFilled = n_event;
        }

        // Generates up to nSeedsMax seed pairs only.
        if(nSeedsFilled>nSeedsMax) nSeedsFilled=nSeedsMax;
        masterRNGEngine->flatArray(nSeedsPerEvent*nSeedsFilled,randDbl); 
        helper->Fill(randDbl,nSeedsFilled,n_event,nSeedsPerEvent);
    }
  }
  
  //Now initialize workers. Check if user defined a WorkerThreadInitialization
  if ( userWorkerThreadInitialization == 0 )
  { userWorkerThreadInitialization = new G4UserWorkerThreadInitialization(); }
    
  //Prepare UI commands for threads
  PrepareCommandsStack();

  //Start worker threads
  CreateAndStartWorkers();
    
  // We need a barrier here. Wait for workers to start event loop.
  //This will return only when all workers have started processing events.
  WaitForReadyWorkers();
}

void G4MTRunManager::RefillSeeds()
{
  G4RNGHelper* helper = G4RNGHelper::GetInstance();
  G4int nFill = 0;
  switch(seedOncePerCommunication)
  {
   case 0:
    nFill = numberOfEventToBeProcessed - nSeedsFilled;
    break;
   case 1:
    nFill = nworkers - nSeedsFilled;
    break;
   case 2:
   default:
    nFill = (numberOfEventToBeProcessed - nSeedsFilled*eventModulo)/eventModulo + 1;
  }
  // Generates up to nSeedsMax seed pairs only.
  if(nFill>nSeedsMax) nFill=nSeedsMax;
  masterRNGEngine->flatArray(nSeedsPerEvent*nFill,randDbl); 
  helper->Refill(randDbl,nFill);
  nSeedsFilled += nFill;
//G4cout<<"helper->Refill() for "<<nFill<<" events."<<G4endl;
}
    
void G4MTRunManager::RunTermination()
{
  //Wait for all worker threads to have finished the run
  //i.e. wait for them to return from RunTermination()
  //This guarantee that userrunaction for workers has been called

  // Wait now for all threads to finish event-loop
  WaitForEndEventLoopWorkers();
  //Now call base-class methof
  G4RunManager::TerminateEventLoop();
  G4RunManager::RunTermination();
}

void G4MTRunManager::ConstructScoringWorlds()
{
    masterScM = G4ScoringManager::GetScoringManagerIfExist();
    //Call base class stuff...
    G4RunManager::ConstructScoringWorlds();
    
    masterWorlds.clear();
    size_t nWorlds = G4TransportationManager::GetTransportationManager()->GetNoWorlds();
    std::vector<G4VPhysicalVolume*>::iterator itrW
                   = G4TransportationManager::GetTransportationManager()->GetWorldsIterator();
    for(size_t iWorld=0;iWorld<nWorlds;iWorld++)
    {
      addWorld(iWorld,*itrW); 
      itrW++;
    }
}

void G4MTRunManager::SetUserInitialization(G4UserWorkerInitialization* userInit)
{
  userWorkerInitialization = userInit;
}

void G4MTRunManager::SetUserInitialization(G4UserWorkerThreadInitialization* userInit)
{
  userWorkerThreadInitialization = userInit;
}

void G4MTRunManager::SetUserInitialization(G4VUserActionInitialization* userInit)
{
  userActionInitialization = userInit;
  userActionInitialization->BuildForMaster();
}

void G4MTRunManager::SetUserInitialization(G4VUserPhysicsList *userPL)
{
    G4RunManager::SetUserInitialization(userPL);
    //Needed for MT, to be moved in kernel 
}

void G4MTRunManager::SetUserInitialization(G4VUserDetectorConstruction *userDC)
{
    G4RunManager::SetUserInitialization(userDC);
}

void G4MTRunManager::SetUserAction(G4UserRunAction* userAction)
{
    G4RunManager::SetUserAction(userAction); 
    if(userAction) userAction->SetMaster();
}

void G4MTRunManager::SetUserAction(G4VUserPrimaryGeneratorAction* /*userAction*/)
{
  G4Exception("G4MTRunManager::SetUserAction()", "Run0123", FatalException,
    "For multi-threaded version, define G4VUserPrimaryGeneratorAction in G4VUserActionInitialization.");
}

void G4MTRunManager::SetUserAction(G4UserEventAction* /*userAction*/)
{
  G4Exception("G4MTRunManager::SetUserAction()", "Run0124", FatalException,
    "For multi-threaded version, define G4UserEventAction in G4VUserActionInitialization.");
}

void G4MTRunManager::SetUserAction(G4UserStackingAction* /*userAction*/)
{
  G4Exception("G4MTRunManager::SetUserAction()", "Run0125", FatalException,
    "For multi-threaded version, define G4UserStackingAction in G4VUserActionInitialization.");
}

void G4MTRunManager::SetUserAction(G4UserTrackingAction* /*userAction*/)
{
  G4Exception("G4MTRunManager::SetUserAction()", "Run0126", FatalException,
    "For multi-threaded version, define G4UserTrackingAction in G4VUserActionInitialization.");
}

void G4MTRunManager::SetUserAction(G4UserSteppingAction* /*userAction*/)
{
  G4Exception("G4MTRunManager::SetUserAction()", "Run0127", FatalException,
    "For multi-threaded version, define G4UserSteppingAction in G4VUserActionInitialization.");
}

void G4MTRunManager::MergeScores(const G4ScoringManager* localScoringManager)
{
  G4AutoLock l(&scorerMergerMutex);
  if(masterScM) masterScM->Merge(localScoringManager); 
}

void G4MTRunManager::MergeRun(const G4Run* localRun)
{
  G4AutoLock l(&runMergerMutex);
  if(currentRun) currentRun->Merge(localRun); 
}

G4bool G4MTRunManager::SetUpAnEvent(G4Event* evt,long& s1,long& s2,long& s3,G4bool reseedRequired)
{
  G4AutoLock l(&setUpEventMutex);
  if( numberOfEventProcessed < numberOfEventToBeProcessed )
  {
    evt->SetEventID(numberOfEventProcessed);
    if(reseedRequired)
    {
      G4RNGHelper* helper = G4RNGHelper::GetInstance();
      G4int idx_rndm = nSeedsPerEvent*nSeedsUsed;
      s1 = helper->GetSeed(idx_rndm);
      s2 = helper->GetSeed(idx_rndm+1);
      if(nSeedsPerEvent==3) s3 = helper->GetSeed(idx_rndm+2);
      nSeedsUsed++;
      if(nSeedsUsed==nSeedsFilled) RefillSeeds();
    }
    numberOfEventProcessed++;
    return true;
  }
  return false;
}

G4int G4MTRunManager::SetUpNEvents(G4Event* evt, G4SeedsQueue* seedsQueue,G4bool reseedRequired)
{
  G4AutoLock l(&setUpEventMutex);
  if( numberOfEventProcessed < numberOfEventToBeProcessed && !runAborted )
  {
    G4int nev = eventModulo;
    if(numberOfEventProcessed + nev > numberOfEventToBeProcessed)
    { nev = numberOfEventToBeProcessed - numberOfEventProcessed; }
    evt->SetEventID(numberOfEventProcessed);
    if(reseedRequired)
    {
      G4RNGHelper* helper = G4RNGHelper::GetInstance();
      G4int nevRnd = nev;
      if(seedOncePerCommunication>0) nevRnd = 1;
      for(int i=0;i<nevRnd;i++)
      {
        seedsQueue->push(helper->GetSeed(nSeedsPerEvent*nSeedsUsed));
        seedsQueue->push(helper->GetSeed(nSeedsPerEvent*nSeedsUsed+1));
        if(nSeedsPerEvent==3)
          seedsQueue->push(helper->GetSeed(nSeedsPerEvent*nSeedsUsed+2));
        nSeedsUsed++;
        if(nSeedsUsed==nSeedsFilled) RefillSeeds();
      }
    }
    numberOfEventProcessed += nev;
    return nev;
  }
  return 0;
}

void G4MTRunManager::TerminateWorkers()
{
    //Force workers to execute (if any) all UI commands left in the stack
    RequestWorkersProcessCommandsStack();
    //Ask workers to exit
    NewActionRequest( WorkerActionRequest::ENDWORKER );
    //Now join threads.
#ifdef G4MULTITHREADED //protect here to prevent warning in compilation
    while ( ! threads.empty() )
    {
        G4Thread* t = * ( threads.begin() );
        threads.pop_front();
        userWorkerThreadInitialization->JoinWorker(t);
        //G4THREADJOIN(*t);
        delete t;
    }
#endif
    threads.clear();
}

void G4MTRunManager::AbortRun(G4bool softAbort)
{
  // This method is valid only for GeomClosed or EventProc state
  G4ApplicationState currentState =
    G4StateManager::GetStateManager()->GetCurrentState();
  if(currentState==G4State_GeomClosed || currentState==G4State_EventProc)
  {
    runAborted = true;
    MTkernel->BroadcastAbortRun(softAbort);
  }
  else
  {
    G4cerr << "Run is not in progress. AbortRun() ignored." << G4endl;
  }
}

void G4MTRunManager::AbortEvent()
{
  // nothing to do in the master thread
}

void G4MTRunManager::WaitForReadyWorkers() {
  beginOfEventLoopBarrier.Wait( GetNumberActiveThreads() );
  endOfEventLoopBarrier.ResetCounter();
  beginOfEventLoopBarrier.ReleaseBarrier();
}

void G4MTRunManager::ThisWorkerReady() {
  beginOfEventLoopBarrier.ThisWorkerReady();
}

void G4MTRunManager::WaitForEndEventLoopWorkers() {
  endOfEventLoopBarrier.Wait( GetNumberActiveThreads() );
  beginOfEventLoopBarrier.ResetCounter();
  endOfEventLoopBarrier.ReleaseBarrier();
}

void G4MTRunManager::ThisWorkerEndEventLoop() {
  endOfEventLoopBarrier.ThisWorkerReady();
}

void G4MTRunManager::NewActionRequest(G4MTRunManager::WorkerActionRequest newRequest) {
   nextActionRequestBarrier.Wait( GetNumberActiveThreads() );
   //nextActionRequest is a shared resource, but there is no
   //data-race thanks to the barrier: all threads are waiting
   nextActionRequest = newRequest;
   nextActionRequestBarrier.ReleaseBarrier();
}

G4MTRunManager::WorkerActionRequest G4MTRunManager::ThisWorkerWaitForNextAction() {
  nextActionRequestBarrier.ThisWorkerReady();
  return nextActionRequest;
}

void G4MTRunManager::RequestWorkersProcessCommandsStack() {
  PrepareCommandsStack();
  NewActionRequest(WorkerActionRequest::PROCESSUI);
  processUIBarrier.SetActiveThreads( GetNumberActiveThreads() );
  processUIBarrier.WaitForReadyWorkers();
}

void G4MTRunManager::ThisWorkerProcessCommandsStackDone() {
  processUIBarrier.ThisWorkerReady();
}
void G4MTRunManager::SetPinAffinity(G4int n)
{
	if ( n == 0 )
	{
		G4Exception("G4MTRunManager::SetPinAffinity",
					"Run0114",FatalException,
					"Pin affinity must be >0 or <0.");
	}
	pinAffinity = n;
	return;
}
