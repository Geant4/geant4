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
#include "G4Timer.hh"
#include "G4ScoringManager.hh"
#include "G4TransportationManager.hh"
#include "G4VUserActionInitialization.hh"
#include "G4UserWorkerInitialization.hh"
#include "G4WorkerThread.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4AutoLock.hh"
#include "G4WorkerRunManager.hh"
#include "G4UserRunAction.hh"
#include "G4Timer.hh"

long* G4MTRunManager::seeds = 0;
G4int G4MTRunManager::seedsnum = 0;
G4ScoringManager* G4MTRunManager::masterScM = 0;
G4MTRunManager::masterWorlds_t G4MTRunManager::masterWorlds = G4MTRunManager::masterWorlds_t();
G4MTRunManager* G4MTRunManager::masterRM = 0;

namespace {
 G4Mutex workesRMMutex = G4MUTEX_INITIALIZER;
 //Mutex to access/manipulate workersRM
 G4Mutex cmdHandlingMutex = G4MUTEX_INITIALIZER;
 //Mutex to access/manipulate uiCmdsForWorker command stack
 G4Mutex scorerMergerMutex = G4MUTEX_INITIALIZER;
 G4Mutex runMergerMutex = G4MUTEX_INITIALIZER;
}

//This is needed to initialize windows conditions
#if defined(WIN32)
namespace {
	void InitializeWindowsConditions();
}
#endif

G4MTRunManager* G4MTRunManager::GetMasterRunManager()
{
    return masterRM;
}

G4RunManagerKernel* G4MTRunManager::GetMasterRunManagerKernel()
{
    return masterRM->kernel;
}

G4MTRunManager::G4MTRunManager() : G4RunManager(false) ,
    nworkers(2),
    masterRNGEngine(0),
    workerG4coutFileName(""),
    workerG4cerrFileName(""),
    workerG4coutAppendFlag(false),
    workerG4cerrAppendFlag(false)
{
    if ( masterRM )
    {
        G4Exception("G4MTRunManager::G4MTRunManager", "Run0035",FatalException,
                    "Another instance of a G4MTRunManager already exists.");
    }
    masterRM = this;
#ifndef G4MULTITHREADED
    G4ExceptionDescription msg;
    msg<<"Geant4 code is compiled without multi-threading support (-DG4MULTITHREADED is set to off).";
    msg<<" This type of RunManager can only be used in mult-threaded applications.";
    G4Exception("G4MTRunManager::G4MTRunManager","Run0035",FatalException,msg);
#endif

    G4UImanager::GetUIpointer()->SetMasterUIManager(true);
    masterScM = G4ScoringManager::GetScoringManagerIfExist();

    //Check if a default RandomNumberGenerator has been created by user,
    // if not create default one
    //Note this call forces creation of defaults if not already there
    //G4Random::getTheEngine(); //User did not specify RNG, create defaults
    //Now remember the master instance of the RNG Engine
    masterRNGEngine = G4Random::getTheEngine();
#if defined (WIN32)
    InitializeWindowsConditions();
#endif
}

G4MTRunManager::~G4MTRunManager()
{
    //TODO: Currently does not work due to concurrent deletion of something that is shared:
    //G4ProcessTable::DeleteMessenger from ~G4RunManager
    //G4cout<<"Destroy MTRunManager"<<G4endl;//ANDREA
    TerminateWorkers();
    delete[] seeds;
    seeds = 0;
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
        G4Exception("G4MTRunManager::SetNumberOfThreads(G4int)","Run0035",JustWarning,"Limitation: Number of threads cannot be changed at this moment (old threads are still alive).");
    }
    else
    {
        nworkers = n;
    }
}

void G4MTRunManager::TerminateEventLoop()
{
    //Nothing to do
}
void G4MTRunManager::ProcessOneEvent(G4int)
{
    //Nothing to do
}
void G4MTRunManager::TerminateOneEvent()
{
    //Nothing to do
}

long G4MTRunManager::GetSeed(G4int i)
{
    if ( i < seedsnum) return seeds[i];
    G4ExceptionDescription msg;
    msg << "Not enough random seeds, requested seed number "<<i<<" but "<<seedsnum<<" seeds available (i>=seedsnum)";
    G4Exception("G4MTRunManager::GetSeed","Run0035", FatalException,msg);
    return -LONG_MAX;
}

void G4MTRunManager::PrepareCommandsStack() {
    G4AutoLock l(&cmdHandlingMutex);
    uiCmdsForWorkers.clear();
    std::vector<G4String>* cmdCopy = G4UImanager::GetUIpointer()->GetCommandStack();
    for ( std::vector<G4String>::const_iterator it = cmdCopy->begin() ;
         it != cmdCopy->end(); ++it )
        uiCmdsForWorkers.push_back(*it);
    delete cmdCopy;
}

std::vector<G4String> G4MTRunManager::GetCommandStack()
{
    G4AutoLock l(&cmdHandlingMutex);
    return uiCmdsForWorkers;
}

void G4MTRunManager::InitializeEventLoop(G4int n_events, const char* macroFile, G4int n_select)
{
    if(verboseLevel>0)
    { timer->Start(); }
    
    n_select_msg = n_select;
    if(macroFile!=0)
    {
        if(n_select_msg<0) n_select_msg = n_events;
        msgText = "/control/execute ";
        msgText += macroFile;
    }
    else
    { n_select_msg = -1; }

    //G4UImanager::GetUIpointer()->IgnoreCmdNotFound(false);
    //User initialize seeds
    
    //If user did not implement InitializeSeeds, use default: 2 seeds per event number
    if ( InitializeSeeds(n_events) == false )
    {
        InitializeSeedsQueue(n_events*2);
        for ( G4int ne = 0 ; ne < n_events*2 ; ++ne)
            AddOneSeed( (long) (100000000L * G4Random::getTheGenerator()->flat()) );
    }
    
    //Now initialize workers. Check if user defined a WorkerInitialization
    if ( userWorkerInitialization == 0 )
    {
////        G4Exception("G4MTRunManager::InitializeEventLoop(G4int,const char*,G4int)",
////                    "Run0035",FatalException,"No G4VUserWorkerInitialization found");
      userWorkerInitialization = new G4UserWorkerInitialization();
    }
    
    //Prepare UI commands for threads
    PrepareCommandsStack();
    //Now loop on requested number of workers
    //This will also start the workers
    //Currently we do not allow to change the
    //number of threads: threads area created once
    //and that's all
    if ( threads.size() == 0 ) {
    for ( G4int nw = 0 ; nw<nworkers; ++nw) {
        //Create a new worker and remember it
        G4WorkerThread* context = new G4WorkerThread;
        context->SetNumberThreads(nworkers);
        context->SetNumberEvents(n_events);
        context->SetThreadId(nw);
        if (! workerG4coutFileName.empty() )
        {
            context->SetOutputFileName(workerG4coutFileName,workerG4coutAppendFlag);
        }
        if (! workerG4cerrFileName.empty() )
        {
            context->SetOutputErrFileName(workerG4cerrFileName,workerG4cerrAppendFlag);
        }
        context->SetOutputUseBuffer( workerG4coutcerrBufferFlag );
        G4Thread* thread = userWorkerInitialization->CreateAndStartWorker(context);
        threads.push_back(thread);
    }
    }
    //Signal to threads they can start a new run
    NewActionRequest(NEXTITERATION);
    
    // We need a barrier here. Wait for workers to start event loop.
    //This will return only when all workers have started processing events.
    WaitForReadyWorkers();
}

void G4MTRunManager::RunTermination()
{
    //Wait for all worker threads to have finished the run
    //i.e. wait for them to return from RunTermination()
    //This guarantee that userrunaction for workers has been called

    // Wait now for all threads to finish event-loop
    WaitForEndEventLoopWorkers();
    //Now call base-class methof
    G4RunManager::RunTermination();
}

void G4MTRunManager::AddOneSeed( long seed )
{
    if ( !seeds )
    {
        G4ExceptionDescription msg;
        msg << "Cannot add a new random seed, call InitializeSeedsQueue first";
        G4Exception("G4MTRunManager::GetSeed","Run0035", FatalException,msg);
    }
    seeds[seedsnum++] = seed;
}

void G4MTRunManager::InitializeSeedsQueue( G4int ns )
{
    if ( seeds ) {
        delete[] seeds;
    }
    seeds = new long[ns];
    seedsnum = 0;
}
void G4MTRunManager::ConstructScoringWorlds()
{
    masterScM = G4ScoringManager::GetScoringManagerIfExist();
    //Call base class stuff...
    G4RunManager::ConstructScoringWorlds();
    
/*************************************************************************************

    //TODO: Understand with Makoto if mapping index<->Pworld is really needed here
    //Now need to remember the master scoring manager, to give access
    //to it to Workers
    masterWorlds.clear();
    if ( masterScM ) {
        G4int counter = 0;
        G4int nPar = masterScM->GetNumberOfMesh();
        for ( G4int iw = 0 ; iw<nPar ; ++iw)
        {
            G4VPhysicalVolume* pWorld = G4TransportationManager::GetTransportationManager()
                ->IsWorldExisting(masterScM->GetWorldName(iw));
            if (!pWorld)
            {
                addWorld(counter++,pWorld);
            }
        }
    }

****************************************************************************************/

    masterWorlds.clear();
    size_t nWorlds = G4TransportationManager::GetTransportationManager()->GetNoWorlds();
    std::vector<G4VPhysicalVolume*>::iterator itrW
                   = G4TransportationManager::GetTransportationManager()->GetWorldsIterator();
    for(size_t iWorld=0;iWorld<nWorlds;iWorld++)
    { addWorld(iWorld,*itrW); itrW++; }

}

void G4MTRunManager::SetUserInitialization(G4UserWorkerInitialization* userInit)
{
  userWorkerInitialization = userInit;
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
    userRunAction = userAction; 
    userRunAction->SetMaster();
}

void G4MTRunManager::SetUserAction(G4VUserPrimaryGeneratorAction* /*userAction*/)
{
  G4Exception("G4MTRunManager::SetUserAction()", "Run3011", FatalException,
    "For multi-threaded version, define G4VUserPrimaryGeneratorAction in G4VUserActionInitialization.");
}

void G4MTRunManager::SetUserAction(G4UserEventAction* /*userAction*/)
{
  G4Exception("G4MTRunManager::SetUserAction()", "Run3011", FatalException,
    "For multi-threaded version, define G4UserEventAction in G4VUserActionInitialization.");
}

void G4MTRunManager::SetUserAction(G4UserStackingAction* /*userAction*/)
{
  G4Exception("G4MTRunManager::SetUserAction()", "Run3011", FatalException,
    "For multi-threaded version, define G4UserStackingAction in G4VUserActionInitialization.");
}

void G4MTRunManager::SetUserAction(G4UserTrackingAction* /*userAction*/)
{
  G4Exception("G4MTRunManager::SetUserAction()", "Run3011", FatalException,
    "For multi-threaded version, define G4UserTrackingAction in G4VUserActionInitialization.");
}

void G4MTRunManager::SetUserAction(G4UserSteppingAction* /*userAction*/)
{
  G4Exception("G4MTRunManager::SetUserAction()", "Run3011", FatalException,
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


void G4MTRunManager::TerminateWorkers()
{
    NewActionRequest( ENDWORKER );
    //G4cout<<"ENDWORKER sent"<<G4endl;//AAAADEBUG
    //Now join threads.
    //TODO: I don't like this since this should go in the thread-model part of the code
    //something like: worker->join or something like this...
#ifdef G4MULTITHREADED //protect here to prevent warning in compilation
    
    //Andrea : We have observed deadlocks here. Under-investigation
    //For the time being we simply comment out everything and we do not
    // join threads: we arrive here at exit of application
    
    //while ( ! threads.empty() )
    //{
        //G4Thread* t = * ( threads.begin() );
        //threads.pop_front();
        
        //G4cout<<"In G4MTRunManager::TeraminteWorkers for "<<t<<G4endl;///AAADEBUG
        //G4THREADJOIN(*t);
        //G4cout<<t<<" Join done"<<G4endl;///AAADEBUG
        //delete t;
    //}
#endif
    threads.clear();

    G4AutoLock l(&workesRMMutex);
    workersRM.clear();
}

void G4MTRunManager::AddWorkerRunManager(G4WorkerRunManager* wrm)
{
    //This operation should be mutexed
    G4AutoLock l(&workesRMMutex);
    workersRM.push_back(wrm);
}


#include "G4IonTable.hh"
#include "G4ParticleTable.hh"
//#include "G4InuclElementaryParticle.hh"
void G4MTRunManager::InitializePhysics()
{
    G4RunManager::InitializePhysics();
    //G4ParticleTable::GetParticleTable()->GetIonTable()->CreateAllIon();
    G4ParticleTable::GetParticleTable()->GetIonTable()->CreateAllIsomer();
    //BERTINI, this is needed to create pseudo-particles, to be removed
    //G4InuclElementaryParticle::Initialize();
}

//Barriers mechanism
// We want to implement two barriers.
// We define a barrier has a point in which threads synchronize.
// When workers threads reach a barrier they wait for the master thread a
// signal that they can continue. The master thread broadcast this signal
// only when all worker threads have reached this point.
// Currently only two points require this sync:
// Just before and just after the for-loop controlling the thread event-loop.
// TODO: If this mechanism is needed in other parts of the code we can provide
// the barrier mechanism as a utility class/functions to the kernel.
// Note: we need a special treatment for WIN32
#ifdef WIN32
#include <windows.h> //For CRITICAL_SECTION objects
#endif
namespace {
    //Acoid compilation warning if squenetial for unused variables
#ifdef G4MULTITHREADED
    //Conditions
    // Condition to signal green light for start of event loop
    G4Condition beginEventLoopCondition = G4CONDITION_INITIALIZER;
    //pthread_cond_t beginEventLoopCondition = PTHREAD_COND_INITIALIZER;
    // Condition to signal green light to finish event loop (actuallyt exit function
    //performing event loop)
    G4Condition endEventLoopCondition = G4CONDITION_INITIALIZER;
    //pthread_cond_t endEventLoopCondition = PTHREAD_COND_INITIALIZER;
    // Condition to signal the num of workers ready for event loop has changed
    G4Condition numWorkersBeginEventLoopChangeCondition = G4CONDITION_INITIALIZER;
    //pthread_cond_t numWorkersBeginEventLoopChangeCondition = PTHREAD_COND_INITIALIZER;
    // Condition to signal the num of workers that terminated event loop has changed
    G4Condition numWorkersEndEventLoopChangedCondition = G4CONDITION_INITIALIZER;
    //pthread_cond_t numWorkersEndEventLoopChangedCondition = PTHREAD_COND_INITIALIZER;
    //Counter/mutex for workers ready to begin event loop
    //
    // This condition is to handle more than one run w/o killing threads
    G4Condition requestChangeActionForThread = G4CONDITION_INITIALIZER;
#endif
    G4Mutex numberOfReadyWorkersMutex = G4MUTEX_INITIALIZER;
    //G4Mutex numberOfReadyWorkersMutex = G4MUTEX_INITIALIZER;
    G4int numberOfReadyWorkers = 0;
    //Counter/mutex for workers with end of event loop
    G4Mutex numberOfEndOfEventLoopWorkersMutex = G4MUTEX_INITIALIZER;
    //G4Mutex numberOfEndOfEventLoopWorkersMutex = G4MUTEX_INITIALIZER;
    G4int numberOfEndOfEventLoopWorkers = 0;
    //
    //More than one 
    G4Mutex changeActionForWorkerMutex = G4MUTEX_INITIALIZER;
#ifdef WIN32
    CRITICAL_SECTION cs1;
    CRITICAL_SECTION cs2;
    CRITICAL_SECTION cs3;
    //Note we need to use two separate counters because
    //we can get a situation in which a thread is much faster then the others
    //(for example if last thread has less events to process.
    //We have the extreme case of some medical applications (moving setups)
    //in which the number of events of a run is ~ number of threads
	void InitializeWindowsConditions()
	{
	#ifdef G4MULTITHREADED
		InitializeConditionVariable( &beginEventLoopCondition );
		InitializeConditionVariable( &endEventLoopCondition );
		InitializeConditionVariable( &numWorkersBeginEventLoopChangeCondition );
		InitializeConditionVariable( &numWorkersEndEventLoopChangedCondition );
        InitializeConditionVariable( &requestChangeStatusForThread );
	#endif
		InitializeCriticalSection( &cs1 );
		InitializeCriticalSection( &cs2 );
        InitializeCriticalSection( &cs3 );
	}
#endif
}

void G4MTRunManager::WaitForReadyWorkers()
{
    while (true) //begin barrier
    {
#ifndef WIN32
        G4AutoLock lockLoop(&numberOfReadyWorkersMutex);
#else
        EnterCriticalSection( &cs1 );
#endif
        
        //Check number of workers ready to begin
        if (numberOfReadyWorkers == nworkers )
        {
            //Ok, interrupt the loop
            break;
        }
        //Wait for the number of workers to be changed
        G4CONDITIONWAIT(&numWorkersBeginEventLoopChangeCondition,
                        &numberOfReadyWorkersMutex);
#ifdef WIN32
        LeaveCriticalSection( &cs1 );
#endif
    }
    //Now number of workers is as expected.
    //Prepare to wait for workers to end eventloop
    //Reset number of workers in "EndOfEventLoop"
    G4AutoLock l(&numberOfEndOfEventLoopWorkersMutex);
    numberOfEndOfEventLoopWorkers = 0;
    
    //signal workers they can start the event-loop
    G4CONDTIONBROADCAST(&beginEventLoopCondition);
}

void G4MTRunManager::ThisWorkerReady()
{
    //Increament number of active worker by 1
#ifndef WIN32
    G4AutoLock lockLoop(&numberOfReadyWorkersMutex);
#else
    EnterCriticalSection( &cs1 );
#endif
    ++numberOfReadyWorkers;
    //Signal the number of workers has changed
    G4CONDTIONBROADCAST(&numWorkersBeginEventLoopChangeCondition);
    //Wait for condition to start eventloop
    G4CONDITIONWAIT(&beginEventLoopCondition,&numberOfReadyWorkersMutex);
#ifdef WIN32
    LeaveCriticalSection( &cs1 );
#endif
}


void G4MTRunManager::WaitForEndEventLoopWorkers()
{
    while (true)
    {
#ifndef WIN32
        G4AutoLock l(&numberOfEndOfEventLoopWorkersMutex);
#else
        EnterCriticalSection( &cs2 );
#endif
        if ( numberOfEndOfEventLoopWorkers == nworkers )
        {
            break;
        }
        G4CONDITIONWAIT(&numWorkersEndEventLoopChangedCondition,
                        &numberOfEndOfEventLoopWorkersMutex);
#ifdef WIN32
        LeaveCriticalSection( &cs2 );
#endif
    }
    //Now number of workers that reached end of event loop is as expected
    //Reset number of workers in ready for work state so a new run can start
    G4AutoLock l(&numberOfEndOfEventLoopWorkersMutex);
    numberOfReadyWorkers = 0;
    //Signal workers they can end event-loop
    G4CONDTIONBROADCAST(&endEventLoopCondition);
}

void G4MTRunManager::ThisWorkerEndEventLoop()
{
    //Increament number of workers in end of evnet loop by 1
#ifndef WIN32
    G4AutoLock l(&numberOfEndOfEventLoopWorkersMutex);
#else
    EnterCriticalSection( &cs2 );
#endif
    ++numberOfEndOfEventLoopWorkers;
    //Signale this number has changed
    G4CONDTIONBROADCAST(&numWorkersEndEventLoopChangedCondition);
    //Wait for condition to exit eventloop
    G4CONDITIONWAIT(&endEventLoopCondition,&numberOfEndOfEventLoopWorkersMutex);
#ifdef WIN32
    LeaveCriticalSection( &cs2 );
#endif
}

void G4MTRunManager::NewActionRequest(G4MTRunManager::WorkerActionRequest newRequest)
{
    //Master thread requests workers to do something new
    //First acquire lock on shared resource
#ifndef WIN32
    G4AutoLock l(&changeActionForWorkerMutex);
#else
    EnterCriticalSection( &cs3 );
#endif
    nextActionRequest = newRequest;

    //Broadcast signal that a new request is available
    G4CONDTIONBROADCAST(&requestChangeActionForThread);
#ifdef WIN32
    LeaveCriticalSection( &cs3 );
#endif
}

namespace  {
    static G4ThreadLocal G4bool processrequest = true;
}

G4MTRunManager::WorkerActionRequest G4MTRunManager::ThisWorkerWaitForNextAction()
{
    WorkerActionRequest result;
    while ( true )
    {
#ifndef WIN32
        G4AutoLock l(&changeActionForWorkerMutex);
#else
        EnterCriticalSection( &cs3 );
#endif
        //The very first time I arrive here,
        //I need to execute the nextRequest whatever it is.
        //The next time I come here I will execute this
        //block only if there is something to do
        if ( processrequest )
        {
            result = nextActionRequest;
            processrequest = false;
            break;
        }
        //Wait for condition from Master thread that a new action is requested.
        //This will block here until the condition is met
        G4CONDITIONWAIT(&requestChangeActionForThread,&changeActionForWorkerMutex);
        //If now I am here it means that the nextaction has changed, so I reset the
        //processrequest flag (that is TLS), so in the next loop I will execute it
        processrequest = true;
#ifdef WIN32
        LeaveCriticalSection( &cs3 );
#endif
    }
    return result;
}

