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
#include "G4VUserWorkerInitialization.hh"
#include "G4WorkerThread.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4AutoLock.hh"
#include "G4WorkerRunManager.hh"

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
}

G4MTRunManager* G4MTRunManager::GetMasterRunManager()
{
    return masterRM;
}

G4RunManagerKernel* G4MTRunManager::GetMasterRunManagerKernel()
{
    return masterRM->kernel;
}

G4MTRunManager::G4MTRunManager() : G4RunManager(false) ,
    nworkers(0),
    masterRNGEngine(0)
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
}

G4MTRunManager::~G4MTRunManager()
{
    //TODO: Currently does not work due to concurrent deletion of something that is shared:
    //G4ProcessTable::DeleteMessenger from ~G4RunManager
    //G4cout<<"Destroy MTRunManager"<<G4endl;//ANDREA
    //DestroyWorkers();
    delete[] seeds;
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
    return -DBL_MAX;
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
    InitializeSeeds(n_events);
    //IF seeds is empty, user did not implement InitializeSeeds, use default: 2 seeds per event number
    if ( NumberOfAvailableSeeds() == 0 )
    {
        InitializeSeedsQueue(n_events*2);
        for ( G4int ne = 0 ; ne < n_events*2 ; ++ne)
            AddOneSeed( (long) (100000000L * CLHEP::HepRandom::getTheGenerator()->flat()) );
    }
    
    //Now initialize workers. Check if user defined a WorkerInitialization
    if ( userWorkerInitialization == 0 )
    {
        //Throw exception
    }
    
    //Prepare UI commands for threads
    PrepareCommandsStack();
    //Now loop on requested number of workers
    //This will also start the workers
    for ( G4int nw = 0 ; nw<nworkers; ++nw) {
        //Create a new worker and remember it
        G4WorkerThread* context = new G4WorkerThread;
        context->SetNumberThreads(nworkers);
        context->SetNumberEvents(n_events);
        context->SetThreadId(nw);
        G4Thread* thread = userWorkerInitialization->CreateAndStartWorker(context);
        threads.push_back(thread);
    }
    //TODO: Should I setup a barrier here?
    //Probably yes: if threads are many many it can take some time before they are ready.
    //Now threads can be signaled with the commands to be executed
    //std::vector<G4String>* cmds = G4UImanager::GetUIpointer()->GetCommandStack();
    //for ( int i = 0 ; i < uiCmdsForWorkers.size() ; ++ i ) G4cout<<"AAA"<<uiCmdsForWorkers[i]<<G4endl;
    //G4String cmdWithError;
    //for ( G4WorkerRunManagerList::iterator it = workersRM.begin();
    //     it != workersRM.end(); ++it )
    //{
        //G4int err = (*it)->NewCommands(cmds, cmdWithError );
        //if ( err != 0 )
        //{
          //  G4ExceptionDescription msg;
            //msg<<"Worker cannot execute command:"<<cmdWithError<<" Return code: "<<err;
          //  G4Exception("G4MTRunManager::InitializeEventLoop","Run0035", FatalException,msg);
        //}
        
    //}
    //Now join threads.
    //I don't like this since this should go in the thread-model part of the code
    //something like: worker->join or something like this...
#ifdef G4MULTITHREADED //protect here to prevent warning in compilation
    for ( G4ThreadsList::iterator tit = threads.begin() ; tit != threads.end() ; ++tit ){
        G4Thread* t = *tit;
        G4THREADJOIN(*t);
    }
#endif
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

void G4MTRunManager::SetUserInitialization(G4VUserWorkerInitialization* userInit)
{
  userWorkerInitialization = userInit;
}

void G4MTRunManager::SetUserInitialization(G4VUserPhysicsList *userPL)
{
    G4RunManager::SetUserInitialization(userPL);
}

void G4MTRunManager::SetUserInitialization(G4VUserDetectorConstruction *userDC)
{
    G4RunManager::SetUserInitialization(userDC);
}

void G4MTRunManager::SetUserAction(G4UserRunAction* userAction)
{ userRunAction = userAction; }

void G4MTRunManager::SetUserAction(G4VUserPrimaryGeneratorAction* /*userAction*/)
{
  G4Exception("G4MTRunManager::SetUserAction()", "Run3011", FatalException,
    "For multi-threaded version, define G4VUserPrimaryGeneratorAction in G4VUserWorkerInitialization.");
}

void G4MTRunManager::SetUserAction(G4UserEventAction* /*userAction*/)
{
  G4Exception("G4MTRunManager::SetUserAction()", "Run3011", FatalException,
    "For multi-threaded version, define G4UserEventAction in G4VUserWorkerInitialization.");
}

void G4MTRunManager::SetUserAction(G4UserStackingAction* /*userAction*/)
{
  G4Exception("G4MTRunManager::SetUserAction()", "Run3011", FatalException,
    "For multi-threaded version, define G4UserStackingAction in G4VUserWorkerInitialization.");
}

void G4MTRunManager::SetUserAction(G4UserTrackingAction* /*userAction*/)
{
  G4Exception("G4MTRunManager::SetUserAction()", "Run3011", FatalException,
    "For multi-threaded version, define G4UserTrackingAction in G4VUserWorkerInitialization.");
}

void G4MTRunManager::SetUserAction(G4UserSteppingAction* /*userAction*/)
{
  G4Exception("G4MTRunManager::SetUserAction()", "Run3011", FatalException,
    "For multi-threaded version, define G4UserSteppingAction in G4VUserWorkerInitialization.");
}

void G4MTRunManager::MergeScores(const G4ScoringManager* localScoringManager)
{ if(masterScM) masterScM->Merge(localScoringManager); }

void G4MTRunManager::MergeRun(const G4Run* localRun)
{ if(currentRun) currentRun->Merge(localRun); }


namespace  {
//predicate implemented as a function
static bool DestroyWorkerRunManager( G4WorkerRunManager* rm ) {
    delete rm;
    return true;
}
}

void G4MTRunManager::DestroyWorkers()
{
    //This operation should be mutexed
    G4AutoLock l(&workesRMMutex);
    workersRM.remove_if(DestroyWorkerRunManager);
}

void G4MTRunManager::AddWorkerRunManager(G4WorkerRunManager* wrm)
{
    //This operation should be mutexed
    G4AutoLock l(&workesRMMutex);
    workersRM.push_back(wrm);
}
