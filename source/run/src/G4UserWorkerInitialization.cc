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
#include "G4UserWorkerInitialization.hh"
#include "G4VUserActionInitialization.hh"
#include "G4UImanager.hh"
#include "G4VUserPhysicsList.hh"
#include "G4AutoLock.hh"
#include "G4ofstreamDestination.hh"
#include "G4coutIdDestination.hh"
#include <sstream>

G4ThreadLocal G4WorkerThread* G4UserWorkerInitialization::wThreadContext = 0;

#ifdef G4MULTITHREADED
G4Thread* G4UserWorkerInitialization::CreateAndStartWorker(G4WorkerThread* wTC)
{
    //Note: this method is called by G4MTRunManager, here we are still sequential
    //Create a new thread/worker structure
    G4Thread* worker = new G4Thread;
    G4THREADCREATE(worker,&G4UserWorkerInitialization::StartThread , wTC );
    return worker;
}
#else
G4Thread* G4UserWorkerInitialization::CreateAndStartWorker(G4WorkerThread*)
{
    return new G4Thread;
}
#endif



void* G4UserWorkerInitialization::StartThread( void* context )
{
    
    //!!!!!!!!!!!!!!!!!!!!!!!!!!
    //!!!!!! IMPORTANT !!!!!!!!!
    //!!!!!!!!!!!!!!!!!!!!!!!!!!
    //This is not sequential anymore and G4UserWorkerInitialization is
    // a shared user initialization class
    // This mean I cannot use data memebers from G4UserWorkerInitialization
    // unless they are invariant ("read-only") and can be safely shared. All the rest that is not
    // invariant should be incapsualted into the context (or, as for wThreadContext be G4ThreadLocal)
    //!!!!!!!!!!!!!!!!!!!!!!!!!!

    //Note: context is of type G4WorkerThread. It is passed by CreateAndStartWorker,
    //that was originally created by G4MTRunManager.

    //wThreadContext is a G4ThreadLocal variable, so this will work!
    wThreadContext = (G4WorkerThread*)context;

    //================
    //Step-0:
    //================
    //Initliazie per-thread stream-output
    //The following line is needed before we actually do IO initialization
    //becasue the constructor of UI manager resets the io destination.
    G4UImanager::GetUIpointer();
    //This creaetes the per-thread instances of cout and cerr buffe3rs
    G4iosInitialization();
    //This creates the output on file
    G4coutDestination *threadcout = 0;
    G4coutDestination *threadcerr = 0;
    if (! wThreadContext->GetOutputFileName().empty() )
    {
        threadcout = new G4CoutToFile;
        std::stringstream fn;
        fn << "G4Worker" <<wThreadContext->GetThreadId()<<"_"<< wThreadContext->GetOutputFileName();
        static_cast<G4CoutToFile*>(threadcout)->SetFileName(fn.str() , wThreadContext->GetOutputFileAppendFlag() );
    }
    else
    {
        //G4coutIdDestination gets both cout/cerr
        threadcout = new G4coutIdDestination( wThreadContext->GetThreadId() );
        static_cast<G4coutIdDestination*>(threadcout)->EnableBuffering( wThreadContext->GetOutputUseBuffer() );
    }
    if (! wThreadContext->GetOutputErrFileName().empty() )
    {
        //If this is created, it will overwrite what was set before for cout
        threadcerr = new G4CerrToFile;
        std::stringstream fn;
        fn << "G4Worker" <<wThreadContext->GetThreadId()<<"_"<< wThreadContext->GetOutputErrFileName();
        static_cast<G4CerrToFile*>(threadcerr)->SetFileName(fn.str() , wThreadContext->GetOutputErrFileAppendFlag() );
        //Force open of file: so the file always exists even if no message is sent there
        static_cast<G4CerrToFile*>(threadcerr)->Open();
    }
    
    //================
    //Step-1:
    //================
    //RNG Engine need to be initialized "cloning" the master one.
    //Since the
    G4MTRunManager* masterRM = G4MTRunManager::GetMasterRunManager(); 
    const CLHEP::HepRandomEngine* masterEngine = masterRM->getMasterRandomEngine();
    masterRM->GetUserWorkerInitialization()->SetupRNGEngine(masterEngine);
    masterRM->GetUserWorkerInitialization()->WorkerInitialize();
    
    //Now initialize worker part of shared objects (geometry/physics)
    wThreadContext->BuildGeometryAndPhysicsVector();
    
    //================
    //Step-2:
    //================
    //Create a G4WorkerRunManager
    G4WorkerRunManager* wrm = new G4WorkerRunManager;

    wrm->SetWorkerThread(wThreadContext);
    
    //================
    //Step-3:
    //================
    // Set the detector and physics list to the worker thread. Share with master
    const G4VUserDetectorConstruction* detector = masterRM->GetUserDetectorConstruction();
    wrm->G4RunManager::SetUserInitialization(const_cast<G4VUserDetectorConstruction*>(detector));
    const G4VUserPhysicsList* physicslist = masterRM->GetUserPhysicsList();
    wrm->SetUserInitialization(const_cast<G4VUserPhysicsList*>(physicslist));
    
    //================
    //Step-4:
    //================
    //Call user method to define user actions and all other stuff
    //Note: Even if thiscontext is per-thread object, the UserWorkerInitialization
    // object is shared
    if(masterRM->GetUserActionInitialization())
    { masterRM->GetUserActionInitialization()->Build(); }
    masterRM->GetUserWorkerInitialization()->WorkerStart();
    //Now initialize run manager
    wrm->Initialize();
    
    //================
    //Step-5:
    //================
    //Now enter a loop to execute requests from Master Thread: get next action
    G4MTRunManager::WorkerActionRequest nextAction = masterRM->ThisWorkerWaitForNextAction();
    while ( nextAction != G4MTRunManager::ENDWORKER )
    {
        if ( nextAction == G4MTRunManager::NEXTITERATION )
        {
            //Execute all stacked commands
            std::vector<G4String> cmds = masterRM->GetCommandStack();
            G4UImanager* uimgr = G4UImanager::GetUIpointer(); //TLS instance
            
            for ( std::vector<G4String>::const_iterator it = cmds.begin()
                 ; it != cmds.end() ; ++it )
            {
                std::cout<<"AAA"<<*it<<std::endl;
                uimgr->ApplyCommand(*it);
            }
            //TODO: when /run/beamOn is not passed, do it here!
            // wrm->BeamOn( wThreadContext->GetNumberEvents() , .... );
            //wrm->BeamOn(1); //AND <-Simulate a second call to /run/beamOn
            //TODO: move this stuff away since it will be allowed to have always the same threads for more runs
            //in interactive mode
            
        }
        else
        {
            G4ExceptionDescription d;
            d<<"Cannot continue, this worker has been requested an unknwon action: "<<nextAction<<" expecting: ENDWORKER(=)"
            <<G4MTRunManager::ENDWORKER<<") or NEXTITERATION(="<<G4MTRunManager::NEXTITERATION<<")";
            G4Exception("G4WorkerInitialization::WorkerStart","Run0035",FatalException,d);
        }
        //Now wait for master thread to signal new action to be performed
        nextAction = masterRM->ThisWorkerWaitForNextAction();
    } //No more actions to perform
    
    //================
    //Step-6:
    //================
    //Ok, now there is nothing else to do , called user defined stuff
    masterRM->GetUserWorkerInitialization()->WorkerStop();
    wThreadContext->DestroyGeometryAndPhysicsVector();

    // Delete worker run manager
    //G4cout<<"Thread ID:"<<wThreadContext->GetThreadId()<<" WorkerRunManager Pointer: "<<wrm<<" PID:"<<wThreadContext->pid<<G4endl;///AAADEBUG
    delete wrm;

    if ( threadcout ) delete threadcout;
    if ( threadcerr ) delete threadcerr;
    G4iosFinalization();
    return static_cast<void*>(0);
}

G4UserWorkerInitialization::G4UserWorkerInitialization()
{;}

G4UserWorkerInitialization::~G4UserWorkerInitialization()
{;}

void G4UserWorkerInitialization::WorkerInitialize() const
{;}

void G4UserWorkerInitialization::WorkerStart() const
{;}

void G4UserWorkerInitialization::WorkerRunStart() const
{;}

void G4UserWorkerInitialization::WorkerRunEnd() const
{;}

void G4UserWorkerInitialization::WorkerStop() const
{;}

void G4UserWorkerInitialization::SetUserAction(G4VUserPrimaryGeneratorAction* action) const
{ G4RunManager::GetRunManager()->SetUserAction(action); } 

void G4UserWorkerInitialization::SetUserAction(G4UserRunAction* action) const
{ G4RunManager::GetRunManager()->SetUserAction(action); } 

void G4UserWorkerInitialization::SetUserAction(G4UserEventAction* action) const
{ G4RunManager::GetRunManager()->SetUserAction(action); } 

void G4UserWorkerInitialization::SetUserAction(G4UserStackingAction* action) const
{ G4RunManager::GetRunManager()->SetUserAction(action); } 

void G4UserWorkerInitialization::SetUserAction(G4UserTrackingAction* action) const
{ G4RunManager::GetRunManager()->SetUserAction(action); } 

void G4UserWorkerInitialization::SetUserAction(G4UserSteppingAction* action) const
{ G4RunManager::GetRunManager()->SetUserAction(action); } 


void G4UserWorkerInitialization::SetupRNGEngine(const CLHEP::HepRandomEngine* aNewRNG) const
{
    //No default available, let's create the instance of random stuff
    //A Call to this just forces the creation to defaults
    G4Random::getTheEngine();
    //Poor man's solution to check which RNG Engine is used in master thread
    if ( dynamic_cast<const CLHEP::HepJamesRandom*>(aNewRNG) ) { G4Random::setTheEngine(new CLHEP::HepJamesRandom); return; }
    if ( dynamic_cast<const CLHEP::RanecuEngine*>(aNewRNG) ) { G4Random::setTheEngine(new CLHEP::RanecuEngine); return; }
    if ( dynamic_cast<const CLHEP::Ranlux64Engine*>(aNewRNG) ) { G4Random::setTheEngine(new CLHEP::Ranlux64Engine); return; }
    if ( dynamic_cast<const CLHEP::MTwistEngine*>(aNewRNG) ) { G4Random::setTheEngine(new CLHEP::MTwistEngine); return; }
    if ( dynamic_cast<const CLHEP::DualRand*>(aNewRNG) ) { G4Random::setTheEngine(new CLHEP::DualRand); return; }
    if ( dynamic_cast<const CLHEP::RanluxEngine*>(aNewRNG) ) { G4Random::setTheEngine(new CLHEP::RanluxEngine); return;}
    if ( dynamic_cast<const CLHEP::RanshiEngine*>(aNewRNG) ) { G4Random::setTheEngine(new CLHEP::RanshiEngine); return; }
}

