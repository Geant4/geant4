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
#include <sstream>

//Will need this for TPMalloc
//#ifdef G4MULTITHREADED
//#define TPMALLOCDEFINESTUB
//#include "tpmalloc/tpmallocstub.h"
//#endif

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




#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVReplica.hh"
#include "G4ParticleDefinition.hh"
#include "G4Region.hh"
#include "G4Material.hh"
#include "G4PhysicsVector.hh"
#include "G4VDecayChannel.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4MaterialTable.hh"
#include "G4PolyconeSide.hh"
#include "G4PolyhedraSide.hh"
#include "G4PVParameterised.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4VModularPhysicsList.hh"




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
//#ifdef G4MULTITHREADED
//    turnontpmalloc();
//#endif
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
    G4int thisID = wThreadContext->GetThreadId();
    G4UImanager::GetUIpointer()->SetUpForAThread(thisID);
    //Set the global thread-ID for this thread,
    //using a global function
    G4Threading::G4SetThreadId(thisID);
    
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
    G4WorkerRunManager* wrm = masterRM->GetUserWorkerInitialization()->CreateWorkerRunManager();

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
    //Note: Even if this context is per-thread object, the UserWorkerInitialization
    // object is shared
    if(masterRM->GetUserActionInitialization())
    { masterRM->GetUserActionInitialization()->Build(); }
    masterRM->GetUserWorkerInitialization()->WorkerStart();
    //Now initialize worker run manager
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
            /**************** TEMPORARY start ***********************
                            A.Dotti 20Jun2013
             The following code is a temporary solution to allow multiple runs changing
             material between runs for 10.0-beta version
             We re-do here for every run except the first the re-initialization of 
             geometry for worker threads. To be put in correct place and re-viewed.
             */
            static G4ThreadLocal G4bool skip = true;
            if ( !skip ) {
                //I need to rememeber SD and Filed Associated with worker, so to re-use it
                // (note that all the stuff after this will reset SD and Field
                typedef std::map<G4LogicalVolume*,std::pair<G4VSensitiveDetector*,G4FieldManager*> > LV2SDFM;
                LV2SDFM lvmap;
                G4PhysicalVolumeStore* mphysVolStore = G4PhysicalVolumeStore::GetInstance();
                for(size_t ip=0; ip<mphysVolStore->size(); ip++)
                {
                    G4VPhysicalVolume* pv = (*mphysVolStore)[ip];
                    G4LogicalVolume *lv = pv->GetLogicalVolume();
                    G4VSensitiveDetector* sd = lv->GetSensitiveDetector();
                    G4FieldManager* fm = lv->GetFieldManager();
                    if ( sd || fm )
                        lvmap[lv] = std::make_pair(lv->GetSensitiveDetector(),lv->GetFieldManager());
                }
                
                //First delete stuff
                const_cast<G4LVManager&>(G4LogicalVolume::GetSubInstanceManager()).FreeSlave();
                const_cast<G4PVManager&>(G4VPhysicalVolume::GetSubInstanceManager()).FreeSlave();
                const_cast<G4PVRManager&>(G4PVReplica::GetSubInstanceManager()).FreeSlave();
                //const_cast<G4PDefManager&>(G4ParticleDefinition::GetSubInstanceManager()).FreeSlave();
                const_cast<G4RegionManager&>(G4Region::GetSubInstanceManager()).FreeSlave();
                //const_cast<G4PVecManager&>(G4PhysicsVector::GetSubInstanceManager()).FreeSlave();
                const_cast<G4PlSideManager&>(G4PolyconeSide::GetSubInstanceManager()).FreeSlave();
                const_cast<G4PhSideManager&>(G4PolyhedraSide::GetSubInstanceManager()).FreeSlave();
                //Now-re-create
                const_cast<G4LVManager&>(G4LogicalVolume::GetSubInstanceManager()).SlaveCopySubInstanceArray();
                const_cast<G4PVManager&>(G4VPhysicalVolume::GetSubInstanceManager()).SlaveCopySubInstanceArray();
                const_cast<G4PVRManager&>(G4PVReplica::GetSubInstanceManager()).SlaveCopySubInstanceArray();
                const_cast<G4RegionManager&>(G4Region::GetSubInstanceManager()).SlaveInitializeSubInstance();
                //G4Material::g4materialSubInstanceManager.SlaveCopySubInstanceArray(); //< Not anymore splitted class
                const_cast<G4PlSideManager&>(G4PolyconeSide::GetSubInstanceManager()).SlaveInitializeSubInstance();
                const_cast<G4PhSideManager&>(G4PolyhedraSide::GetSubInstanceManager()).SlaveInitializeSubInstance();
            
                G4PhysicalVolumeStore* physVolStore = G4PhysicalVolumeStore::GetInstance();
                for(size_t ip=0; ip<physVolStore->size(); ip++)
                {
                    G4VPhysicalVolume* physVol = (*physVolStore)[ip];
                    G4LogicalVolume *g4LogicalVolume = physVol->GetLogicalVolume();
                    //use shadow pointer
                    G4VSolid *g4VSolid = g4LogicalVolume->GetMasterSolid();
                    G4PVReplica *g4PVReplica = 0;
                    g4PVReplica =  dynamic_cast<G4PVReplica*>(physVol);
                    if (g4PVReplica)
                    {
                        //g4PVReplica->SlaveG4PVReplica(g4PVReplica);
                        g4PVReplica->InitialiseWorker(g4PVReplica);
                        G4PVParameterised *g4PVParameterised = 0;
                        g4PVParameterised =  dynamic_cast<G4PVParameterised*>(physVol);
                        if (g4PVParameterised)
                        {
                            //01.25.2009 Xin Dong: For a G4PVParameterised instance, assoicated a
                            //cloned solid for each worker thread. If all solids support this clone
                            //method, we do not need to dynamically cast to solids that support this
                            //clone method. Before all solids support this clone method, we do similar
                            //thing here to dynamically cast and then get the clone method.
                        
                            //Threads may clone some solids simultaneously. Those cloned solids will be
                            //Registered into a shared solid store (C++ container). Need a lock to
                            //guarantee thread safety
                            //G4AutoLock aLock(&solidclone);
                            G4VSolid *slaveg4VSolid = g4VSolid->Clone();
                            //aLock.unlock();
                            //g4LogicalVolume->SlaveG4LogicalVolume(g4LogicalVolume, slaveg4VSolid, 0);
                            g4LogicalVolume->InitialiseWorker(g4LogicalVolume,slaveg4VSolid,0);
                        }
                        else
                        {
                            //g4LogicalVolume->SlaveG4LogicalVolume(g4LogicalVolume, g4VSolid, 0);
                            g4LogicalVolume->InitialiseWorker(g4LogicalVolume,g4VSolid,0);
                        }
                    }
                    else
                    {
                        //g4LogicalVolume->SlaveG4LogicalVolume(g4LogicalVolume, g4VSolid, 0);
                        g4LogicalVolume->InitialiseWorker(g4LogicalVolume,g4VSolid,0);
                    }
                }
                // Now set back SD and FM pointers to logical volumes
                for ( LV2SDFM::const_iterator it = lvmap.begin() ; it != lvmap.end() ; ++it )
                {
                    G4LogicalVolume* lv      = it->first;
                    G4VSensitiveDetector* sd = (it->second).first;
                    G4FieldManager* fm       = (it->second).second;
                    lv->SetFieldManager(fm, false); //What should be the second parameter? We use always false for MT mode
                    lv->SetSensitiveDetector(sd);
                }
            } else {
                skip = false;
            }
            /**************** TEMPORARY end ***********************
             */
            
            //Execute all stacked commands
            std::vector<G4String> cmds = masterRM->GetCommandStack();
            G4UImanager* uimgr = G4UImanager::GetUIpointer(); //TLS instance
            
            for ( std::vector<G4String>::const_iterator it = cmds.begin()
                 ; it != cmds.end() ; ++it )
            {
/////////                std::cout<<"AAA"<<*it<<std::endl;
                uimgr->ApplyCommand(*it);
            }
            //TODO: when /run/beamOn is not passed, do it here!
            ///wrm->BeamOn( 000000wThreadContext->GetNumberEvents() , .... );
            //wrm->BeamOn(1); //AND <-Simulate a second call to /run/beamOn
            //TODO: move this stuff away since it will be allowed to have always the same threads for more runs
            //in interactive mode
            
        }
        else
        {
            G4ExceptionDescription d;
            d<<"Cannot continue, this worker has been requested an unknwon action: "
             <<nextAction<<" expecting: ENDWORKER(="
             <<G4MTRunManager::ENDWORKER<<") or NEXTITERATION(="<<G4MTRunManager::NEXTITERATION<<")";
            G4Exception("G4UserWorkerInitialization::StartThread","Run0035",FatalException,d);
        }
        //Now wait for master thread to signal new action to be performed
        nextAction = masterRM->ThisWorkerWaitForNextAction();
    } //No more actions to perform
    
    //================
    //Step-6:
    //================
    //Ok, now there is nothing else to do , called user defined stuff
    masterRM->GetUserWorkerInitialization()->WorkerStop();
    delete wrm;
    wThreadContext->DestroyGeometryAndPhysicsVector();

    //G4cout<<"Thread ID:"<<wThreadContext->GetThreadId()<<" WorkerRunManager Pointer: "<<wrm<<G4endl;///AAADEBUG


//#ifdef G4MULTITHREADED
//    turnofftpmalloc();
//#endif
    return static_cast<void*>(0);
}

//Avoid compilation warning in sequential
#ifdef G4MULTITHREADED
void G4UserWorkerInitialization::JoinWorker(G4Thread* aThread)
{
    G4THREADJOIN(*aThread);
}
#else
void G4UserWorkerInitialization::JoinWorker(G4Thread*)
{
}
#endif

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

namespace {
    G4Mutex rngCreateMutex = G4MUTEX_INITIALIZER;
}

void G4UserWorkerInitialization::SetupRNGEngine(const CLHEP::HepRandomEngine* aNewRNG) const
{
    //No default available, let's create the instance of random stuff
    //A Call to this just forces the creation to defaults
    G4Random::getTheEngine();
    //Poor man's solution to check which RNG Engine is used in master thread
    // Need to make these calls thread safe
    G4AutoLock l(&rngCreateMutex);
    if ( dynamic_cast<const CLHEP::HepJamesRandom*>(aNewRNG) ) { G4Random::setTheEngine(new CLHEP::HepJamesRandom); return; }
    if ( dynamic_cast<const CLHEP::RanecuEngine*>(aNewRNG) ) { G4Random::setTheEngine(new CLHEP::RanecuEngine); return; }
    if ( dynamic_cast<const CLHEP::Ranlux64Engine*>(aNewRNG) ) { G4Random::setTheEngine(new CLHEP::Ranlux64Engine); return; }
    if ( dynamic_cast<const CLHEP::MTwistEngine*>(aNewRNG) ) { G4Random::setTheEngine(new CLHEP::MTwistEngine); return; }
    if ( dynamic_cast<const CLHEP::DualRand*>(aNewRNG) ) { G4Random::setTheEngine(new CLHEP::DualRand); return; }
    if ( dynamic_cast<const CLHEP::RanluxEngine*>(aNewRNG) ) { G4Random::setTheEngine(new CLHEP::RanluxEngine); return;}
    if ( dynamic_cast<const CLHEP::RanshiEngine*>(aNewRNG) ) { G4Random::setTheEngine(new CLHEP::RanshiEngine); return; }
}

G4WorkerRunManager* G4UserWorkerInitialization::CreateWorkerRunManager() const
{
    return new G4WorkerRunManager();
}

