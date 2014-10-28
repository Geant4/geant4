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
#include "G4WorkerThread.hh"
#include "G4WorkerRunManager.hh"
#include "G4MTRunManager.hh"

void G4WorkerThread::SetThreadId(G4int tid)
{
    threadId = tid;
}

G4int G4WorkerThread::GetThreadId() const
{
    return threadId;
}

void G4WorkerThread::SetNumberThreads(G4int nw)
{
    numThreads = nw;
}

G4int G4WorkerThread::GetNumberThreads() const
{
    return numThreads;
}

#include "G4GeometryWorkspace.hh"
#include "G4GeometryWorkspacePool.hh"
#include "G4SolidsWorkspace.hh"
#include "G4SolidsWorkspacePool.hh"
#include "G4ParticlesWorkspace.hh"
#include "G4PhysicsListWorkspace.hh"

void G4WorkerThread::BuildGeometryAndPhysicsVector()
{
    // Initialise all split classes in the geometry with copy of data from master thread
    G4GeometryWorkspacePool::GetInstance()->CreateAndUseWorkspace();
    G4SolidsWorkspacePool::GetInstance()->CreateAndUseWorkspace();

    G4ParticlesWorkspace::GetPool()->CreateAndUseWorkspace();
    G4PhysicsListWorkspace::GetPool()->CreateAndUseWorkspace();
//    const_cast<G4VUPLManager&>(G4VUserPhysicsList::GetSubInstanceManager()).NewSubInstances();
//    const_cast<G4VPCManager&>(G4VPhysicsConstructor::GetSubInstanceManager()).NewSubInstances();
//    const_cast<G4VMPLManager&>(G4VModularPhysicsList::GetSubInstanceManager()).WorkerCopySubInstanceArray();
}

void G4WorkerThread::DestroyGeometryAndPhysicsVector()
{
#if 0
    // Manage Geometry Workspace explicitly
    fGeometryWrk->ReleaseAndDestroyWorkspace();
    delete fGeometryWrk;
    fGeometryWrk=0;
#else
    // Alternative:
    // Initialise all split classes in the geometry with copy of data from master thread
    // G4GeometryWorkspacePool::GetInstance()->ReleaseAndDestroyMyWorkspace();
    G4GeometryWorkspacePool::GetInstance()->GetWorkspace()->DestroyWorkspace();
    G4SolidsWorkspacePool::GetInstance()->GetWorkspace()->DestroyWorkspace();
    G4ParticlesWorkspace::GetPool()->GetWorkspace()->DestroyWorkspace();
    
    G4PhysicsListWorkspace::GetPool()->GetWorkspace()->DestroyWorkspace();
#endif

//    const_cast<G4VUPLManager&>(G4VUserPhysicsList::GetSubInstanceManager()).FreeWorker();
//    const_cast<G4VPCManager&>(G4VPhysicsConstructor::GetSubInstanceManager()).FreeWorker();
//    const_cast<G4VMPLManager&>(G4VModularPhysicsList::GetSubInstanceManager()).FreeWorker();
}

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Region.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"


void G4WorkerThread::UpdateGeometryAndPhysicsVectorFromMaster()
{
    //=================================================
    //Step-0: keep sensitive detector and field manager
    //=================================================
    // First remember SD and Filed Associated with worker
    //  in order to re-use it
    // (note that all the stuff after this will reset SD and Field
    typedef std::map<G4LogicalVolume*,std::pair<G4VSensitiveDetector*,G4FieldManager*> > LV2SDFM;
    LV2SDFM lvmap;
//    typedef std::map<G4LogicalVolume*,std::pair<G4Region*,G4bool> > LV2Region;
//    LV2Region lv2rmap;
    typedef std::map<G4Region*,std::pair<G4FastSimulationManager*,G4UserSteppingAction*> > R2FSM;
    R2FSM rgnmap;
    G4LogicalVolumeStore* mLogVolStore = G4LogicalVolumeStore::GetInstance();
    for(size_t ip=0; ip<mLogVolStore->size(); ip++)
    {
        G4LogicalVolume *lv = (*mLogVolStore)[ip];
        //The following needs an explanation.
        //Consider the case in which the user adds one LogVolume between the runs. The problem is that the thread-local part
        //(split class) of the G4LogicalVolume object is not initialized for workers because the initialization is done once when the
        //thread starts (see G4MTRunManagerKernel::StartThread Step-2 that calls G4WorkerThread::BuildGeometryAndPhysicsVector in this class)
        //The problem is that pointers of SD and FM for these newly added LV
        //may be invalid pointers (because never initialized, we have seen this behavior in our testing). If now we remember
        //them and re-use them in Step-4 below we set invalid pointers to LV for this thread.
        //Thus we need a way to know if for a given LV we need to remember or not the SD and FM pointers.
        //To solve this problem: We assume that the ConstructSDandField is called also by Master thread
        //thus for newly added LV the shadow pointers of SD and Fields are correct.
        // (LIMITATION: this assumption may be too stringent, a user to save memory could instantiate SD only
        // for workers, but we require this not to happen!).
        // Thus is a SD and FieldMgr are needed for this particular LV, and shadow are !=0 it means that user
        // wants an SD and FM to be associated with LV, we get the values and we remember them.
        G4VSensitiveDetector* sd = 0;
        G4FieldManager* fmgr = 0;
        if ( lv->GetMasterSensitiveDetector() != 0 ) sd = lv->GetSensitiveDetector();
        if ( lv->GetMasterFieldManager() != 0 ) fmgr = lv->GetFieldManager();
        if ( sd || fmgr ) lvmap[lv] = std::make_pair(sd,fmgr);
//        G4Region* rgn = lv->GetRegion();
//        G4bool isRoot = lv->IsRootRegion();
//        if ( rgn || isRoot ) lv2rmap[lv] = std::make_pair(rgn,isRoot);
    }
    G4RegionStore* mRegStore = G4RegionStore::GetInstance();
    for(size_t ir=0; ir<mRegStore->size(); ir++)
    {
        G4Region* reg = (*mRegStore)[ir];
        G4FastSimulationManager* fsm = reg->GetFastSimulationManager();
        G4UserSteppingAction* usa = reg->GetRegionalSteppingAction();
        if ( reg || usa ) rgnmap[reg] = std::make_pair(fsm,usa);
    }

    //===========================
    //Step-1: Clean the workspace
    //===========================
    G4GeometryWorkspace* geomWorkspace= G4GeometryWorkspacePool::GetInstance()->GetWorkspace();
    geomWorkspace->DestroyWorkspace();
    G4SolidsWorkspace* solidWorkspace = G4SolidsWorkspacePool::GetInstance()->GetWorkspace();
    solidWorkspace->DestroyWorkspace();
    
    //===========================
    //Step-2: Re-create and initialize workspace
    //===========================
    geomWorkspace->InitialiseWorkspace();
    solidWorkspace->InitialiseWorkspace();

    // Alternative
    //
    // To wipe, release and get a new one ...
    // G4GeometryWorkspacePool *fWorkspaceMgr=  G4GeometryWorkspaceOutlet::GetInstance();
    // fWorkspaceMgr->ReleaseAndDestroyMyWorkspace();
    // Now re-create
    // fWorkspaceMgr->CreateAndUseWorkspace();
    
    
    //===================================================
    //Step-4: Restore sensitive detector and field manaer
    //===================================================
    for ( LV2SDFM::const_iterator it = lvmap.begin() ; it != lvmap.end() ; ++it )
    {
        G4LogicalVolume* lv      = it->first;
        G4VSensitiveDetector* sd = (it->second).first;
        G4FieldManager* fmgr       = (it->second).second;
        if(fmgr) lv->SetFieldManager(fmgr, false); //What should be the second parameter? We use always false for MT mode
        if(sd) lv->SetSensitiveDetector(sd);
    }
    
//    for ( LV2Region::const_iterator it2 = lv2rmap.begin() ; it2 != lv2rmap.end() ; it2++ )
//    {
//        G4LogicalVolume* lv2 = it2->first;
//        G4Region* rgn = (it2->second).first;
//        if(rgn) lv2->SetRegion(rgn);
//        G4bool isRoot = (it2->second).second;
//        lv2->SetRegionRootFlag(isRoot);
//    }

    
    for ( R2FSM::const_iterator it3 = rgnmap.begin() ; it3 != rgnmap.end() ; it3++ )
    {
        G4Region* reg = it3->first;
        G4FastSimulationManager* fsm = (it3->second).first;
        if(fsm) reg->SetFastSimulationManager(fsm);
        G4UserSteppingAction* usa = (it3->second).second;
        if(usa) reg->SetRegionalSteppingAction(usa);
    }
}

