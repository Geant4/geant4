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
//#include "G4VUserWorkerInitialization.hh"

//void G4WorkerThread::SetUserWorkerInitialization(G4VUserWorkerInitialization *userWorkerInit)
//{
//    uWorkerInit = userWorkerInit;
//}

//G4VUserWorkerInitialization* G4WorkerThread::GetUserWorkerInitialization() const
//{
//    return uWorkerInit;
//}

void G4WorkerThread::SetThreadId(G4int tid)
{
    threadId = tid;
}

G4int G4WorkerThread::GetThreadId() const
{
    return threadId;
}

void G4WorkerThread::SetNumberEvents(G4int nevt)
{
    totalNumEvents = nevt;
}

G4int G4WorkerThread::GetNumberEvents() const
{
    return totalNumEvents;
}

void G4WorkerThread::SetNumberThreads(G4int nw)
{
    numThreads = nw;
}

G4int G4WorkerThread::GetNumberThreads() const
{
    return numThreads;
}

//TODO: refactorize? Parts of these delegated to correct classes?
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

namespace  {
    G4Mutex solidclone = G4MUTEX_INITIALIZER;
}

void G4WorkerThread::BuildGeometryAndPhysicsVector()
{        
    //Geometry related, split classes mechanism: instantiate sub-instance for this thread
    const_cast<G4LVManager&>(G4LogicalVolume::GetSubInstanceManager()).SlaveCopySubInstanceArray();
    const_cast<G4PVManager&>(G4VPhysicalVolume::GetSubInstanceManager()).SlaveCopySubInstanceArray();
    const_cast<G4PVRManager&>(G4PVReplica::GetSubInstanceManager()).SlaveCopySubInstanceArray();
    const_cast<G4RegionManager&>(G4Region::GetSubInstanceManager()).SlaveInitializeSubInstance();
    //G4Material::g4materialSubInstanceManager.SlaveCopySubInstanceArray(); //< Not anymore splitted class
    const_cast<G4PlSideManager&>(G4PolyconeSide::GetSubInstanceManager()).SlaveInitializeSubInstance();
    const_cast<G4PhSideManager&>(G4PolyhedraSide::GetSubInstanceManager()).SlaveInitializeSubInstance();
    //Physics related
    //const_cast<G4PVecManager&>(G4PhysicsVector::GetSubInstanceManager()).SlaveInitializeSubInstance();
    const_cast<G4DecayChannelManager&>(G4VDecayChannel::GetSubInstanceManager()).NewSubInstances();
    const_cast<G4PDefManager&>(G4ParticleDefinition::GetSubInstanceManager()).NewSubInstances();
    const_cast<G4VUPLManager&>(G4VUserPhysicsList::GetSubInstanceManager()).NewSubInstances();
    const_cast<G4VPCManager&>(G4VPhysicsConstructor::GetSubInstanceManager()).NewSubInstances();
    const_cast<G4VMPLManager&>(G4VModularPhysicsList::GetSubInstanceManager()).SlaveCopySubInstanceArray();
    
    G4PhysicalVolumeStore* physVolStore = G4PhysicalVolumeStore::GetInstance();
    for (size_t ip=0; ip<physVolStore->size(); ip++)
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
                G4AutoLock aLock(&solidclone);
                G4VSolid *slaveg4VSolid = g4VSolid->Clone();
                aLock.unlock();
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
    
    //const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
    
    //size_t nmat = theMaterialTable->size();
    //size_t i;
    //for(i=0; i<nmat; i++) {
    //    ((*theMaterialTable)[i])->SlaveG4Material();
    //}
}

void G4WorkerThread::DestroyGeometryAndPhysicsVector()
{
    G4PhysicalVolumeStore* physVolStore = G4PhysicalVolumeStore::GetInstance();
    for (size_t ip=0; ip<physVolStore->size(); ip++)
    {
        G4VPhysicalVolume* physVol = (*physVolStore)[ip];
        G4LogicalVolume *g4LogicalVolume = physVol->GetLogicalVolume();
        //    G4VSolid *g4VSolid = g4LogicalVolume->fSolid;
        G4PVReplica *g4PVReplica = 0;
        g4PVReplica =  dynamic_cast<G4PVReplica*>(physVol);
        if (g4PVReplica)
        {
            //g4PVReplica->DestroySlaveG4PVReplica(g4PVReplica);
            g4PVReplica->TerminateWorker(g4PVReplica);
            G4PVParameterised *g4PVParameterised = 0;
            g4PVParameterised =  dynamic_cast<G4PVParameterised*>(physVol);
            if (g4PVParameterised)
            {
                //        G4VSolid *slaveg4VSolid = g4VSolid->Clone();
                //g4LogicalVolume->DestroySlaveG4LogicalVolume(g4LogicalVolume);
                g4LogicalVolume->TerminateWorker(g4LogicalVolume);
                //        delete slaveg4VSolid;
            }
            else
            {
                //g4LogicalVolume->DestroySlaveG4LogicalVolume(g4LogicalVolume);
                g4LogicalVolume->TerminateWorker(g4LogicalVolume);
            }
        }
        else
        {
            //g4LogicalVolume->DestroySlaveG4LogicalVolume(g4LogicalVolume);
            g4LogicalVolume->TerminateWorker(g4LogicalVolume);
        }
    }
    
    //const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
    
    //size_t nmat = theMaterialTable->size();
    //size_t i;
    //for(i=0; i<nmat; i++) {
    //    ((*theMaterialTable)[i])->DestroySlaveG4Material();
    //}
    
    const_cast<G4LVManager&>(G4LogicalVolume::GetSubInstanceManager()).FreeSlave();
    const_cast<G4PVManager&>(G4VPhysicalVolume::GetSubInstanceManager()).FreeSlave();
    const_cast<G4PVRManager&>(G4PVReplica::GetSubInstanceManager()).FreeSlave();
    const_cast<G4PDefManager&>(G4ParticleDefinition::GetSubInstanceManager()).FreeSlave();
    const_cast<G4RegionManager&>(G4Region::GetSubInstanceManager()).FreeSlave();
    //const_cast<G4PVecManager&>(G4PhysicsVector::GetSubInstanceManager()).FreeSlave();
    const_cast<G4DecayChannelManager&>(G4VDecayChannel::GetSubInstanceManager()).FreeSlave();
    const_cast<G4PlSideManager&>(G4PolyconeSide::GetSubInstanceManager()).FreeSlave();
    const_cast<G4PhSideManager&>(G4PolyhedraSide::GetSubInstanceManager()).FreeSlave();
    const_cast<G4VUPLManager&>(G4VUserPhysicsList::GetSubInstanceManager()).FreeSlave();
    const_cast<G4VPCManager&>(G4VPhysicsConstructor::GetSubInstanceManager()).FreeSlave();
    const_cast<G4VMPLManager&>(G4VModularPhysicsList::GetSubInstanceManager()).FreeSlave();
}

