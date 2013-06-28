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
#include <pthread.h>
#include "tbbhelper.hh"

#include "G4PVReplica.hh"
#include "G4PVParameterised.hh"
#include "G4VPVParameterisation.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VSolid.hh"
#include "G4Material.hh"
#include "G4VDecayChannel.hh"
#include "G4PolyconeSide.hh"
#include "G4PolyhedraSide.hh"
#include "G4LogicalVolume.hh"
#include "G4Region.hh"
#include "G4ParticleDefinition.hh"

#include "G4VUserPhysicsList.hh"
#include "G4VPhysicsConstructor.hh" 
#include "G4VModularPhysicsList.hh"

#include "G4AutoLock.hh"

G4Mutex solidclone= G4MUTEX_INITIALIZER;

void tbbSlaveBuildGeometryAndPhysicsVector()
{
  TBBMSG("Copying geometry end physics vector");

  // Implementation copied from
  // void G4WorkerThread::BuildGeometryAndPhysicsVector()
  //   John Apostolakis, May 30, 2013
  
    //Geometry related, split classes mechanism: 
    //   instantiate sub-instances for this thread
    const_cast<G4LVManager&>(G4LogicalVolume::GetSubInstanceManager()).
       SlaveCopySubInstanceArray();
    const_cast<G4PVManager&>(G4VPhysicalVolume::GetSubInstanceManager()).
       SlaveCopySubInstanceArray();
    const_cast<G4PVRManager&>(G4PVReplica::GetSubInstanceManager()).
       SlaveCopySubInstanceArray();
    const_cast<G4RegionManager&>(G4Region::GetSubInstanceManager()).
       SlaveInitializeSubInstance();
    const_cast<G4PlSideManager&>(G4PolyconeSide::GetSubInstanceManager()).
       SlaveInitializeSubInstance();
    const_cast<G4PhSideManager&>(G4PolyhedraSide::GetSubInstanceManager()).
       SlaveInitializeSubInstance();

    //Physics related
    const_cast<G4DecayChannelManager&>(G4VDecayChannel::
                                       GetSubInstanceManager())
                                      .NewSubInstances();
    // const_cast<G4PVecManager&>(G4PhysicsVector::GetSubInstanceManager()).
       // SlaveInitializeSubInstance();

     // Particle definition - has pointer to the process manager
    const_cast<G4PDefManager&>(G4ParticleDefinition::GetSubInstanceManager()).
       NewSubInstances();

    // Physics list related
    const_cast<G4VUPLManager&>(G4VUserPhysicsList::GetSubInstanceManager()).
       NewSubInstances();
    const_cast<G4VPCManager&>(G4VPhysicsConstructor::GetSubInstanceManager()).
       NewSubInstances();
    const_cast<G4VMPLManager&>(G4VModularPhysicsList::GetSubInstanceManager()).
       SlaveCopySubInstanceArray();
    
    G4PhysicalVolumeStore* physVolStore = G4PhysicalVolumeStore::GetInstance();
    for (size_t ip=0; ip<physVolStore->size(); ip++)
    {
        G4VPhysicalVolume* physVol = (*physVolStore)[ip];
        G4LogicalVolume *g4LogicalVolume = physVol->GetLogicalVolume();
        //use shadow pointer
        G4VSolid *g4VSolid = g4LogicalVolume->GetMasterSolid();
        G4PVReplica *g4PVReplica = 0;
        g4PVReplica =  dynamic_cast<G4PVReplica*>(physVol);
        if (!g4PVReplica)
        {
            // Placement volume
          
            g4LogicalVolume->InitialiseWorker(g4LogicalVolume,g4VSolid,0);
        }
        else
        {          
            g4PVReplica->InitialiseWorker(g4PVReplica);
            if( ! g4PVReplica->IsParameterised() )
            {
               g4LogicalVolume->InitialiseWorker(g4LogicalVolume,
                                                 g4VSolid,0);
            }
            else
            {
               G4PVParameterised *g4PVParameterised
                    =  dynamic_cast<G4PVParameterised*>(physVol);
               if (!g4PVParameterised)
               {
                 G4Exception("tbbSlaveBuildGeometryAndPhysicsVector", 
                   "Runtime Error PV01",
                   FatalException,
                  "Cannot find Parameterisation for G4PVParameterised object."
                    );

                 // Each G4PVParameterised instance, has associated with it
                 //  at least one solid for each worker thread.
                 // *Simple* Parameterisations have a single type of solid, 
                 // and the pointer points to the same instance of a solid 
                 // during the simulation.
                 // For this case, it is possible to adapt automatically to
                 // multi-threading, simply by cloning the solid - so long
                 // as all solids support the Clone() method.

                // Check whether it is a simple parameterisation or not
                G4VPVParameterisation *param= 
                   g4PVParameterised->GetParameterisation();
                unsigned int numCopies= g4PVParameterised->GetMultiplicity();
                unsigned int numDifferent= 0;
                for( unsigned int i=0; i< numCopies; i++)
                {
                   G4VSolid *solidChk= param->ComputeSolid(i,  physVol);
                   if( solidChk != g4VSolid)
                   {
                      numDifferent++;
                   }
                }
                if( param->IsNested() || numDifferent>0 )
                {
                   G4ExceptionDescription ed;
                    if( param->IsNested() )
                       ed << " Parameterisation is Nested. " << G4endl;
                    if( numDifferent > 0)
                    {
                       ed << " Parameterisation is implemented using "
                          << " several instances of Solids "
                          << " - potentially to support different " 
                          << " types of solids. " << G4endl;
                    }
                    ed << "  The current implementation does not support "
                      << " this type of Parameterisation" << G4endl;
                    G4Exception("tbbSlaveBuildGeometryAndPhysicsVector", 
                                "GeometryNotSupportedInMT-01", 
                                FatalException, ed);
                }
                  
                // Threads may clone some solids simultaneously. Those cloned 
                // solids will be registered into a shared solid store
                //  (which is a C++ container). 
                // Need a lock to guarantee thread safety. 
                G4AutoLock aLock(&solidclone);
                G4VSolid *slaveg4VSolid = g4VSolid->Clone();
                aLock.unlock();
                if( slaveg4VSolid )
                {
                  g4LogicalVolume->InitialiseWorker(g4LogicalVolume,
                                                    slaveg4VSolid,0);
                }else{
                  // In the case that not all solids support(ed) the Clone()
                  // method, we do similar thing here to dynamically cast
                  // and then get the clone method.
                  G4ExceptionDescription ed;
                  ed << " A solid lacks the Clone() method - or Clone() failed."
                     << G4endl;
                  ed << "   Type of solid: " << g4VSolid->GetEntityType() 
                     << G4endl;
                  ed << "   Parameters: " << *g4VSolid << G4endl;
                  G4Exception("tbbSlaveBuildGeometryAndPhysicsVector", 
                              "TBB-Build-Geometry001",
                              FatalException, ed);
                }
               }
        }
    }

    TBBMSG("Copying geometry end physics vector. Done!"); 
  }
}

void tbbSlaveDestroyGeometryAndPhysicsVector()
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
    
    const_cast<G4LVManager&>(G4LogicalVolume::GetSubInstanceManager()).
       FreeSlave();
    const_cast<G4PVManager&>(G4VPhysicalVolume::GetSubInstanceManager()).
       FreeSlave();
    const_cast<G4PVRManager&>(G4PVReplica::GetSubInstanceManager()).
       FreeSlave();
    const_cast<G4PDefManager&>(G4ParticleDefinition::
                               GetSubInstanceManager()).
       FreeSlave();
    const_cast<G4RegionManager&>(G4Region::GetSubInstanceManager()).
       FreeSlave();
    // const_cast<G4PVecManager&>(G4PhysicsVector::
                               // GetSubInstanceManager()).FreeSlave();
    const_cast<G4DecayChannelManager&>(G4VDecayChannel::
                                       GetSubInstanceManager()).
       FreeSlave();
    const_cast<G4PlSideManager&>(G4PolyconeSide::
                                 GetSubInstanceManager()).FreeSlave();
    const_cast<G4PhSideManager&>(G4PolyhedraSide::
                                 GetSubInstanceManager()).FreeSlave();
    const_cast<G4VUPLManager&>(G4VUserPhysicsList::
                                 GetSubInstanceManager()).FreeSlave();
    const_cast<G4VPCManager&>(G4VPhysicsConstructor::
                              GetSubInstanceManager()).FreeSlave();
    const_cast<G4VMPLManager&>(G4VModularPhysicsList::
                               GetSubInstanceManager()).
       FreeSlave();
}

#include <iostream>
// #include "G4MTGetTid.hh"

#ifndef  GETTID_hh
#include <sys/types.h>
extern "C" {
pid_t gettid();
}
#endif

// #include "G4AutoLock.hh"
G4Mutex debugouttbb = G4MUTEX_INITIALIZER;
void tbbdebugmsg( const char* file, int where, const char* msg) {
  G4AutoLock l(&debugouttbb);
  std::cout<<"tid:"<< gettid() <<";"<<file<<":"<<where<<" -> "<<msg<<std::endl;
}

#include <sys/syscall.h>
pid_t gettid()
{
  pid_t mytid = syscall(SYS_thread_selfid);
  return mytid;
}
