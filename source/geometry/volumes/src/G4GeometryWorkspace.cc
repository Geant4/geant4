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

#include "G4GeometryWorkspace.hh"

#include "G4PVReplica.hh"
#include "G4PVParameterised.hh"
#include "G4VPVParameterisation.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVReplica.hh"
#include "G4Region.hh"

#include "G4AutoLock.hh"


G4GeometryWorkspace::G4GeometryWorkspace()
{
  
  fpLogicalVolumeSIM=
      &const_cast<G4LVManager&>(G4LogicalVolume::GetSubInstanceManager());
  fpPhysicalVolumeSIM=
      &const_cast<G4PVManager&>( G4VPhysicalVolume::GetSubInstanceManager() );
  fpReplicaSIM=
      &const_cast<G4PVRManager&>(G4PVReplica::GetSubInstanceManager());
  fpRegionSIM=
      &const_cast<G4RegionManager&>(G4Region::GetSubInstanceManager());

  // Create a work area for Logical Volumes in this thread - then capture its address
  InitialiseWorkspace();
  
  // fpLogicalVolumeSIM->SlaveCopySubInstanceArray();
  fLogicalVolumeOffset= fpLogicalVolumeSIM->GetOffset();

  // fpPhysicalVolumeSIM->SlaveCopySubInstanceArray();
  fPhysicalVolumeOffset= fpPhysicalVolumeSIM->GetOffset();

  // fpReplicaSIM->SlaveCopySubInstanceArray();
  fReplicaOffset= fpReplicaSIM->GetOffset();

  // fpRegionSIM->SlaveInitializeSubInstance();
  fRegionOffset= fpRegionSIM->GetOffset();
}

G4GeometryWorkspace::~G4GeometryWorkspace()
{
  
}

//  Static methods 
//      For with current (original) G4WorkerThread -- which uses static methods

void
G4GeometryWorkspace::UseWorkspace()
{
  G4cout << "G4GeometryWorkspace::UseWorkspace: Copying geometry - Start " << G4endl;

  // Implementation copied from  G4WorkerThread::BuildGeometryAndPhysicsVector()
  //  and improved for G4PVParamaterised
  //   John Apostolakis, May 30, 2013
  
  //Geometry related, split classes mechanism: instantiate sub-instance for this thread
  fpLogicalVolumeSIM->UseWorkArea(fLogicalVolumeOffset);
  fpPhysicalVolumeSIM->UseWorkArea(fPhysicalVolumeOffset);
  
  fpReplicaSIM->UseWorkArea(fReplicaOffset);
  fpRegionSIM->UseWorkArea(fRegionOffset);

  InitialisePhysicalVolumes();  // Not sure if this is needed - done already ?
}


void G4GeometryWorkspace::ReleaseWorkspace()
//  The opposite of Use Workspace - let go of it.
{
  fpLogicalVolumeSIM->UseWorkArea(0);
  fpPhysicalVolumeSIM->UseWorkArea(0);
  
  fpReplicaSIM->UseWorkArea(0);
  fpRegionSIM->UseWorkArea(0);

  InitialisePhysicalVolumes();  // Not sure if this is needed - done already ?
}

G4Mutex solidclone= G4MUTEX_INITIALIZER;

void G4GeometryWorkspace::InitialisePhysicalVolumes()
  {
    G4PhysicalVolumeStore* physVolStore = G4PhysicalVolumeStore::GetInstance();
    for (size_t ip=0; ip<physVolStore->size(); ip++)
    {
        G4VPhysicalVolume* physVol = (*physVolStore)[ip];
        G4LogicalVolume *logicalVol = physVol->GetLogicalVolume();
        //use shadow pointer
        G4VSolid *solid = logicalVol->GetMasterSolid();
        G4PVReplica *g4PVReplica = 0;
        g4PVReplica =  dynamic_cast<G4PVReplica*>(physVol);
        if (!g4PVReplica)
        {
            // Placement volume
            logicalVol->InitialiseWorker(logicalVol,solid,0);
        }
        else
        {          
            //g4PVReplica->SlaveG4PVReplica(g4PVReplica);
            g4PVReplica->InitialiseWorker(g4PVReplica);
            if( ! g4PVReplica->IsParameterised() )
            {
              //logicalVol->SlavelogicalVol(logicalVol, solid, 0);
               logicalVol->InitialiseWorker(logicalVol,solid,0);
            }
            else
            {
               G4PVParameterised *paramVol
                   =  dynamic_cast<G4PVParameterised*>(physVol);
               if (!paramVol)
               {
                  G4Exception("G4GeometryWorkspace::CreateAndUseWorkspace", "Runtime Error PV01",
                              FatalException,
                              "Cannot find Parameterisation for G4PVParameterised object.");
               }
               CloneParameterisedVolume( paramVol );
            }
        }
    }
    G4cout << "G4GeometryWorkspace::CreateAndUseWorkspace: "
           << "Copying geometry - Done!" << G4endl;
}


G4bool G4GeometryWorkspace::CloneParameterisedVolume( G4PVParameterised *paramVol )
{
  // Each G4PVParameterised instance, has associated with it at least one
  // solid for each worker thread.
  // *Simple* Parameterisations have a single type of solid, and the
  // pointer points to the same instance of a solid during the  simulation.
  // For this case, it is possible to adapt automatically to
  // multi-threading, simply by cloning the solid - so long
  // as all solids support the Clone() method.
  
  // Check whether it is a simple parameterisation or not
  G4VPVParameterisation *param= paramVol->GetParameterisation();
  unsigned int numCopies= paramVol->GetMultiplicity();
  unsigned int numDifferent= 0;

  G4LogicalVolume *logicalV= paramVol->GetLogicalVolume();
  // assert( logicalV != 0);
  G4VSolid *solid= logicalV->GetSolid();
  
  for( unsigned int i=0; i< numCopies; i++)
  {
    G4VSolid *solidChk= param->ComputeSolid(i, paramVol);
    if( solidChk != solid)
    {
      numDifferent++;
    }
  }
  if( param->IsNested() || (numDifferent>0) )
  {
    G4ExceptionDescription ed;
    if( param->IsNested() )
      ed << " Parameterisation is Nested. " << G4endl;
    if( numDifferent > 0)
    {
      ed << " Parameterisation is implemented using several instances of Solids "
      << " - potentially to support different types of solids. " << G4endl;
    }
    ed << "  The current implementation does not support "
    << " this type of Parameterisation" << G4endl;
    G4Exception("G4GeometryWorkspace::CloneParameterisedVolume", "GeometryNotSupportedInMT-01",
                FatalException, ed);
  }
  
  // Threads may attempt to clone a solids simultaneously. Those cloned solids will be
  // registered into a shared solid store (C++ container). Need a lock to
  // guarantee thread safety
  G4AutoLock aLock(&solidclone);
  G4VSolid *workerSolid = solid->Clone();
  aLock.unlock();
  if( workerSolid )
  {
    logicalV->InitialiseWorker(logicalV,workerSolid,0);
  }else{
    // In the case that not all solids support(ed) the Clone()
    // method, we do similar thing here to dynamically cast
    // and then get the clone method.
    G4ExceptionDescription ed;
    ed << " ERROR: Unable to initialise geometry for worker node." << G4endl;
    ed << " A solid lacks the Clone() method - or Clone() failed." << G4endl;
    ed << "   Type of solid: " << solid->GetEntityType() << G4endl;
    ed << "   Parameters: " << *solid << G4endl;
    G4Exception(" G4GeometryWorkspace::CloneParameterisedVolume", "MT-BuildGeometry001",
                FatalException, ed);
  }
  return true; // It Worked
}


void
G4GeometryWorkspace::InitialiseWorkspace()
{
  G4cout << "G4GeometryWorkspace::InitialiseWorkspace: Copying geometry - Start " << G4endl;
  
  // Implementation copied from  G4WorkerThread::BuildGeometryAndPhysicsVector()
  //  and improved for G4PVParamaterised
  //   John Apostolakis, May 30, 2013
  
  //Geometry related, split classes mechanism:
  //   Do *NOT* instantiate sub-instance for this thread,
  //     just copy the contents !!
  fpLogicalVolumeSIM->SlaveCopySubInstanceArray();
  fpPhysicalVolumeSIM->SlaveCopySubInstanceArray();
  fpReplicaSIM->SlaveCopySubInstanceArray();
  fpRegionSIM->SlaveInitializeSubInstance();

  InitialisePhysicalVolumes();
  
  G4cout << "G4GeometryWorkspace::CreateAndUseWorkspace: "
  << "Copying geometry - Done!" << G4endl;
}

void G4GeometryWorkspace::DestroyWorkspace()
{
  G4PhysicalVolumeStore* physVolStore = G4PhysicalVolumeStore::GetInstance();
  for (size_t ip=0; ip<physVolStore->size(); ip++)
  {
    G4VPhysicalVolume* physVol = (*physVolStore)[ip];
    G4LogicalVolume *logicalVol = physVol->GetLogicalVolume();
    G4PVReplica *g4PVReplica = 0;
    g4PVReplica =  dynamic_cast<G4PVReplica*>(physVol);
    if (g4PVReplica)
    {
      g4PVReplica->TerminateWorker(g4PVReplica);
      G4PVParameterised *paramVol = 0;
      paramVol =  dynamic_cast<G4PVParameterised*>(physVol);
      if (paramVol)
      {
        // G4VSolid *solid = logicalVol->fSolid;
        logicalVol->TerminateWorker(logicalVol);
        // if( solid->IsClone() ) delete solid;
      }
      else
      {
        logicalVol->TerminateWorker(logicalVol);
      }
    }
    else
    {
      logicalVol->TerminateWorker(logicalVol);
    }
  }
  
  fpLogicalVolumeSIM->FreeSlave();
  fpPhysicalVolumeSIM->FreeSlave();
  fpReplicaSIM->FreeSlave();
  fpRegionSIM->FreeSlave();
}
