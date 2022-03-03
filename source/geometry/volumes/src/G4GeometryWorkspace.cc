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
// G4GeometryWorkspace - implementation
//
// Authors: John Apostolakis (CERN), Andrea Dotti (SLAC), July 2013
// --------------------------------------------------------------------

#include "G4GeometryWorkspace.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4AutoLock.hh"

namespace
{
  G4Mutex mutex_init = G4MUTEX_INITIALIZER;
  G4GeometryWorkspace::pool_type thePool;
}

// ----------------------------------------------------------------------
//
G4GeometryWorkspace::pool_type* G4GeometryWorkspace::GetPool()
{
  return &thePool;
}

// ----------------------------------------------------------------------
//
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

  // Create a work area for Logical Volumes in this thread
  // then capture its address
  //
  InitialiseWorkspace();
  
  fLogicalVolumeOffset = fpLogicalVolumeSIM->GetOffset();

  fPhysicalVolumeOffset = fpPhysicalVolumeSIM->GetOffset();

  fReplicaOffset = fpReplicaSIM->GetOffset();

  fRegionOffset = fpRegionSIM->GetOffset();
}

// ----------------------------------------------------------------------
//
G4GeometryWorkspace::~G4GeometryWorkspace()
{  
}

// ----------------------------------------------------------------------
//
void
G4GeometryWorkspace::UseWorkspace()
{
  // Geometry related, split classes mechanism: instantiate sub-instance
  // for this thread
  //
  fpLogicalVolumeSIM->UseWorkArea(fLogicalVolumeOffset);
  fpPhysicalVolumeSIM->UseWorkArea(fPhysicalVolumeOffset);
  
  fpReplicaSIM->UseWorkArea(fReplicaOffset);
  fpRegionSIM->UseWorkArea(fRegionOffset);

  // When recycling a workspace 
  //   - it must be a lightweight operation, to reuse a valid work area
  //   - so it must NOT Initialise anything!
  // Do not call InitialisePhysicalVolumes();
}

// ----------------------------------------------------------------------
//
void G4GeometryWorkspace::ReleaseWorkspace()
{
  fpLogicalVolumeSIM->UseWorkArea(nullptr);
  fpPhysicalVolumeSIM->UseWorkArea(nullptr);
  
  fpReplicaSIM->UseWorkArea(nullptr);
  fpRegionSIM->UseWorkArea(nullptr);
}

// ----------------------------------------------------------------------
//
void G4GeometryWorkspace::InitialisePhysicalVolumes()
{
  G4PhysicalVolumeStore* physVolStore = G4PhysicalVolumeStore::GetInstance();
  for (std::size_t ip=0; ip<physVolStore->size(); ++ip)
  {
    G4VPhysicalVolume* physVol = (*physVolStore)[ip];
    G4LogicalVolume *logicalVol = physVol->GetLogicalVolume();
 
    // Use shadow pointer
    //
    G4VSolid* solid = logicalVol->GetMasterSolid();
    G4PVReplica* g4PVReplica = dynamic_cast<G4PVReplica*>(physVol);
    if (g4PVReplica == nullptr)
    {
      // Placement volume
      //
      logicalVol->InitialiseWorker(logicalVol,solid,nullptr);
    }
    else
    {          
      g4PVReplica->InitialiseWorker(g4PVReplica);
      logicalVol->InitialiseWorker(logicalVol,solid,nullptr);

      // If the replica's solid (in LV) is changed during navigation,
      // it must be thread-private
      //
      CloneReplicaSolid( g4PVReplica );
    }
  }
}

// ----------------------------------------------------------------------
// Create a clone of the solid for this replica in this thread
//
G4bool G4GeometryWorkspace::CloneReplicaSolid( G4PVReplica* replicaPV )
{
  G4LogicalVolume* logicalV = replicaPV->GetLogicalVolume();
  G4VSolid* solid = logicalV->GetSolid();

  G4AutoLock aLock(&mutex_init);
  G4VSolid* workerSolid = solid->Clone();
  aLock.unlock();

  if( workerSolid != nullptr )
  {
    logicalV->InitialiseWorker(logicalV,workerSolid,nullptr);
  }
  else
  {
    // In the case that not all solids support(ed) the Clone()
    // method, we do similar thing here to dynamically cast
    // and then get the clone method
    //
    G4ExceptionDescription ed;
    ed << "ERROR - Unable to initialise geometry for worker node." << "\n"
       << "A solid lacks the Clone() method - or Clone() failed." << "\n"
       << "   Type of solid: " << solid->GetEntityType() << "\n"
       << "   Parameters: " << *solid;
    G4Exception("G4GeometryWorkspace::CloneReplicaSolid()",
                "GeomVol0003", FatalException, ed);
    return false; 
  }
  return true; // It Worked
}

// ----------------------------------------------------------------------
//
void G4GeometryWorkspace::InitialiseWorkspace()
{
  // Geometry related, split classes mechanism:
  // Do *NOT* instantiate sub-instance for this thread, just copy the contents!
  //
  fpLogicalVolumeSIM->SlaveCopySubInstanceArray();
  fpPhysicalVolumeSIM->SlaveCopySubInstanceArray();
  fpReplicaSIM->SlaveCopySubInstanceArray();
  fpRegionSIM->SlaveInitializeSubInstance();

  InitialisePhysicalVolumes();
}

// ----------------------------------------------------------------------
//
void G4GeometryWorkspace::DestroyWorkspace()
{
  G4PhysicalVolumeStore* physVolStore = G4PhysicalVolumeStore::GetInstance();
  for (std::size_t ip=0; ip<physVolStore->size(); ++ip)
  {
    G4VPhysicalVolume* physVol = (*physVolStore)[ip];
    G4LogicalVolume* logicalVol = physVol->GetLogicalVolume();
    G4PVReplica* g4PVReplica = dynamic_cast<G4PVReplica*>(physVol);
    if (g4PVReplica != nullptr)
    {
      g4PVReplica->TerminateWorker(g4PVReplica);
    }
    logicalVol->TerminateWorker(logicalVol);
  }

  // Threads may attempt to free memory simultaneously.
  // Need a lock to guarantee thread safety
  //
  G4AutoLock aLock(&mutex_init);
  fpLogicalVolumeSIM->FreeSlave();
  fpPhysicalVolumeSIM->FreeSlave();
  fpReplicaSIM->FreeSlave();
  fpRegionSIM->FreeSlave();
  aLock.unlock();
}
