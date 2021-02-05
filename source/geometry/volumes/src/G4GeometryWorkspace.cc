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
  if( fVerbose ) 
  { 
    G4cout << "G4GeometryWorkspace::UseWorkspace: Start " << G4endl;
  }

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

  if( fVerbose )
  { 
     G4cout << "G4GeometryWorkspace::UseWorkspace:  End " << G4endl;
  }
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
      logicalVol->InitialiseWorker(logicalVol,solid,0);
    }
    else
    {          
      g4PVReplica->InitialiseWorker(g4PVReplica);
      if( !g4PVReplica->IsParameterised() )
      {
        logicalVol->InitialiseWorker(logicalVol,solid,0);

        // If the replica's solid (in LV) is changed during navigation,
        // it must be thread-private
        //
        CloneReplicaSolid( g4PVReplica );
      }
      else
      {
        G4PVParameterised* paramVol = dynamic_cast<G4PVParameterised*>(physVol);
        if (paramVol == nullptr)
        {
          G4Exception("G4GeometryWorkspace::CreateAndUseWorkspace()",
                      "GeomVol0003", FatalException,
                      "Cannot find Parameterisation for parameterised volume.");
        }
        CloneParameterisedSolids( paramVol );
      }
    }
  }
  if( fVerbose )
  {
    G4cout << "G4GeometryWorkspace::InitialisePhysicalVolumes: "
           << "Copying geometry - Done!" << G4endl;
  }
}

// ----------------------------------------------------------------------
// Create a clone of the solid for this replica in this thread
//
G4bool G4GeometryWorkspace::CloneReplicaSolid( G4PVReplica *replicaPV )
{
  G4LogicalVolume* logicalV = replicaPV->GetLogicalVolume();
  G4VSolid* solid = logicalV->GetSolid();

  G4AutoLock aLock(&mutex_init);
  G4VSolid* workerSolid = solid->Clone();
  aLock.unlock();

  if( workerSolid != nullptr )
  {
    logicalV->InitialiseWorker(logicalV,workerSolid,0);
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
    G4Exception("G4GeometryWorkspace::CloneParameterisedVolume()",
                "GeomVol0003", FatalException, ed);
    return false; 
  }
  return true; // It Worked
}

// ----------------------------------------------------------------------
// Each G4PVParameterised instance has associated with it at least one
// solid for each worker thread.
// *Simple* Parameterisations have a single type of solid, and the
// pointer points to the same instance of a solid during the simulation.
// For this case, it is possible to adapt automatically to
// multi-threading, simply by cloning the solid - so long
// as all solids support the Clone() method.
//
G4bool G4GeometryWorkspace::
CloneParameterisedSolids( G4PVParameterised* paramVol )
{
  // Check whether it is a simple parameterisation or not
  //
  // G4VPVParameterisation *param= paramVol->GetParameterisation();
  // unsigned int numCopies= paramVol->GetMultiplicity();
  // unsigned int numDifferent= 0;

  G4LogicalVolume* logicalV= paramVol->GetLogicalVolume();
  G4VSolid* solid= logicalV->GetSolid();
  
  //  for( unsigned int i=0; i< numCopies; ++i)
  //  {
  //    G4VSolid *solidChk= param->ComputeSolid(i, paramVol);
  //    if( solidChk != solid)
  //    {
  //      ++numDifferent;
  //    }
  //  }
  //  if( numDifferent>0 )
  //  {
  //    G4ExceptionDescription ed;
  //    ed << "ERROR - Parameterisation using several instances of Solids \n"
  //       << "potentially to support different types of solids. \n"
  //       << "Geant4-MT currently does not support this type of \n"
  //       << "parameterisation, sorry !";
  //    G4Exception("G4GeometryWorkspace::CloneParameterisedVolume()",
  //                "GeomVol0001", FatalException, ed);
  //  }
  
  // Threads may attempt to clone a solid simultaneously.
  // Those cloned solids will be registered into a shared solid
  // store (C++ container). Need a lock to guarantee thread safety
  //
  G4AutoLock aLock(&mutex_init);
  G4VSolid *workerSolid = solid->Clone();
  aLock.unlock();
  if( workerSolid != nullptr )
  {
    logicalV->InitialiseWorker(logicalV, workerSolid, nullptr);
  }
  else
  {
    // In the case that not all solids support(ed) the Clone()
    // method, we do similar thing here to dynamically cast
    // and then get the clone method
    //
    G4ExceptionDescription ed;
    ed << "ERROR - Unable to initialise geometry for worker node. \n"
       << "A solid lacks the Clone() method - or Clone() failed. \n"
       << "   Type of solid: " << solid->GetEntityType() << "\n"
       << "   Parameters: " << *solid;
    G4Exception("G4GeometryWorkspace::CloneParameterisedVolume()",
                "GeomVol0003", FatalException, ed);
  }
  return true; // It Worked
}

// ----------------------------------------------------------------------
//
void G4GeometryWorkspace::InitialiseWorkspace()
{
  if( fVerbose )
  {
    G4cout << "G4GeometryWorkspace::InitialiseWorkspace():"
           << " Copying geometry - Start " << G4endl;
  }
  
  // Geometry related, split classes mechanism:
  // Do *NOT* instantiate sub-instance for this thread, just copy the contents!
  //
  fpLogicalVolumeSIM->SlaveCopySubInstanceArray();
  fpPhysicalVolumeSIM->SlaveCopySubInstanceArray();
  fpReplicaSIM->SlaveCopySubInstanceArray();
  fpRegionSIM->SlaveInitializeSubInstance();

  InitialisePhysicalVolumes();
  
  if( fVerbose )
  {
    G4cout << "G4GeometryWorkspace::InitialiseWorkspace: "
           << "Copying geometry - Done!" << G4endl;
  }
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
