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
// G4GeometryWorkspace
//
// Class Description:
//
// Class managing the per-thread state of the geometry, spanning those
// which have a per-thread state and their dependents.
// In particular it:
//       - owns the arrays that implement 'split' classes 
//       - owns classes/objects which are owned by the split classes.
// The classes/objects affected are:
//       - 'split' classes part of its state is per-thread,
//       - per-thread objects, in particular those which are owned 
//         by the split classes.
// Goal: Take ownership and control of per-thread state of 
//       classes to work with multi-threading. 

// Authors: John Apostolakis (CERN), Andrea Dotti (SLAC), July 2013
// --------------------------------------------------------------------
#ifndef G4GEOMETRYWORKSPACE_HH
#define G4GEOMETRYWORKSPACE_HH

#include "G4TWorkspacePool.hh"
#include "G4PVReplica.hh"
#include "G4PVParameterised.hh"
#include "G4VPVParameterisation.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4VSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4Region.hh"

class G4GeometryWorkspace
{
  public:
 
    using pool_type = G4TWorkspacePool<G4GeometryWorkspace>;

    G4GeometryWorkspace();
   ~G4GeometryWorkspace();

    void UseWorkspace();     // Take ownership
    void ReleaseWorkspace(); // Release ownership
    void DestroyWorkspace(); // Release ownership and destroy

    void InitialiseWorkspace();
      // To be called at start of each run (especially 2nd and further runs)

    inline void   SetVerbose(G4bool v) { fVerbose = v; } 
    inline G4bool GetVerbose()  { return fVerbose; } 

    static pool_type* GetPool();
  
  protected:  // Implementation methods

    void InitialisePhysicalVolumes();
    G4bool CloneParameterisedSolids( G4PVParameterised* paramVol );
    G4bool CloneReplicaSolid( G4PVReplica* );
  
  private:    // Helper pointers - can be per instance or shared

    G4LVManager*     fpLogicalVolumeSIM;
    G4PVManager*     fpPhysicalVolumeSIM;
    G4PVRManager*    fpReplicaSIM;
    G4RegionManager* fpRegionSIM;
  
    // Per Instance variables
    // NOTE: the ownership of the Data Arrays is IN this object
    // Store SubInstanceManager object pointers (SIM pointers)

    G4LVData*      fLogicalVolumeOffset;
    G4PVData*      fPhysicalVolumeOffset;
    G4ReplicaData* fReplicaOffset;       
    G4RegionData*  fRegionOffset;        

    G4bool         fVerbose = false;
};

#endif
