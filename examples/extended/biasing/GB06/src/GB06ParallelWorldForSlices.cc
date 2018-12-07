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
/// \file GB06/src/GB06ParallelWorldForSlices.cc
/// \brief Implementation of the GB06ParallelWorldForSlices class
//
//
#include "GB06ParallelWorldForSlices.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

#include "GB06BOptrSplitAndKillByImportance.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB06ParallelWorldForSlices::GB06ParallelWorldForSlices(G4String worldName)
  : G4VUserParallelWorld(worldName)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB06ParallelWorldForSlices::~GB06ParallelWorldForSlices()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB06ParallelWorldForSlices::Construct()
{
  // -- Inform about construction:
  // -- (fWorldName is a protected data member of the base parallel world class)
  G4cout << "Parallel World `" << fWorldName << "' constructed." << G4endl;

  // -------------------------
  //  Build parallel geometry:
  // -------------------------
  
  // -- Obtain clone of mass geometry world from GetWorld() base class utility:
  G4VPhysicalVolume* physicalParallelWorld = GetWorld();
  G4LogicalVolume*    logicalParallelWorld = physicalParallelWorld->GetLogicalVolume();


  // -- We overlay a sliced geometry on top of the block of concrete in the mass geometry
  // -- (ie, in the detector construction class), using the same dimensions.
  // -- [Note that this is a choice : we can use different dimensions and shapes, creating
  // -- a new solid for that.]
  // -- For this we:
  // --     - 1) get back the solid used to create the concrete shield;
  // --     - 2) create a new logical volume of same shape than the shield and we place
  // --          inside the slices
  // --     - 3) place the sliced structure, using the placement of the physical volume of
  // --          the concrete shield
  // -- In all this construction, no materials are used, as only the volumes boundaries
  // -- are of interest. Note that the absence of materials is only possible in parallel
  // -- geometries.

  
  // -- 1) get back the solid used to create the concrete shield:
  //       ------------------------------------------------------
  
  // -- get back the logical volume of the shield, using its name:
  G4LogicalVolume* shieldLogical =
    G4LogicalVolumeStore::GetInstance()->GetVolume("shield.logical");
  
  // -- get back the solid, a G4box in this case. We cast the pointer to access later on
  // -- the G4Box class specific methods:
  G4Box* shieldSolid = (G4Box*) shieldLogical->GetSolid();
  
  // -- we now re-create a logical volume for the mother volume of the slices:
  G4LogicalVolume* motherForSlicesLogical =
    new G4LogicalVolume(shieldSolid,                 // its solid
                        nullptr,                     // no material
                        "motherForSlices.logical");  // its name


  
  // -- 2) new logical volume of same shape than the shield and place inside the slices:
  //       -----------------------------------------------------------------------------
  
  // -- We create now the slices; we choose 20 slices:
  const G4int nSlices(20);
  // -- the solid for slices:
  G4double halfSliceZ = shieldSolid->GetZHalfLength() / nSlices;
  G4Box*             sliceSolid = new G4Box("slice.solid",
                                            shieldSolid->GetXHalfLength(),
                                            shieldSolid->GetYHalfLength(),
                                            halfSliceZ                    );
  
  // -- the logical volume for slices:
  G4LogicalVolume* sliceLogical = new G4LogicalVolume(sliceSolid,        // its solid
                                                      nullptr,           // no material
                                                      "slice.logical");  // its name
  
  // -- we use a replica, to place the 20 slices in one go, along the Z axis:
  new G4PVReplica( "slice.physical",        // its name
                   sliceLogical,            // its logical volume
                   motherForSlicesLogical,  // its mother volume
                   kZAxis,                  // axis of replication
                   nSlices,                 // number of replica
                   2*halfSliceZ);           // width of replica
  
  
  // -- 3) place the sliced structure, using the concrete shield placement:
  //       ----------------------------------------------------------------
  
  // -- get back the physical volume of the shield, using its name:
  // -- (note that we know we have only one physical volume with this name. If we had
  // -- several, we should loop by ourselves on the store which is of
  // -- std::vector<G4VPhysicalVolume*> type.)
  G4VPhysicalVolume*
    shieldPhysical = G4PhysicalVolumeStore::GetInstance()->GetVolume("shield.physical");
  
  // -- get back the translation
  // -- (we don't try to get back the rotation, we know we used nullptr):
  G4ThreeVector translation = shieldPhysical->GetObjectTranslation();
  
  // -- finally, we place the sliced structure:
  new G4PVPlacement( nullptr,                           // no rotation
                     translation,                       // translate as for the shield
                     motherForSlicesLogical,            // its logical volume
                     "motherForSlices.physical",        // its name
                     logicalParallelWorld,              // its mother  volume
                     false,                             // no boolean operation
                     0);                                // copy number

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB06ParallelWorldForSlices::ConstructSD()
{
  // -- Create the biasing operator:
  auto biasingOperator = new GB06BOptrSplitAndKillByImportance("neutron","parallelOptr");
  // -- Tell it it is active for this parallel geometry, passing the world
  // -- volume of this geometry :
  biasingOperator->SetParallelWorld( GetWorld() );

  // -- Attach to the logical volume where the biasing has to be applied:
  auto slice = G4LogicalVolumeStore::GetInstance()->GetVolume("slice.logical");
  biasingOperator->AttachTo(slice);

  // -- Create a simple "volume importance" map, linking replica numbers to importances:
  //    --------------------------------------------------------------------------------
  // -- we define the map as going from an importance to 2*importance when going from
  // -- a slice to the next one, in the Z direction.
  // -- Get back the replica of slices:
  G4PVReplica* slicePhysical =
    (G4PVReplica*)(G4PhysicalVolumeStore::GetInstance()->GetVolume("slice.physical"));
  G4int nReplica = slicePhysical->GetMultiplicity();
  // -- We use and fill the map we defined in the biasing operator:
  G4int importance = 1;
  for ( G4int iReplica = 0 ; iReplica < nReplica ; iReplica++ )
    {
      (biasingOperator->GetImportanceMap())[ iReplica ] = importance;
      importance *= 2;
    }
  
}
