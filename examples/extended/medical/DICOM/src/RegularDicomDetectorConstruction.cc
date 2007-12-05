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
// History:
//	Pedro Arce  
//
//*******************************************************

#include "globals.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"

#include "RegularDicomDetectorConstruction.hh"
#include "G4PhantomParameterisation.hh"

RegularDicomDetectorConstruction::RegularDicomDetectorConstruction() : DicomDetectorConstruction()
{
}

RegularDicomDetectorConstruction::~RegularDicomDetectorConstruction()
{
}


//-------------------------------------------------------------
void RegularDicomDetectorConstruction::ConstructPatient()
{
  //---- Extract number of voxels and voxel dimensions
  G4int nVoxelX = fZSliceHeaderMerged->GetNoVoxelX();
  G4int nVoxelY = fZSliceHeaderMerged->GetNoVoxelY();
  G4int nVoxelZ = fZSliceHeaderMerged->GetNoVoxelZ();

  G4double voxelDimX = fZSliceHeaderMerged->GetVoxelHalfX();
  G4double voxelDimY = fZSliceHeaderMerged->GetVoxelHalfY();
  G4double voxelDimZ = fZSliceHeaderMerged->GetVoxelHalfZ();
#ifdef G4VERBOSE
  G4cout << " nVoxelX " << nVoxelX << " voxelDimX " << voxelDimX <<G4endl;
  G4cout << " nVoxelY " << nVoxelY << " voxelDimY " << voxelDimY <<G4endl;
  G4cout << " nVoxelZ " << nVoxelZ << " voxelDimZ " << voxelDimZ <<G4endl;
  G4cout << " totalPixels " << nVoxelX*nVoxelY*nVoxelZ <<  G4endl;
#endif

  //----- Create parameterisation 
  G4PhantomParameterisation* param = new G4PhantomParameterisation();

  //----- Set voxel dimensions
  param->SetVoxelDimensions( voxelDimX, voxelDimY, voxelDimZ );

  //----- Set number of voxels 
  param->SetNoVoxel( nVoxelX, nVoxelY, nVoxelZ );

  //----- Set list of materials
  param->SetMaterials( fMaterials ); 

  //----- Set list of material indices: for each voxel it is a number that correspond to the index of its material in the vector of materials defined above
  param->SetMaterialIndices( fMateIDs );

  //----- Define voxel logical volume
  G4Box* voxel_solid = new G4Box( "Voxel", voxelDimX, voxelDimY, voxelDimZ);
  G4LogicalVolume* voxel_logic = new G4LogicalVolume(voxel_solid,fMaterials[0],"VoxelLogical",0,0,0); // material is not relevant, it will be changed by the ComputeMaterial method of the parameterisation

  //----- Define the volume that contains all the voxels
  G4Box* container_solid = new G4Box("PhantomContainer",nVoxelX*voxelDimX,nVoxelY*voxelDimY,nVoxelZ*voxelDimZ);
  G4LogicalVolume* container_logic = 
    new G4LogicalVolume( container_solid, 
			 fMaterials[0],  //the material is not important, it will be fully filled by the voxels
			 "PhantomContainer", 
			 0, 0, 0 );
  //--- Place it on the world
  G4double offsetX = (fZSliceHeaderMerged->GetMaxX() + fZSliceHeaderMerged->GetMinX() ) /2.;
  G4double offsetY = (fZSliceHeaderMerged->GetMaxY() + fZSliceHeaderMerged->GetMinY() ) /2.;
  G4double offsetZ = (fZSliceHeaderMerged->GetMaxZ() + fZSliceHeaderMerged->GetMinZ() ) /2.;
  G4ThreeVector posCentreVoxels(offsetX,offsetY,offsetZ);
#ifdef G4VERBOSE
  G4cout << " placing voxel container volume at " << posCentreVoxels << G4endl;
#endif
  G4VPhysicalVolume * container_phys = 
    new G4PVPlacement(0,  // rotation
		      posCentreVoxels,
		      container_logic,     // The logic volume
		      "PhantomContainer",  // Name
		      world_logic,  // Mother
		      false,           // No op. bool.
		      1);              // Copy number

  //--- This physical volume should be assigned as the container volume of the parameterisation
  param->BuildContainerSolid(container_phys);

  //--- Assure yourself that the voxels are completely filling the container volume
  param->CheckVoxelsFillContainer( container_solid->GetXHalfLength(), 
                                   container_solid->GetYHalfLength(), 
                                   container_solid->GetZHalfLength() );


  //----- The G4PVParameterised object that uses the created parameterisation should be placed in the container logical volume
  G4PVParameterised * patient_phys = new G4PVParameterised("Patient",voxel_logic,container_logic,
			kXAxis, nVoxelX*nVoxelY*nVoxelZ, param);
  // if axis is set as kUndefined instead of kXAxis, GEANT4 will do an smart voxel optimisation (not needed if G4RegularNavigation is used)

  //----- Set this physical volume as having a regular structure of type 1, so that G4RegularNavigation is used
  patient_phys->SetRegularStructureId(1); // if not set, G4VoxelNavigation will be used instead 
}

