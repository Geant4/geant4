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
//
/// \file DicomRegularDetectorConstruction.cc
/// \brief Implementation of the DicomRegularDetectorConstruction clas
//
// History:
//        Pedro Arce  
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

#include "DicomRegularDetectorConstruction.hh"
#include "DicomPhantomParameterisationColour.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomRegularDetectorConstruction::DicomRegularDetectorConstruction() 
 : DicomDetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DicomRegularDetectorConstruction::~DicomRegularDetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DicomRegularDetectorConstruction::ConstructPhantom()
{
#ifdef G4VERBOSE
  G4cout << "DicomRegularDetectorConstruction::ConstructPhantom " << G4endl;
#endif

  //----- Create parameterisation 
  DicomPhantomParameterisationColour* param = 
    new DicomPhantomParameterisationColour();

  //----- Set voxel dimensions
  param->SetVoxelDimensions( fVoxelHalfDimX, fVoxelHalfDimY, fVoxelHalfDimZ );

  //----- Set number of voxels 
  param->SetNoVoxel( fNVoxelX, fNVoxelY, fNVoxelZ );

  //----- Set list of materials
  param->SetMaterials( fMaterials ); 

  //----- Set list of material indices: for each voxel it is a number that
  // correspond to the index of its material in the vector of materials
  // defined above
  param->SetMaterialIndices( fMateIDs );

  //----- Define voxel logical volume
  G4Box* voxel_solid = 
    new G4Box( "Voxel", fVoxelHalfDimX, fVoxelHalfDimY, fVoxelHalfDimZ);
  G4LogicalVolume* voxel_logic = 
    new G4LogicalVolume(voxel_solid,fMaterials[0],"VoxelLogical",
                             0,0,0);
  // material is not relevant, it will be changed by the
  // ComputeMaterial method of the parameterisation

    voxel_logic->SetVisAttributes(
                     new G4VisAttributes(G4VisAttributes::GetInvisible()));
    
  //--- Assign the fContainer volume of the parameterisation
  param->BuildContainerSolid(fContainer_phys);

  //--- Assure yourself that the voxels are completely filling the 
  // fContainer volume
  param->CheckVoxelsFillContainer( fContainer_solid->GetXHalfLength(), 
                                   fContainer_solid->GetYHalfLength(), 
                                   fContainer_solid->GetZHalfLength() );

  //----- The G4PVParameterised object that uses the created parameterisation
  // should be placed in the fContainer logical volume
  G4PVParameterised * phantom_phys = 
    new G4PVParameterised("phantom",voxel_logic,fContainer_logic,
                          kXAxis, fNVoxelX*fNVoxelY*fNVoxelZ, param);
  // if axis is set as kUndefined instead of kXAxis, GEANT4 will 
  //  do an smart voxel optimisation 
  // (not needed if G4RegularNavigation is used)

  //----- Set this physical volume as having a regular structure of type 1, 
  // so that G4RegularNavigation is used
  phantom_phys->SetRegularStructureId(1); // if not set, G4VoxelNavigation
  //will be used instead 

  SetScorer(voxel_logic);
}
