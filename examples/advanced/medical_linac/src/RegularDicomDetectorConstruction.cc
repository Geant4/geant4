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
//
// The code is part of the DICOM extended example and it was modified by :
//	^Claudio Andenna  claudio.andenna@ispesl.it, claudio.andenna@iss.infn.it
//      *Barbara Caccia barbara.caccia@iss.it
//      with the support of Pablo Cirrone (LNS, INFN Catania Italy)
//	with the contribute of Alessandro Occhigrossi*
//
// ^INAIL DIPIA - ex ISPESL and INFN Roma, gruppo collegato Sanità, Italy
// *Istituto Superiore di Sanità and INFN Roma, gruppo collegato Sanità, Italy
//  Viale Regina Elena 299, 00161 Roma (Italy)
//  tel (39) 06 49902246
//  fax (39) 06 49387075
//
// more information:
// http://g4advancedexamples.lngs.infn.it/Examples/medical-linac
//
//*******************************************************//


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
#include "DicomPhantomParameterisationColour.hh"

RegularDicomDetectorConstruction::RegularDicomDetectorConstruction() : DicomDetectorConstruction()
{
}

RegularDicomDetectorConstruction::~RegularDicomDetectorConstruction()
{
}

//-------------------------------------------------------------
void RegularDicomDetectorConstruction::ConstructPatient()
{
#ifdef G4VERBOSE
  G4cout << "RegularDicomDetectorConstruction::ConstructPatient " << G4endl;
#endif

  //----- Create parameterisation 
  DicomPhantomParameterisationColour* param = new DicomPhantomParameterisationColour(this->getDicomDirectory(), this->getDicomColorMap());

  //----- Set voxel dimensions
  param->SetVoxelDimensions( voxelHalfDimX, voxelHalfDimY, voxelHalfDimZ );

  //----- Set number of voxels 
  param->SetNoVoxel( nVoxelX, nVoxelY, nVoxelZ );

  //----- Set list of materials
  param->SetMaterials( fMaterials ); 

  //----- Set list of material indices: for each voxel it is a number that correspond to the index of its material in the vector of materials defined above
  param->SetMaterialIndices( fMateIDs );

  //----- Define voxel logical volume
  G4Box* voxel_solid = new G4Box( "Voxel", voxelHalfDimX, voxelHalfDimY, voxelHalfDimZ);
  G4LogicalVolume* voxel_logic = new G4LogicalVolume(voxel_solid,fMaterials[0],"VoxelLogical",0,0,0); // material is not relevant, it will be changed by the ComputeMaterial method of the parameterisation

  //--- Assign the container volume of the parameterisation
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

  SetScorer(voxel_logic);
		G4VisAttributes* simpleAlSVisAtt= new G4VisAttributes(G4Colour::White());
		simpleAlSVisAtt->SetVisibility(false);
		simpleAlSVisAtt->SetForceWireframe(false);
// 		simpleAlSVisAtt->SetVisibility(true);
// 		simpleAlSVisAtt->SetForceWireframe(true);
		container_logic->SetVisAttributes(simpleAlSVisAtt);

}

