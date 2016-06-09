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
// Authors: S. Guatelli and M. G. Pia, INFN Genova, Italy
// 
// Based on code developed by the undergraduate student G. Guerrieri 
// Note: this is a preliminary beta-version of the code; an improved 
// version will be distributed in the next Geant4 public release, compliant
// with the design in a forthcoming publication, and subject to a 
// design and code review.
//
#include "G4ParameterisedBreast.hh"
#include "globals.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Tubs.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4HumanPhantomMaterial.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4HumanPhantomROGeometry.hh"
#include "G4HumanPhantomSDBreast.hh"
G4ParameterisedBreast::G4ParameterisedBreast():breastSD(0)
{
}

G4ParameterisedBreast::~G4ParameterisedBreast()
{


}

G4VPhysicalVolume* G4ParameterisedBreast::ConstructBreast(G4VPhysicalVolume* mother, G4String sex, G4bool sensitivity)
{
  
 G4cout << "ConstructParameterisedBreast for " << sex << G4endl;
 
 G4HumanPhantomMaterial* material = new G4HumanPhantomMaterial();
 G4Material* adipose = material -> GetMaterial("adipose");
 G4Material* adipose_glandular = material -> GetMaterial("adipose_glandular");
 delete material;

  G4double rmin = 0.* cm;
  G4double rmax = 6. * cm;
  G4double zz = 6. * cm;
  G4double startPhi = 0.* degree;
  G4double spanningPhi = 180. * degree;

  G4Tubs* breast = new G4Tubs("out_breast", rmin, rmax,
			      zz/2., startPhi, spanningPhi);

  G4LogicalVolume* breast_log =  new G4LogicalVolume(breast,
						   adipose,
						   "Breast",
						   0, 0, 0);  
  rmax = 5.5 *cm;
  zz = 5. *cm;
  G4Tubs* innerBreast = new G4Tubs("inner_breast",
				   rmin, rmax,
				   zz/2., startPhi, spanningPhi);
  
  G4LogicalVolume* innerBreast_log = new G4LogicalVolume(innerBreast,
							adipose_glandular,
							"InnerBreast",
							 0, 0, 0);

  G4RotationMatrix* matrix = new G4RotationMatrix();
  matrix -> rotateX(-90.* degree);
  matrix -> rotateY(180.* degree);
  matrix -> rotateZ(-18. * degree);


  G4VPhysicalVolume* physBreast = new G4PVPlacement(matrix,G4ThreeVector(8.63*cm, 46.87* cm, 8.4854 *cm),
     			       "physicalParameterisedBreast",
  		       breast_log,
  		       mother,
  		       false,
  		       0);

 G4VPhysicalVolume* physInnerBreast = new G4PVPlacement(0,0,
      			       "innerBreast",
  			       innerBreast_log,
			       physBreast,
			       false,
			       0);


 // Parameterisation with Replicas
 /*
  G4double rmin_voxelz = 0.* cm;
  G4double rmax_voxelz = 5.5 * cm;
  G4double zz_voxels = 5. * cm;
  G4double startPhi_voxelz = 0.* degree;
  G4double spanningPhi_voxelz = 180. * degree; 
  G4int nslice = 10;

  G4Tubs* voxelz = new G4Tubs("voxel_z", rmin_voxelz, rmax_voxelz, zz_voxels/(2*nslice),startPhi_voxelz, spanningPhi_voxelz);
  G4LogicalVolume* voxelz_log = new G4LogicalVolume(voxelz, adipose_glandular, "voxelz_log",0 , 0, 0);  
  G4VPhysicalVolume* voxelz_phys = new G4PVReplica("LinearArray", voxelz_log, physInnerBreast,
						   kZAxis, nslice, zz_voxels/nslice);

  G4double voxel_height = (zz_voxels/nslice);

  G4int n_voxels = 10;

  G4double voxel_phi = spanningPhi_voxelz/n_voxels;

  G4cout <<"voxel_phi: " <<voxel_phi/degree << " deg"<< G4endl;
 
  G4Tubs* voxel = new G4Tubs("voxel", rmin_voxelz, rmax_voxelz, voxel_height/2.,
			     startPhi_voxelz,voxel_phi);
			     
  G4LogicalVolume* voxel_log = new G4LogicalVolume(voxel, adipose_glandular, "voxel_log", 0,0,0);

  G4VPhysicalVolume* voxel_phys = new G4PVReplica("RZPhiSlices",
 						    voxel_log,
 						    voxelz_phys,
						  kPhi, n_voxels, voxel_phi,- (voxel_phi/2.)); 
						  
 
 */
  G4VisAttributes* BreastVisAtt = new G4VisAttributes(G4Colour(1.0,0.41,0.71));
  BreastVisAtt -> SetForceSolid(false);
  breast_log -> SetVisAttributes(BreastVisAtt);

  G4VisAttributes* innerBreastVisAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.));
  innerBreastVisAtt -> SetForceSolid(false);
  innerBreast_log -> SetVisAttributes(innerBreastVisAtt);
  // voxelz_log -> SetVisAttributes(innerBreastVisAtt);
  // voxel_log -> SetVisAttributes(innerBreastVisAtt);

  G4cout << "Parameterised Breast created !!!!!!" << G4endl;
  
  // Testing Breast Volume
  G4double BreastVol = breast_log->GetSolid()->GetCubicVolume();
  G4cout << "Volume of Breast = " << BreastVol/cm3 << " cm^3" << G4endl;
  
  // Testing Breast Material
  G4String BreastMat = breast_log->GetMaterial()->GetName();
  G4cout << "Material of Breast = " << BreastMat << G4endl;
  
  G4String innerBreastMat = innerBreast_log->GetMaterial()->GetName();
  G4cout << "Material of the Inner Breast = " << BreastMat << G4endl;
  
  // Testing Density
  //G4double BreastDensity = logicBreast->GetMaterial()->GetDensity();
  //G4cout << "Density of Material = " << BreastDensity*cm3/g << " g/cm^3" << G4endl;

  // Testing Mass
  //G4double BreastMass = (BreastVol)*BreastDensity;
  //G4cout << "Mass of Breast = " << BreastMass/gram << " g" << G4endl;

 
  // Sensitive detector
 // Sensitive Body Part
  if (sensitivity==true)
  { 
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    
    if(!breastSD)
      {
       G4String sensitiveDetectorName = "Parameterised_breast";
       
       breastSD = new G4HumanPhantomSDBreast(sensitiveDetectorName);
   
       G4String ROGeometryName = "PhantomROGeometry";
       
        G4HumanPhantomROGeometry* phantomROGeometry = new G4HumanPhantomROGeometry(ROGeometryName);
        phantomROGeometry -> BuildROGeometry();
	breastSD -> SetROgeometry(phantomROGeometry);
        SDman -> AddNewDetector(breastSD);
        innerBreast_log -> SetSensitiveDetector(breastSD);
      }
  }
 return physInnerBreast;

}
