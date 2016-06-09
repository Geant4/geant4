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
// This is the *BASIC* version of Hadrontherapy, a Geant4-based application
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy
//
// Visit the Hadrontherapy web site (http://www.lns.infn.it/link/Hadrontherapy) to request 
// the *COMPLETE* version of this program, together with its documentation;
// Hadrontherapy (both basic and full version) are supported by the Italian INFN
// Institute in the framework of the MC-INFN Group
//

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4Colour.hh"
#include "G4UserLimits.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4NistManager.hh"

#include "HadrontherapyDetectorConstruction.hh"
#include "HadrontherapyDetectorROGeometry.hh"
#include "HadrontherapyDetectorMessenger.hh"
#include "HadrontherapyDetectorSD.hh"
#include "HadrontherapyMatrix.hh"
#include "HadrontherapyAnalysisManager.hh"

#include <cmath>

/////////////////////////////////////////////////////////////////////////////
HadrontherapyDetectorConstruction::HadrontherapyDetectorConstruction(G4VPhysicalVolume* physicalTreatmentRoom)
  : motherPhys(physicalTreatmentRoom), // pointer to WORLD volume 
    detectorSD(0), detectorROGeometry(0), matrix(0), 
    phantom(0), detector(0),
    phantomLogicalVolume(0), detectorLogicalVolume(0), 
    phantomPhysicalVolume(0), detectorPhysicalVolume(0),
    aRegion(0)
{
  HadrontherapyAnalysisManager::GetInstance();

  /* NOTE! that the HadrontherapyDetectorConstruction class
   * does NOT inherit from G4VUserDetectorConstruction G4 class
   * So the Construct() mandatory virtual method is inside another geometric class
   * like the passiveProtonBeamLIne, ...
   */

  // Messenger to change parameters of the phantom/detector geometry
  detectorMessenger = new HadrontherapyDetectorMessenger(this);

  // Default detector voxels size
  // 200 slabs along the beam direction (X)
  sizeOfVoxelAlongX = 200 *CLHEP::um; 
  sizeOfVoxelAlongY = 4 *CLHEP::cm; 
  sizeOfVoxelAlongZ = 4 *CLHEP::cm; 

  // Define here the material of the water phantom and of the detector
  SetPhantomMaterial("G4_WATER"); 
  // Construct geometry (messenger commands)
  SetDetectorSize(4.*CLHEP::cm, 4.*CLHEP::cm, 4.*CLHEP::cm);
  SetPhantomSize(40. *CLHEP::cm, 40. *CLHEP::cm, 40. *CLHEP::cm);
  SetPhantomPosition(G4ThreeVector(20. *CLHEP::cm, 0. *CLHEP::cm, 0. *CLHEP::cm));
  SetDetectorToPhantomPosition(G4ThreeVector(0. *CLHEP::cm, 18. *CLHEP::cm, 18. *CLHEP::cm));


  // Write virtual parameters to the real ones and check for consistency      
  UpdateGeometry();
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyDetectorConstruction::~HadrontherapyDetectorConstruction()
{ 
    delete detectorROGeometry;  
    delete matrix;  
    delete detectorMessenger;
}

/////////////////////////////////////////////////////////////////////////////
// ConstructPhantom() is the method that reconstuct a water box (called phantom 
// (or water phantom) in the usual Medical physicists slang). 
// A water phantom can be considered a good
// approximation of a an human body. 
void HadrontherapyDetectorConstruction::ConstructPhantom()
{
    // Definition of the solid volume of the Phantom
    phantom = new G4Box("Phantom", 
			phantomSizeX/2, 
			phantomSizeY/2, 
			phantomSizeZ/2);
    
// Definition of the logical volume of the Phantom
    phantomLogicalVolume = new G4LogicalVolume(phantom,	
					     phantomMaterial, 
					     "phantomLog", 0, 0, 0);
  
    // Definition of the physics volume of the Phantom
    phantomPhysicalVolume = new G4PVPlacement(0,
	                                    phantomPosition,
					    "phantomPhys",
					    phantomLogicalVolume,
					    motherPhys,
					    false,
					    0);

// Visualisation attributes of the phantom
    red = new G4VisAttributes(G4Colour(255/255., 0/255. ,0/255.));
    red -> SetVisibility(true);
    red -> SetForceSolid(true);
    red -> SetForceWireframe(true);
    phantomLogicalVolume -> SetVisAttributes(red);
}

/////////////////////////////////////////////////////////////////////////////
// ConstructDetector() it the method the reconstruct a detector region 
// inside the water phantom. It is a volume, located inside the water phantom
// and with two coincident faces:
//
//           **************************
//           *   water phantom        *
//           *                        *
//           *                        *
//           *---------------         *
//  Beam     *              -         *
//  ----->   * detector     -         *
//           *              -         *
//           *---------------         *
//           *                        *
//           *                        *
//           *                        *
//           **************************
//
// The detector is the volume that can be dived in slices or voxelized
// and in it we can collect a number of usefull information:
// dose distribution, fluence distribution, LET and so on
void HadrontherapyDetectorConstruction::ConstructDetector()
{

    // Definition of the solid volume of the Detector
    detector = new G4Box("Detector", 
			 detectorSizeX/2, 
			 detectorSizeY/2, 
			 detectorSizeZ/2);
    
    // Definition of the logic volume of the Phantom
    detectorLogicalVolume = new G4LogicalVolume(detector,
						detectorMaterial,
						"DetectorLog",
						0,0,0);
// Definition of the physical volume of the Phantom 
    detectorPhysicalVolume = new G4PVPlacement(0, 
					       detectorPosition, // Setted by displacement 
					       "DetectorPhys", 
					       detectorLogicalVolume, 
					       phantomPhysicalVolume, 
					       false,0);

// Visualisation attributes of the detector 
    skyBlue = new G4VisAttributes( G4Colour(135/255. , 206/255. ,  235/255. ));
    skyBlue -> SetVisibility(true);
    skyBlue -> SetForceSolid(true);
    //skyBlue -> SetForceWireframe(true);
    detectorLogicalVolume -> SetVisAttributes(skyBlue);

  // **************
  // Cut per Region
  // **************
  
  // A smaller cut is fixed in the phantom to calculate the energy deposit with the
  // required accuracy 
    if (!aRegion)
    {
	aRegion = new G4Region("DetectorLog");
	detectorLogicalVolume -> SetRegion(aRegion);
	aRegion -> AddRootLogicalVolume(detectorLogicalVolume);
    }

}

/////////////////////////////////////////////////////////////////////////////

void  HadrontherapyDetectorConstruction::ConstructSensitiveDetector(G4ThreeVector detectorToWorldPosition)
{  
    // Install new Sensitive Detector and ROGeometry 
    delete detectorROGeometry; // this should be safe in C++ also if we have a NULL pointer
    //if (detectorSD) detectorSD->PrintAll();
    //delete detectorSD;
    // Sensitive Detector and ReadOut geometry definition
    G4SDManager* sensitiveDetectorManager = G4SDManager::GetSDMpointer();

    static G4String sensitiveDetectorName = "Detector"; 
    if (!detectorSD)
	{
	    // The sensitive detector is instantiated
	    detectorSD = new HadrontherapyDetectorSD(sensitiveDetectorName);
	}
    // The Read Out Geometry is instantiated
    static G4String ROGeometryName = "DetectorROGeometry";
    detectorROGeometry = new HadrontherapyDetectorROGeometry(ROGeometryName,
							    detectorToWorldPosition,
							    detectorSizeX/2,
							    detectorSizeY/2,
							    detectorSizeZ/2,
							    numberOfVoxelsAlongX,
							    numberOfVoxelsAlongY,
							    numberOfVoxelsAlongZ);

    G4cout << "Instantiating new Read Out Geometry \"" << ROGeometryName << "\""<< G4endl;
    // This will invoke Build() HadrontherapyDetectorROGeometry virtual method 
    detectorROGeometry -> BuildROGeometry();
    // Attach ROGeometry to SDetector
    detectorSD -> SetROgeometry(detectorROGeometry);
    //sensitiveDetectorManager -> Activate(sensitiveDetectorName, true);
    if (!sensitiveDetectorManager -> FindSensitiveDetector(sensitiveDetectorName, false))
	{
	    G4cout << "Registering new DetectorSD \"" << sensitiveDetectorName << "\""<< G4endl;
	    // Register user SD
	    sensitiveDetectorManager -> AddNewDetector(detectorSD);
	    // Attach SD to detector logical volume
	    detectorLogicalVolume -> SetSensitiveDetector(detectorSD);
	}
}
void  HadrontherapyDetectorConstruction::ParametersCheck()
{
    // Check phantom/detector sizes & relative position
    if (!IsInside(detectorSizeX, 
		detectorSizeY, 
		detectorSizeZ,
		phantomSizeX,
		phantomSizeY,
		phantomSizeZ,
		detectorToPhantomPosition
		))
      G4Exception("HadrontherapyDetectorConstruction::ParametersCheck()", "Hadrontherapy0001", FatalException, "Error: Detector is not fully inside Phantom!");

    // Check Detector sizes respect to the voxel ones

    if ( detectorSizeX < sizeOfVoxelAlongX) {
      G4Exception("HadrontherapyDetectorConstruction::ParametersCheck()", "Hadrontherapy0002", FatalException, "Error:  Detector X size must be bigger or equal than that of Voxel X!");
    }
    if ( detectorSizeY < sizeOfVoxelAlongY) {
      G4Exception(" HadrontherapyDetectorConstruction::ParametersCheck()", "Hadrontherapy0003", FatalException, "Error:  Detector Y size must be bigger or equal than that of Voxel Y!");
    }
    if ( detectorSizeZ < sizeOfVoxelAlongZ) {
      G4Exception(" HadrontherapyDetectorConstruction::ParametersCheck()", "Hadrontherapy0004", FatalException, "Error:  Detector Z size must be bigger or equal than that of Voxel Z!");
    }

}
/////////////////
// MESSENGERS //
////////////////

G4bool HadrontherapyDetectorConstruction::SetPhantomMaterial(G4String material)
{

    if (G4Material* pMat = G4NistManager::Instance()->FindOrBuildMaterial(material, false) )
    {
	phantomMaterial  = pMat;
	detectorMaterial = pMat;
	if (detectorLogicalVolume && phantomLogicalVolume) 
	{
	    detectorLogicalVolume -> SetMaterial(pMat); 
	    phantomLogicalVolume ->  SetMaterial(pMat);

	    G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
	    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
	    G4cout << "The material of Phantom/Detector has been changed to " << material << G4endl;
	}
    }
    else
    {
	G4cout << "WARNING: material \"" << material << "\" doesn't exist in NIST elements/materials"
	    " table [located in $G4INSTALL/source/materials/src/G4NistMaterialBuilder.cc]" << G4endl; 
	G4cout << "Use command \"/parameter/nist\" to see full materials list!" << G4endl; 
	return false;
    }

    return true;
}
/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorConstruction::SetPhantomSize(G4double sizeX, G4double sizeY, G4double sizeZ)
{
    if (sizeX > 0.) phantomSizeX = sizeX;
    if (sizeY > 0.) phantomSizeY = sizeY;
    if (sizeZ > 0.) phantomSizeZ = sizeZ;
}
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorConstruction::SetDetectorSize(G4double sizeX, G4double sizeY, G4double sizeZ)
{
    if (sizeX > 0.) {detectorSizeX = sizeX;}
    if (sizeY > 0.) {detectorSizeY = sizeY;}
    if (sizeZ > 0.) {detectorSizeZ = sizeZ;}
    SetVoxelSize(sizeOfVoxelAlongX, sizeOfVoxelAlongY, sizeOfVoxelAlongZ);
}
/////////////////////////////////////////////////////////////////////////////

void HadrontherapyDetectorConstruction::SetVoxelSize(G4double sizeX, G4double sizeY, G4double sizeZ)
{
    if (sizeX > 0.) {sizeOfVoxelAlongX = sizeX;}
    if (sizeY > 0.) {sizeOfVoxelAlongY = sizeY;}
    if (sizeZ > 0.) {sizeOfVoxelAlongZ = sizeZ;}
}
void HadrontherapyDetectorConstruction::SetPhantomPosition(G4ThreeVector pos)
{
    phantomPosition = pos;
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorConstruction::SetDetectorToPhantomPosition(G4ThreeVector displ)
{
    detectorToPhantomPosition = displ;
}
/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorConstruction::UpdateGeometry()
{
    /* 
     * Check parameters consistency
     */
    ParametersCheck();

    G4GeometryManager::GetInstance() -> OpenGeometry();
    if (phantom)
    {
	phantom -> SetXHalfLength(phantomSizeX/2);
	phantom -> SetYHalfLength(phantomSizeY/2);
	phantom -> SetZHalfLength(phantomSizeZ/2);
	phantomPhysicalVolume -> SetTranslation(phantomPosition);
    }
    else   ConstructPhantom();

    // Get the center of the detector 
    SetDetectorPosition();
    if (detector)
    {
	detector -> SetXHalfLength(detectorSizeX/2);
	detector -> SetYHalfLength(detectorSizeY/2);
	detector -> SetZHalfLength(detectorSizeZ/2);
	detectorPhysicalVolume -> SetTranslation(detectorPosition);
    }
    else    ConstructDetector();

    // Round to nearest integer number of voxel 
    numberOfVoxelsAlongX = G4lrint(detectorSizeX / sizeOfVoxelAlongX);
    sizeOfVoxelAlongX = ( detectorSizeX / numberOfVoxelsAlongX );

    numberOfVoxelsAlongY = G4lrint(detectorSizeY / sizeOfVoxelAlongY);
    sizeOfVoxelAlongY = ( detectorSizeY / numberOfVoxelsAlongY );

    numberOfVoxelsAlongZ = G4lrint(detectorSizeZ / sizeOfVoxelAlongZ);
    sizeOfVoxelAlongZ = ( detectorSizeZ / numberOfVoxelsAlongZ );

    ConstructSensitiveDetector(GetDetectorToWorldPosition());

    volumeOfVoxel = sizeOfVoxelAlongX * sizeOfVoxelAlongY * sizeOfVoxelAlongZ;
    massOfVoxel = detectorMaterial -> GetDensity() * volumeOfVoxel;
    //  This will clear the existing matrix (together with all data inside it)! 
    matrix = HadrontherapyMatrix::GetInstance(numberOfVoxelsAlongX, 
	    numberOfVoxelsAlongY,
	    numberOfVoxelsAlongZ,
	    massOfVoxel);

    // Initialize analysis
#ifdef G4ANALYSIS_USE_ROOT
        HadrontherapyAnalysisManager* analysis = HadrontherapyAnalysisManager::GetInstance();
    analysis -> flush();     // Finalize the root file 
    analysis -> book();
#endif
    // Inform the kernel about the new geometry
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4RunManager::GetRunManager() -> PhysicsHasBeenModified();

    PrintParameters();
}

void HadrontherapyDetectorConstruction::PrintParameters()
{

    G4cout << "The (X,Y,Z) dimensions of the phantom are : (" << 
	G4BestUnit( phantom -> GetXHalfLength()*2., "Length") << ',' << 
	G4BestUnit( phantom -> GetYHalfLength()*2., "Length") << ',' << 
	G4BestUnit( phantom -> GetZHalfLength()*2., "Length") << ')' << G4endl; 
    
    G4cout << "The (X,Y,Z) dimensions of the detector are : (" << 
	G4BestUnit( detector -> GetXHalfLength()*2., "Length") << ',' << 
	G4BestUnit( detector -> GetYHalfLength()*2., "Length") << ',' << 
	G4BestUnit( detector -> GetZHalfLength()*2., "Length") << ')' << G4endl; 

    G4cout << "Displacement between Phantom and World is: "; 
    G4cout << "DX= "<< G4BestUnit(phantomPosition.getX(),"Length") << 
	"DY= "<< G4BestUnit(phantomPosition.getY(),"Length") << 
	"DZ= "<< G4BestUnit(phantomPosition.getZ(),"Length") << G4endl;

    G4cout << "The (X,Y,Z) sizes of the Voxels are: (" << 
	G4BestUnit(sizeOfVoxelAlongX, "Length")  << ',' << 
	G4BestUnit(sizeOfVoxelAlongY, "Length")  << ',' << 
	G4BestUnit(sizeOfVoxelAlongZ, "Length") << ')' << G4endl;

    G4cout << "The number of Voxels along (X,Y,Z) is: (" << 
	numberOfVoxelsAlongX  << ',' <<
	numberOfVoxelsAlongY  <<','  <<
	numberOfVoxelsAlongZ  << ')' << G4endl;

}
