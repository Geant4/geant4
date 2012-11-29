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
// This is the *BASIC* version of IORT, a Geant4-based application
//
// Main Authors: G.Russo(a,b), C.Casarino*(c), G.C. Candiano(c), G.A.P. Cirrone(d), F.Romano(d)
// Contributor Authors: S.Guatelli(e)
// Past Authors: G.Arnetta(c), S.E.Mazzaglia(d)
//    
//   (a) Fondazione Istituto San Raffaele G.Giglio, Cefalù, Italy
//   (b) IBFM-CNR , Segrate (Milano), Italy
//   (c) LATO (Laboratorio di Tecnologie Oncologiche), Cefalù, Italy
//   (d) Laboratori Nazionali del Sud of the INFN, Catania, Italy
//   (e) University of Wallongong, Australia
//
//   *Corresponding author, email to carlo.casarino@polooncologicocefalu.it
//////////////////////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
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
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4Colour.hh"
#include "G4UserLimits.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4NistManager.hh"
#include "IORTDetectorConstruction.hh"
#include "IORTDetectorROGeometry.hh"
#include "IORTDetectorMessenger.hh"
#include "IORTDetectorSD.hh"
#include "IORTMatrix.hh"
#include "IORTAnalysisManager.hh"
#include "G4Tubs.hh"

/////////////////////////////////////////////////////////////////////////////
IORTDetectorConstruction::IORTDetectorConstruction(G4VPhysicalVolume* physicalTreatmentRoom)
  : motherPhys(physicalTreatmentRoom), // pointer to WORLD volume 
    detectorSD(0), detectorROGeometry(0), matrix(0), 
    phantom(0), detector(0),
    phantomLogicalVolume(0), detectorLogicalVolume(0), 
    phantomPhysicalVolume(0), detectorPhysicalVolume(0),
    aRegion(0),

    solidDiscoIORT0(0),
    logicDiscoIORT0(0),
    physiDiscoIORT0(0),

    solidDiscoIORT(0),
    logicDiscoIORT(0),
    physiDiscoIORT(0),

    solidDiscoIORT1(0),
    logicDiscoIORT1(0),
    physiDiscoIORT1(0)

{
  IORTAnalysisManager::GetInstance();

  /* NOTE! that the IORTDetectorConstruction class
   * does NOT inherit from G4VUserDetectorConstruction G4 class
   * So the Construct() mandatory virtual method is inside another geometric class
   * like the collimatorXXBeamLIne, ...
   */

  // Messenger to change parameters of the phantom/detector geometry
  detectorMessenger = new IORTDetectorMessenger(this);

  // Default detector voxels size
  // 200 slabs along the beam direction (X)
  sizeOfVoxelAlongX = 0.5 *CLHEP::mm; // 
  sizeOfVoxelAlongY = 0.5 *CLHEP::mm; //  
  sizeOfVoxelAlongZ = 0.5 *CLHEP::mm; // 

  // Define here the material of the water phantom and of the detector
  SetPhantomMaterial("G4_WATER"); 
  // Construct geometry (messenger commands)
  SetDetectorSize(7.*CLHEP::cm, 15.*CLHEP::cm, 15.*CLHEP::cm);    
  SetPhantomSize(20. *CLHEP::cm, 20. *CLHEP::cm, 20. *CLHEP::cm);   
  SetPhantomPosition(G4ThreeVector(4.5 *CLHEP::cm, 0. *CLHEP::cm, 0. *CLHEP::cm)); 
  SetDetectorToPhantomPosition(G4ThreeVector(0. *CLHEP::cm, 2.5 *CLHEP::cm, 2.5 *CLHEP::cm));  

  // Default protection disc geometry and materials
  SetOuterRadiusDiscoIORT (40. *CLHEP::mm);  
  SetinnerRadiusDiscoIORT (0.*CLHEP::mm);   
  SetheightDiscoIORT (2.0*CLHEP::mm);        
  SetDiscoXPositionIORT (-11.0*CLHEP::mm);
  SetDiscoYPositionIORT (0.0*CLHEP::mm);
  SetDiscoZPositionIORT (0.0*CLHEP::mm);
  SetDiscoMaterialIORT("G4_WATER");   

  SetOuterRadiusDiscoIORT1 (40. *CLHEP::mm);  
  SetinnerRadiusDiscoIORT1 (0.*CLHEP::mm);   
  SetheightDiscoIORT1 (1.0*CLHEP::mm);        
  SetDiscoXPositionIORT1 (-8.0*CLHEP::mm);
  SetDiscoMaterialIORT1("G4_WATER");

  SetAngleDiscoIORT0 (90.0 *CLHEP::deg);    

  // Write virtual parameters to the real ones and check for consistency      
  UpdateGeometry();
}

/////////////////////////////////////////////////////////////////////////////
IORTDetectorConstruction::~IORTDetectorConstruction()
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
void IORTDetectorConstruction::ConstructPhantom()
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
    //red -> SetForceSolid(true);
    //red -> SetForceWireframe(true);
    phantomLogicalVolume -> SetVisAttributes(red); 
    //phantomLogicalVolume -> SetVisAttributes(G4VisAttributes::Invisible);
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
void IORTDetectorConstruction::ConstructDetector()
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
    //skyBlue = new G4VisAttributes( G4Colour(135/255. , 206/255. ,  235/255. ));
    G4VisAttributes * skyBlue1 = new G4VisAttributes( G4Colour(135/255. , 206/255. ,  235/255. ));
    //skyBlue1 -> SetForceWireframe(true);
    //skyBlue1 -> SetForceSolid(true);
    //skyBlue -> SetVisibility(true);
    //skyBlue -> SetForceSolid(true);
    //skyBlue -> SetForceWireframe(false);
    //detectorLogicalVolume -> SetVisAttributes(skyBlue);
    detectorLogicalVolume -> SetVisAttributes(skyBlue1);
    
   // detectorLogicalVolume -> SetVisAttributes(G4VisAttributes::Invisible);


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

void IORTDetectorConstruction::ConstructDisc()
{
// ---------------------------------------------------------------//
  //                    6.0 mm Protection Discs Volume          //
    // ---------------------------------------------------------------//
  const G4double startAngleDiscoIORT0 = 0.*CLHEP::deg;
  const G4double spanningAngleDiscoIORT0 = 360.*CLHEP::deg;

  //G4double phi0 = 180. *CLHEP::deg; // messenger    

  // Matrix definition for a rotation (deg).       
  G4RotationMatrix rm0;               
  rm0.rotateY(AngleDiscoIORT0);

  
  solidDiscoIORT0 = new G4Tubs("DiscoIORT0", innerRadiusDiscoIORT, 
				    OuterRadiusDiscoIORT,
				    (heightDiscoIORT + heightDiscoIORT1), 
				    startAngleDiscoIORT0, 
				    spanningAngleDiscoIORT0);

  G4LogicalVolume* logDiscoIORT0 = new G4LogicalVolume(solidDiscoIORT0, 
							      detectorMaterial, "DiscoIORT0Log", 0, 0, 0);

  physiDiscoIORT0 = new G4PVPlacement(G4Transform3D(rm0, G4ThreeVector((DiscoXPositionIORT + heightDiscoIORT1),DiscoYPositionIORT,DiscoZPositionIORT)),
					   "DiscoIORT0Phys", logDiscoIORT0, detectorPhysicalVolume, false, 0); 

  white = new G4VisAttributes( G4Colour());
  white -> SetVisibility(true);
  // white -> SetForceSolid(true);
  logDiscoIORT0 -> SetVisAttributes(white);

// ---------------------------------------------------------------//
  //                    4.0 mm Aluminium Protection Disc          //
    // ---------------------------------------------------------------//
  //G4bool isotopes = false;
 // G4Material* leadNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb", isotopes);
 // DiscoMaterialIORT = leadNist;   // messenger
  gray = new G4VisAttributes( G4Colour(0.5, 0.5, 0.5 ));
  gray-> SetVisibility(true);
  //gray -> SetForceWireframe(true);
  //gray-> SetForceSolid(true);

  gray1 = new G4VisAttributes( G4Colour(0.7, 0.7, 0.7 ));
  gray1-> SetVisibility(true);
  //gray1 -> SetForceWireframe(true);
  //gray1-> SetForceSolid(true);
 // const G4double OuterRadiusDiscoIORT = 35. *CLHEP::mm;  // messenger
 // const G4double innerRadiusDiscoIORT = 0.*CLHEP::mm;   // messenger
 //  const G4double heightDiscoIORT = 3.0*CLHEP::mm;        // messenger
  const G4double startAngleDiscoIORT = 0.*CLHEP::deg;
  const G4double spanningAngleDiscoIORT = 360.*CLHEP::deg;
 // const G4double DiscoXPositionIORT = -14.0*CLHEP::mm; // messenger 

//G4Material* DiscoMaterialIORT = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLEXIGLASS", isotopes);// messenger 
  
   
  G4double phi = 0. *CLHEP::deg;     

  // Matrix definition for a 90 deg rotation. Also used for other volumes       
  G4RotationMatrix rm;               
  rm.rotateY(phi);

  
  solidDiscoIORT = new G4Tubs("DiscoIORT", innerRadiusDiscoIORT, 
				    OuterRadiusDiscoIORT,
				    heightDiscoIORT, 
				    startAngleDiscoIORT, 
				    spanningAngleDiscoIORT);

  G4LogicalVolume* logDiscoIORT = new G4LogicalVolume(solidDiscoIORT, 
							      DiscoMaterialIORT, "DiscoIORTLog", 0, 0, 0);

  physiDiscoIORT = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector(0.,0.,(- heightDiscoIORT1))),
					   "DiscoIORTPhys", logDiscoIORT, physiDiscoIORT0, false, 0); 
  
  logDiscoIORT -> SetVisAttributes(gray1);
  
  
      // ---------------------------------------------------------------//
      //             2.0 mm Lead Protection Disc                          //
      // ---------------------------------------------------------------//
  
 // const G4double OuterRadiusDiscoIORT1 = 35. *CLHEP::mm;
 // const G4double innerRadiusDiscoIORT1 = 0.*CLHEP::mm;
 // const G4double heightDiscoIORT1 = 0.5*CLHEP::mm;
  const G4double startAngleDiscoIORT1 = 0.*CLHEP::deg;
  const G4double spanningAngleDiscoIORT1 = 360.*CLHEP::deg;
//  const G4double DiscoXPositionIORT1 = -10.5*CLHEP::mm; messenger
//  G4Material* DiscoMaterialIORT1 = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu", isotopes);// messenger  
  
       
  
  solidDiscoIORT1 = new G4Tubs("DiscoIORT1", innerRadiusDiscoIORT1, 
				    OuterRadiusDiscoIORT1,
				    heightDiscoIORT1, 
				    startAngleDiscoIORT1, 
				    spanningAngleDiscoIORT1);

  G4LogicalVolume* logDiscoIORT1 = new G4LogicalVolume(solidDiscoIORT1, 
							      DiscoMaterialIORT1, "DiscoIORTLog1", 0, 0, 0);

  physiDiscoIORT1 = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector(0.,0.,heightDiscoIORT)),
					   "DiscoIORTPhys1", logDiscoIORT1, physiDiscoIORT0, false, 0); 
  white = new G4VisAttributes( G4Colour());
  white -> SetVisibility(true);
  white -> SetForceSolid(true);
  logDiscoIORT1 -> SetVisAttributes(gray);

}
/////////////////////////////////////////////////////////////////////////////

void  IORTDetectorConstruction::ConstructSensitiveDetector(G4ThreeVector detectorToWorldPosition)
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
	    detectorSD = new IORTDetectorSD(sensitiveDetectorName);
	}
    // The Read Out Geometry is instantiated
    static G4String ROGeometryName = "DetectorROGeometry";
    detectorROGeometry = new IORTDetectorROGeometry(ROGeometryName,
							    detectorToWorldPosition,
							    detectorSizeX/2,  // controllare che sia necessario /2
							    detectorSizeY/2,  // CONFERMATO!!! ci vuole!!!
							    detectorSizeZ/2,
							    numberOfVoxelsAlongX,
							    numberOfVoxelsAlongY,
							    numberOfVoxelsAlongZ);

    G4cout << "Instantiating new Read Out Geometry \"" << ROGeometryName << "\""<< G4endl;
    // This will invoke Build() IORTDetectorROGeometry virtual method 
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
void  IORTDetectorConstruction::ParametersCheck()
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
      G4Exception("IORTDetectorConstruction::ParametersCheck()", "IORT0001", FatalException, "Error: Detector is not fully inside Phantom!");

    // Check Detector sizes respect to the voxel ones

    if ( detectorSizeX < sizeOfVoxelAlongX) {
      G4Exception("IORTDetectorConstruction::ParametersCheck()", "IORT0002", FatalException, "Error: Detector X size must be bigger or equal than that of Voxel X");
    }
    if ( detectorSizeY < sizeOfVoxelAlongY) {
      G4Exception("IORTDetectorConstruction::ParametersCheck()", "IORT0003", FatalException, "Error: Detector X size must be bigger or equal than that of Voxel Y");	
    }
    if ( detectorSizeZ < sizeOfVoxelAlongZ) {
      G4Exception("IORTDetectorConstruction::ParametersCheck()", "IORT0004", FatalException, "Error: Detector X size must be bigger or equal than that of Voxel Z");
    }

}
/////////////////
// MESSENGERS //
////////////////

G4bool IORTDetectorConstruction::SetPhantomMaterial(G4String material)
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

G4bool IORTDetectorConstruction::SetDiscoMaterialIORT(G4String material)
{

    if (G4Material* dMat = G4NistManager::Instance()->FindOrBuildMaterial(material, false) )
    {
	DiscoMaterialIORT  = dMat;
	
	if (logicDiscoIORT) 
	{
	    logicDiscoIORT -> SetMaterial(dMat); 
	    
	    G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
	    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
	    G4cout << "The material of Protection disc 1 has been changed to " << material << G4endl;
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

G4bool IORTDetectorConstruction::SetDiscoMaterialIORT1(G4String material)
{

    if (G4Material* d1Mat = G4NistManager::Instance()->FindOrBuildMaterial(material, false) )
    {
	DiscoMaterialIORT1  = d1Mat;
	
	if (logicDiscoIORT1) 
	{
	    logicDiscoIORT1 -> SetMaterial(d1Mat); 
	    
	    G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
	    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
	    G4cout << "The material of Protection disc 2 has been changed to " << material << G4endl;
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
void IORTDetectorConstruction::SetPhantomSize(G4double sizeX, G4double sizeY, G4double sizeZ)
{
    if (sizeX > 0.) phantomSizeX = sizeX;
    if (sizeY > 0.) phantomSizeY = sizeY;
    if (sizeZ > 0.) phantomSizeZ = sizeZ;
}
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
void IORTDetectorConstruction::SetDetectorSize(G4double sizeX, G4double sizeY, G4double sizeZ)
{
    if (sizeX > 0.) {detectorSizeX = sizeX;}
    if (sizeY > 0.) {detectorSizeY = sizeY;}
    if (sizeZ > 0.) {detectorSizeZ = sizeZ;}
    SetVoxelSize(sizeOfVoxelAlongX, sizeOfVoxelAlongY, sizeOfVoxelAlongZ);
}
/////////////////////////////////////////////////////////////////////////////

void IORTDetectorConstruction::SetVoxelSize(G4double sizeX, G4double sizeY, G4double sizeZ)
{
    if (sizeX > 0.) {sizeOfVoxelAlongX = sizeX;}
    if (sizeY > 0.) {sizeOfVoxelAlongY = sizeY;}
    if (sizeZ > 0.) {sizeOfVoxelAlongZ = sizeZ;}
}
void IORTDetectorConstruction::SetPhantomPosition(G4ThreeVector pos)
{
    phantomPosition = pos;
}

/////////////////////////////////////////////////////////////////////////////
void IORTDetectorConstruction::SetDetectorToPhantomPosition(G4ThreeVector displ)
{
    detectorToPhantomPosition = displ;
}
/////////////////////////////////protection disc///////////////////////////

void IORTDetectorConstruction::SetOuterRadiusDiscoIORT(G4double outerr)
{
    if (outerr > 0.) {OuterRadiusDiscoIORT = outerr;}
    
}

void IORTDetectorConstruction::SetinnerRadiusDiscoIORT(G4double innerr)
{
    if (innerr > 0.) {innerRadiusDiscoIORT = innerr;}
    
}

void IORTDetectorConstruction::SetheightDiscoIORT(G4double height)
{
    if (height > 0.) {heightDiscoIORT = height;}
    
}

void IORTDetectorConstruction::SetDiscoXPositionIORT(G4double xpos)
{
    
    DiscoXPositionIORT = xpos;
        
}

void IORTDetectorConstruction::SetDiscoYPositionIORT(G4double ypos)
{
    
    DiscoYPositionIORT = ypos;
        
}

void IORTDetectorConstruction::SetDiscoZPositionIORT(G4double zpos)
{
    
    DiscoZPositionIORT = zpos;
        
}

void IORTDetectorConstruction::SetOuterRadiusDiscoIORT1(G4double outerr)
{
    if (outerr > 0.) {OuterRadiusDiscoIORT1 = outerr;}
    
}

void IORTDetectorConstruction::SetinnerRadiusDiscoIORT1(G4double innerr)
{
    if (innerr > 0.) {innerRadiusDiscoIORT1 = innerr;}
    
}

void IORTDetectorConstruction::SetheightDiscoIORT1(G4double height)
{
    if (height > 0.) {heightDiscoIORT1 = height;}
    
}

void IORTDetectorConstruction::SetDiscoXPositionIORT1(G4double xpos)
{
    
    DiscoXPositionIORT1 = xpos;
}

void IORTDetectorConstruction::SetAngleDiscoIORT0(G4double phi0)
{
    
    AngleDiscoIORT0 = phi0;
}

/////////////////////////////////protection disc///////////end///////////


////////////////////////////////////////////////////////////////////////////////
void IORTDetectorConstruction::UpdateGeometry()
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

    // update disc function
    delete solidDiscoIORT0;
    delete logicDiscoIORT0;
    delete physiDiscoIORT0;
    delete solidDiscoIORT;
    delete logicDiscoIORT;
    delete physiDiscoIORT;
    delete solidDiscoIORT1;
    delete logicDiscoIORT1;
    delete physiDiscoIORT1;
    ConstructDisc();
    

    // Round to nearest integer number of voxel 
    numberOfVoxelsAlongX = G4lrint(detectorSizeX / sizeOfVoxelAlongX);
    sizeOfVoxelAlongX = ( detectorSizeX / numberOfVoxelsAlongX );

    numberOfVoxelsAlongY = G4lrint(detectorSizeY / sizeOfVoxelAlongY);
    sizeOfVoxelAlongY = ( detectorSizeY / numberOfVoxelsAlongY );

    numberOfVoxelsAlongZ = G4lrint(detectorSizeZ / sizeOfVoxelAlongZ);
    sizeOfVoxelAlongZ = ( detectorSizeZ / numberOfVoxelsAlongZ );

    //G4cout << "*************** DetectorToWorldPosition " << GetDetectorToWorldPosition()/cm << "\n";
    ConstructSensitiveDetector(GetDetectorToWorldPosition());

    volumeOfVoxel = sizeOfVoxelAlongX * sizeOfVoxelAlongY * sizeOfVoxelAlongZ;
    massOfVoxel = detectorMaterial -> GetDensity() * volumeOfVoxel;
    //  This will clear the existing matrix (together with all data inside it)! 
    matrix = IORTMatrix::GetInstance(numberOfVoxelsAlongX, 
	    numberOfVoxelsAlongY,
	    numberOfVoxelsAlongZ,
	    massOfVoxel);

    // Initialize analysis
/*
    IORTAnalysisManager* analysis = IORTAnalysisManager::GetInstance();
#ifdef G4ANALYSIS_USE_ROOT
    analysis -> flush();     // Finalize the root file 
    analysis -> book();
#endif
*/
    // Inform the kernel about the new geometry
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4RunManager::GetRunManager() -> PhysicsHasBeenModified();

    PrintParameters();
}

void IORTDetectorConstruction::DeleteDisc()
{
    delete solidDiscoIORT0;
    delete logicDiscoIORT0;
    delete physiDiscoIORT0;
    delete solidDiscoIORT;
    delete logicDiscoIORT;
    delete physiDiscoIORT;
    delete solidDiscoIORT1;
    delete logicDiscoIORT1;
    delete physiDiscoIORT1;
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4RunManager::GetRunManager() -> PhysicsHasBeenModified();
}


void IORTDetectorConstruction::PrintParameters()
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
