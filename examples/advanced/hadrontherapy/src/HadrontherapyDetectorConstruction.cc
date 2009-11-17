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
// HadrontherapyDetectorConstruction.cc
//
// See more at: http://workgroup.lngs.infn.it/geant4lns/
//
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the INFN, Catania, Italy
// (b) INFN Section of Genova
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------
#include "G4SDManager.hh"
#include "G4RunManager.hh"
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
#include "HadrontherapyDetectorROGeometry.hh"
#include "HadrontherapyDetectorMessenger.hh"
#include "HadrontherapyDetectorSD.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "HadrontherapyMatrix.hh"
/////////////////////////////////////////////////////////////////////////////
HadrontherapyDetectorConstruction::HadrontherapyDetectorConstruction(G4VPhysicalVolume* physicalTreatmentRoom)
  : motherPhys(physicalTreatmentRoom),
    detectorSD(0), detectorROGeometry(0), matrix(0), 
    phantomPhysicalVolume(0), 
    detectorLogicalVolume(0), detectorPhysicalVolume(0),
    phantomSizeX(20.*cm), phantomSizeY(20.*cm), phantomSizeZ(20.*cm), // Default half dimensions 
    detectorSizeX(2.*cm), detectorSizeY(2.*cm), detectorSizeZ(2.*cm),
    phantomPosition(20.*cm, 0.*cm, 0.*cm),
    detectorToPhantomPosition(0.*cm,18.*cm,18.*cm)// Default displacement of the detector respect to the phantom
{
  // NOTE! that the HadrontherapyDetectorConstruction class
  // does NOT inherit from G4VUserDetectorConstruction G4 class
  // So the Construct() mandatory virtual method is inside another geometric class
  // (like the passiveProtonBeamLIne, ...)

  // Messenger to change parameters of the phantom/detector geometry
  detectorMessenger = new HadrontherapyDetectorMessenger(this);

  // Default detector voxels size
  // 200 slabs along the beam direction (X)
  sizeOfVoxelAlongX = 200 *um;  
  sizeOfVoxelAlongY = 2 * detectorSizeY;
  sizeOfVoxelAlongZ = 2 * detectorSizeZ;

  // Calculate (and eventually set) detector position by displacement, phantom size and detector size
  SetDetectorPosition();

  // Build phantom and associated detector 
  ConstructPhantom();
  ConstructDetector();

  // Set number of the detector voxels along X Y and Z directions.  
  // This will construct also the sensitive detector, the ROGeometry 
  // and the matrix where the energy deposited is collected!
  SetNumberOfVoxelBySize(sizeOfVoxelAlongX, sizeOfVoxelAlongY, sizeOfVoxelAlongZ);
}

/////////////////////////////////////////////////////////////////////////////
HadrontherapyDetectorConstruction::~HadrontherapyDetectorConstruction()
{ 
    if (detectorROGeometry) delete detectorROGeometry;  
    if (matrix) delete matrix;  
    delete detectorMessenger;
}

void HadrontherapyDetectorConstruction::ConstructPhantom()
{
  //----------------------------------------
  // Phantom:
  // A box used to approximate tissues
  //----------------------------------------

    G4bool isotopes =  false; 
    G4Material* waterNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER", isotopes);
    phantom = new G4Box("Phantom",phantomSizeX, phantomSizeY, phantomSizeZ);
    phantomLogicalVolume = new G4LogicalVolume(phantom,	
					     waterNist, 
					     "phantomLog", 0, 0, 0);
  
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
//red -> SetForceWireframe(true);
    phantomLogicalVolume -> SetVisAttributes(red);
}

/////////////////////////////////////////////////////////////////////////////
void HadrontherapyDetectorConstruction::ConstructDetector()
{
  //-----------
  // Detector
  //-----------
    G4bool isotopes =  false; 
    G4Material* waterNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER", isotopes);
    detector = new G4Box("Detector",detectorSizeX,detectorSizeY,detectorSizeZ);
    detectorLogicalVolume = new G4LogicalVolume(detector,
						waterNist,
						"DetectorLog",
						0,0,0);
// Detector is attached by default to the phantom face directly exposed to the beam 
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
  
}
/////////////////////////////////////////////////////////////////////////////
void  HadrontherapyDetectorConstruction::ConstructSensitiveDetector(G4ThreeVector detectorToWorldPosition)
{  
    // Install new Sensitive Detector and ROGeometry 
    delete detectorROGeometry; // this should be safe in C++ also if we have a NULL pointer

    // Sensitive Detector and ReadOut geometry definition
    G4SDManager* sensitiveDetectorManager = G4SDManager::GetSDMpointer();

    G4String sensitiveDetectorName = "Detector"; 
    if (!detectorSD)
	{
	    // The sensitive detector is instantiated
	    detectorSD = new HadrontherapyDetectorSD(sensitiveDetectorName);
	}
    // The Read Out Geometry is instantiated
    G4String ROGeometryName = "DetectorROGeometry";
    detectorROGeometry = new HadrontherapyDetectorROGeometry(ROGeometryName,
							    detectorToWorldPosition,
							    detectorSizeX,
							    detectorSizeY,
							    detectorSizeZ,
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

/////////////////
// MESSENGERS //
////////////////
G4bool HadrontherapyDetectorConstruction::SetNumberOfVoxelBySize(G4double sizeX, G4double sizeY, G4double sizeZ)
{
    // Only change positive dimensions 
    // XXX numberOfVoxels must be an integer, warn the user

    if (sizeX > 0)
	{
	  if (sizeX > 2*detectorSizeX)
	   {
		   G4cout << "WARNING: Voxel X size must be smaller or equal than that of detector X" << G4endl;
		   return false;
	   }
	 // Round to zero 
	  numberOfVoxelsAlongX = (G4int) (2 * detectorSizeX / sizeX);
          sizeOfVoxelAlongX = (2 * detectorSizeX / numberOfVoxelsAlongX );
	}

    if (sizeY > 0)
	{
	  if (sizeY > 2*detectorSizeY)
	   {
		   G4cout << "WARNING: Voxel Y size must be smaller or equal than that of detector Y" << G4endl;
		   return false;
	   }
	 numberOfVoxelsAlongY = (G4int) (2 * detectorSizeY / sizeY);
	 sizeOfVoxelAlongY = (2 * detectorSizeY / numberOfVoxelsAlongY );
	}
	  if (sizeZ > 0)
	{
	  if (sizeZ > 2*detectorSizeZ)
	   {
		G4cout << "WARNING: Voxel Z size must be smaller or equal than that of detector Z" << G4endl;
		return false;
	   }
	 numberOfVoxelsAlongZ = (G4int) (2 * detectorSizeZ / sizeZ);
	 sizeOfVoxelAlongZ = (2 * detectorSizeZ / numberOfVoxelsAlongZ );
	}

    G4cout << "The (X,Y,Z) sizes of the Voxels are: (" << 
		G4BestUnit(sizeOfVoxelAlongX, "Length")  << ',' << 
		G4BestUnit(sizeOfVoxelAlongY, "Length")  << ',' << 
		G4BestUnit(sizeOfVoxelAlongZ, "Length") << ')' << G4endl;

    G4cout << "The number of Voxels along (X,Y,Z) is: (" << 
		numberOfVoxelsAlongX  << ',' <<
	        numberOfVoxelsAlongY  <<','  <<
		numberOfVoxelsAlongZ  << ')' << G4endl;

    //  This will clear the existing matrix (together with data inside it)! 
    matrix = HadrontherapyMatrix::getInstance(numberOfVoxelsAlongX, 
					      numberOfVoxelsAlongY,
	   			 	      numberOfVoxelsAlongZ);

    // Here construct the Sensitive Detector and Read Out Geometry
    ConstructSensitiveDetector(GetDetectorToWorldPosition());
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    return true;
}
/////////////////////////////////////////////////////////////////////////////
G4bool HadrontherapyDetectorConstruction::SetDetectorSize(G4double sizeX, G4double sizeY, G4double sizeZ)
{
// Check that the detector stay inside the phantom
	if (sizeX > 0 && sizeX < sizeOfVoxelAlongX) {G4cout << "WARNING: Detector X size must be bigger than that of Voxel X" << G4endl; return false;}
	if (sizeY > 0 && sizeY < sizeOfVoxelAlongY) {G4cout << "WARNING: Detector Y size must be bigger than that of Voxel Y" << G4endl; return false;}
	if (sizeZ > 0 && sizeZ < sizeOfVoxelAlongZ) {G4cout << "WARNING: Detector Z size must be bigger than that of Voxel Z" << G4endl; return false;}

	if (!IsInside(sizeX/2, 
		      sizeY/2, 
		      sizeZ/2, 
		      phantomSizeX, 
		      phantomSizeY,
		      phantomSizeZ,
		      detectorToPhantomPosition))
	{return false;}
// Negative or null values mean don't change it!
    if (sizeX > 0) {
	                detectorSizeX = sizeX/2;
		        detector -> SetXHalfLength(detectorSizeX);
		   }

    if (sizeY > 0) {
			detectorSizeY = sizeY/2;
		        detector -> SetYHalfLength(detectorSizeY);
                   }

    if (sizeZ > 0) {
			detectorSizeZ = sizeZ/2;
			detector -> SetZHalfLength(detectorSizeZ);
		    }


    G4cout << "The (X,Y,Z) dimensions of the detector are : (" << 
	    	  G4BestUnit( detector -> GetXHalfLength()*2., "Length") << ',' << 
	    	  G4BestUnit( detector -> GetYHalfLength()*2., "Length") << ',' << 
	    	  G4BestUnit( detector -> GetZHalfLength()*2., "Length") << ')' << G4endl; 
// Adjust detector position
    SetDetectorPosition();
// Adjust voxels number accordingly to new detector geometry 
// Matrix will be re-instantiated!
// Voxels and ROGeometry must follow the detector!
    SetNumberOfVoxelBySize(sizeOfVoxelAlongX, sizeOfVoxelAlongY, sizeOfVoxelAlongZ); 
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    return true;
}

/////////////////////////////////////////////////////////////////////////////
G4bool HadrontherapyDetectorConstruction::SetPhantomSize(G4double sizeX, G4double sizeY, G4double sizeZ)
{

	if (!IsInside(detectorSizeX, 
				  detectorSizeY, 
				  detectorSizeZ,
				  sizeX/2,//method parameters
				  sizeY/2,
				  sizeZ/2,
				  detectorToPhantomPosition
				  ))
	return false;

// Only change positive dimensions 
    if (sizeX > 0) {
                     phantomSizeX = sizeX/2;
                     phantom -> SetXHalfLength(phantomSizeX);
		   }
    if (sizeY > 0) {
                     phantomSizeY = sizeY/2;
                     phantom -> SetYHalfLength(phantomSizeY);
		   }

    if (sizeZ > 0) {
                     phantomSizeZ = sizeZ/2;
                     phantom -> SetZHalfLength(phantomSizeZ);
		   }
 

    G4cout << "The (X,Y,Z) dimensions of the phantom are : (" << 
	    	  G4BestUnit( phantom -> GetXHalfLength()*2., "Length") << ',' << 
	    	  G4BestUnit( phantom -> GetYHalfLength()*2., "Length") << ',' << 
	    	  G4BestUnit( phantom -> GetZHalfLength()*2., "Length") << ')' << G4endl; 
//G4cout << '\n' << "Coordinate volume: " << phantomPhysicalVolume -> GetTranslation() << G4endl; 
// Adjust detector position inside phantom
    SetDetectorPosition();

    ConstructSensitiveDetector(GetDetectorToWorldPosition());
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    return true;
}
/////////////////////////////////////////////////////////////////////////////
G4bool HadrontherapyDetectorConstruction::SetPhantomPosition(G4ThreeVector displacement)
{
// Set Phantom position respect to the World 
// TODO check limits!
    phantomPosition = displacement;
    if (phantomPhysicalVolume) 
	{
	    phantomPhysicalVolume -> SetTranslation(phantomPosition);
	    G4cout << "Displacement between Phantom and World is: "; 
	    G4cout << "DX= "<< G4BestUnit(phantomPosition.getX(),"Length") << 
		      "DY= "<< G4BestUnit(phantomPosition.getY(),"Length") << 
		      "DZ= "<< G4BestUnit(phantomPosition.getZ(),"Length") << G4endl;

// Redraw ROGeometry!
	    ConstructSensitiveDetector(GetDetectorToWorldPosition());
	    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
	}
    return true;
}

/////////////////////////////////////////////////////////////////////////////
G4bool HadrontherapyDetectorConstruction::SetDetectorToPhantomPosition(G4ThreeVector displacement)
{
// Ignore negative values
    if (displacement[0] < 0.) displacement[0] = detectorToPhantomPosition[0];
    if (displacement[1] < 0.) displacement[1] = detectorToPhantomPosition[1];
    if (displacement[2] < 0.) displacement[2] = detectorToPhantomPosition[2];

    if (!IsInside(detectorSizeX, 
		  detectorSizeY, 
		  detectorSizeZ,
		  phantomSizeX,
		  phantomSizeY,
		  phantomSizeZ,
		  displacement // method parameter!
		  ))
    {return false;}
    detectorToPhantomPosition = displacement;

// Adjust detector position inside phantom
    SetDetectorPosition();

    ConstructSensitiveDetector(GetDetectorToWorldPosition());
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    return true;
}
