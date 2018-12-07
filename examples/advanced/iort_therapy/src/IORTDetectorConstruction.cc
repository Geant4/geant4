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
//   (e) University of Wollongong, Australia
//
//   *Corresponding author, email to carlo.casarino@polooncologicocefalu.it
//////////////////////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "G4SystemOfUnits.hh"
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
#include "G4Tubs.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSDoseDeposit3D.hh"
#include "IORTDetectorMessenger.hh"

/////////////////////////////////////////////////////////////////////////////
IORTDetectorConstruction::IORTDetectorConstruction(G4VPhysicalVolume* physicalTreatmentRoom)
  : motherPhys(physicalTreatmentRoom), // pointer to WORLD volume 
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

  /* NOTE! that the IORTDetectorConstruction class
   * does NOT inherit from G4VUserDetectorConstruction G4 class
   * So the Construct() mandatory virtual method is inside another geometric class
   * like the collimatorXXBeamLIne, ...
   */

  // Messenger to change parameters of the phantom/detector geometry
  detectorMessenger = new IORTDetectorMessenger(this);

  // Define here the material of the water phantom and of the detector
  SetPhantomMaterial("G4_WATER"); 

  // Construct geometry (messenger commands)

  // Detector
  // Default detector sizes
  detectorSizeX = 7.* cm;
  detectorSizeY = 15.* cm;
  detectorSizeZ = 15.* cm;

  SetDetectorSize(detectorSizeX,  detectorSizeY,  detectorSizeZ); 
 
  // Phantom 
  SetPhantomSize(20. *cm, 20. *cm, 20. *cm);   
  SetPhantomPosition(G4ThreeVector(4.5 *cm, 0. *cm, 0. *cm)); 
  SetDetectorToPhantomPosition(G4ThreeVector(0. *cm, 2.5 *cm, 2.5 *cm));  

  // Default protection disc geometry and materials
  SetOuterRadiusDiscoIORT (40. *mm);  
  SetinnerRadiusDiscoIORT (0.*mm);   
  SetheightDiscoIORT (2.0*mm);        
  SetDiscoXPositionIORT (-11.0*mm);
  SetDiscoYPositionIORT (0.0*mm);
  SetDiscoZPositionIORT (0.0*mm);
  SetDiscoMaterialIORT("G4_WATER");   

  SetOuterRadiusDiscoIORT1 (40. *mm);  
  SetinnerRadiusDiscoIORT1 (0.*mm);   
  SetheightDiscoIORT1 (1.0*mm);        
  SetDiscoXPositionIORT1 (-8.0*mm);
  SetDiscoMaterialIORT1("G4_WATER");

  SetAngleDiscoIORT0 (90.0 *deg);    

  // Write virtual parameters to the real ones and check for consistency      
  UpdateGeometry(); 
}

/////////////////////////////////////////////////////////////////////////////
IORTDetectorConstruction::~IORTDetectorConstruction()
{ 
    delete detectorMessenger;
}

/////////////////////////////////////////////////////////////////////////////
// ConstructPhantom() is the method that reconstuct a water box (called phantom 
// (or water phantom)). 
// A water phantom can be considered a good
// approximation of a an human body. 
////////////////////////////////////////////////////////////////////////////

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
    //phantomLogicalVolume -> SetVisAttributes(G4VisAttributes::GetInvisible());
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
   
    G4VisAttributes * skyBlue1 = new G4VisAttributes( G4Colour(135/255. , 206/255. ,  235/255. ));
    detectorLogicalVolume -> SetVisAttributes(skyBlue1);
   
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
 G4cout << "The Detector has been built --- Add a scoring mesh for it  in the GUI if appropriate (similar to the phantom one)" << G4endl;
 
}

void IORTDetectorConstruction::ConstructDisc()
{
// ---------------------------------------------------------------//
  //                    6.0 mm Protection Discs Volume          //
    // ---------------------------------------------------------------//
  const G4double startAngleDiscoIORT0 = 0.*deg;
  const G4double spanningAngleDiscoIORT0 = 360.*deg;

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

  const G4double startAngleDiscoIORT = 0.*deg;
  const G4double spanningAngleDiscoIORT = 360.*deg;
   
  G4double phi = 0. *deg;     

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
 
  const G4double startAngleDiscoIORT1 = 0.*deg;
  const G4double spanningAngleDiscoIORT1 = 360.*deg;    
  
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
}
/////////////////////////////////////////////////////////////////////////////

void IORTDetectorConstruction::SetVoxelSize(G4double , G4double , G4double)
{
    G4cout<< "SetVoxelSize method is not needed anymore " << G4endl;
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
    if (innerr >= 0.) {innerRadiusDiscoIORT = innerr;}
    
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
    if (innerr >= 0.) {innerRadiusDiscoIORT1 = innerr;}
    
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
   if (physiDiscoIORT1) delete physiDiscoIORT1;
   if (logicDiscoIORT1) delete logicDiscoIORT1;
   if (solidDiscoIORT1) delete solidDiscoIORT1;
   
   if (physiDiscoIORT) delete physiDiscoIORT;
   if (logicDiscoIORT) delete logicDiscoIORT;
   if (solidDiscoIORT) delete solidDiscoIORT;

   if (physiDiscoIORT0) delete physiDiscoIORT0;
   if (logicDiscoIORT0)  delete logicDiscoIORT0;
   if (solidDiscoIORT0) delete solidDiscoIORT0;
     
    ConstructDisc();

    // Inform the kernel about the new geometry
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
    G4RunManager::GetRunManager() -> PhysicsHasBeenModified();

    PrintParameters();
}

void IORTDetectorConstruction::DeleteDisc()
{ 
  if (physiDiscoIORT1) delete physiDiscoIORT1;
  if (logicDiscoIORT1) delete logicDiscoIORT1;
  if (solidDiscoIORT1) delete solidDiscoIORT1;
  
  if (physiDiscoIORT) delete physiDiscoIORT;
  if (logicDiscoIORT) delete logicDiscoIORT;
  if (solidDiscoIORT) delete solidDiscoIORT;
  
  if (physiDiscoIORT0) delete physiDiscoIORT0;
  if (logicDiscoIORT0) delete logicDiscoIORT0;
  if (solidDiscoIORT0) delete solidDiscoIORT0;
   
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
}

