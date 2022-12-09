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
/// \file FlashDetectorConstruction.cc
/// \brief Implementation of the FlashDetectorConstruction class


#include "FlashDetectorConstruction.hh"

#include "G4RunManager.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Region.hh"
#include "G4SDManager.hh"

#include "G4AutoDelete.hh"
#include "G4Box.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4Tubs.hh"

#include "G4GeometryManager.hh"
#include "G4GeometryTolerance.hh"

#include "G4UserLimits.hh"

#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "G4SystemOfUnits.hh"

#include "Applicator.hh"


#include "G4MaterialPropertiesTable.hh"

#include "G4PSEnergyDeposit.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4VisAttributes.hh"
#include "FlashDetectorMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FlashDetectorConstruction::FlashDetectorConstruction()
    : G4VUserDetectorConstruction(), physicalTreatmentRoom(0),logicTreatmentRoom(0), Collimator(0), fPhantom(0),
fPhantomLogicalVolume(0),fPhant_phys(0),
      fCheckOverlaps(true) {

  DefineMaterials();
  fDetectorMessenger = new FlashDetectorMessenger(this);

  SetPhantomSize(30. *cm, 30. *cm, 30. *cm);
  SetAirGap(0*cm); // Set the air gap between the water phantom and the end of the applicator
  SetDetectorThickness(10*um);
  SetDetector_subThickness(370*um);
  SetDetectorWidth(5*mm);
  SetAirGap_water_detector(0*cm); // Set the air gap between the end of the water phantom and the entrance of the detector

  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FlashDetectorConstruction::~FlashDetectorConstruction() {

  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FlashDetectorConstruction::DefineMaterials() {
  nist = G4NistManager::Instance();
//write here a function to define custom materials
G4bool isotopes = false;
 Si = nist->FindOrBuildElement("Si", isotopes);
 C = nist->FindOrBuildElement("C", isotopes);

  }


//
G4VPhysicalVolume *
FlashDetectorConstruction::ConstructPhantom(G4double CollPos) {
//This function creates a cubic phantom with the point Collpos on the surface of the cube.

 fPhantomMaterial = nist->FindOrBuildMaterial("G4_WATER");

 fPosition_coefficient = CollPos;
 
 fPhantom_coordinateX = (fPosition_coefficient * mm + fPhantomSizeX / 2);

 fPhantomPosition =  G4ThreeVector(fPhantom_coordinateX, 0. * mm, 0. * mm); //phantom is constructed with the entrance surface attached to the applicator 
  

  // Definition of the solid volume of the Phantom
  fPhantom = new G4Box("Phantom", fPhantomSizeX / 2, fPhantomSizeY / 2,
                      fPhantomSizeZ / 2);

  // Definition of the logical volume of the Phantom
  fPhantomLogicalVolume =
      new G4LogicalVolume(fPhantom, fPhantomMaterial, "phantomLog", 0, 0, 0);

  // Definition of the physical volume of the Phantom
  fPhant_phys =
      new G4PVPlacement(0, fPhantomPosition, "phantomPhys", fPhantomLogicalVolume,
                        physicalTreatmentRoom, false, 0);
//define the region to set cuts in FlashPhysicsList.cc and step limit
  G4Region *PhantomRegion = new G4Region("Phantom_reg");
  fPhantomLogicalVolume->SetRegion(PhantomRegion);
  PhantomRegion->AddRootLogicalVolume(fPhantomLogicalVolume);

  // Visualisation attributes of the phantom
  red = new G4VisAttributes(G4Colour(0 / 255., 255 / 255., 0 / 255.));
  red->SetVisibility(true);

  blue = new G4VisAttributes(G4Colour(0 / 255., 0. / 255., 255. / 255.));
  blue->SetVisibility(true);

  fPhantomLogicalVolume->SetVisAttributes(red);
//set step limit in phantom
  G4double maxStep = 0.1 * mm;
  fStepLimit = new G4UserLimits(maxStep);
  fPhantomLogicalVolume->SetUserLimits(fStepLimit);

  return fPhant_phys;
}

G4VPhysicalVolume *
FlashDetectorConstruction::ConstructDetector(){
 //Detector


  G4double fDensity_SiC=3.22*g/cm3;

  SiC=new G4Material("SiC", fDensity_SiC,2);
  SiC->AddElement(Si,1);
  SiC->AddElement(C,1);
 
 fDetectorMaterial=SiC;


 fDetectorPosition=fPhantom_coordinateX+fAirGap+fPhantomSizeX/2+fDet_thickness/2+fAirGap_phantom_det;
 
 fDet_box = new G4Box("Detector",fDet_thickness/2,fDet_width/2,fDet_width/2);
 
  // Definition of the logical volume of the Detector
  fDetLogicalVolume =
      new G4LogicalVolume(fDet_box, fDetectorMaterial, "DetectorLog", 0, 0, 0);
  fDet_phys = new G4PVPlacement(0,G4ThreeVector(fDetectorPosition, 0. * mm, 0. * mm), "DetPhys",fDetLogicalVolume,physicalTreatmentRoom,false, 0, fCheckOverlaps);

  
  fDet_sub = new G4Box("Det_sub",fDet_sub_thickness/2,fDet_width/2,fDet_width/2);
 
  // Definition of the logical volume of the Detector 
  fDet_sub_LogicalVolume =
      new G4LogicalVolume(fDet_sub, fDetectorMaterial, "Det_sub_Log", 0, 0, 0);
  fDet_sub_phys = new G4PVPlacement(0,G4ThreeVector(fDetectorPosition+fDet_thickness+fDet_sub_thickness/2, 0. * mm, 0. * mm), "Det_sub_Phys",fDet_sub_LogicalVolume,physicalTreatmentRoom,false, 0, fCheckOverlaps);


    return fDet_phys;

}

G4VPhysicalVolume *FlashDetectorConstruction::Construct() {
  // -----------------------------
  // Treatment room - World volume
  //------------------------------
  // Treatment room sizes
  const G4double worldX = 400.0 * cm;
  const G4double worldY = 400.0 * cm;
  const G4double worldZ = 400.0 * cm;
  G4bool isotopes = false;

  airNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);
  // Air
  //

  G4Box *treatmentRoom = new G4Box("TreatmentRoom", worldX, worldY, worldZ);
  logicTreatmentRoom = new G4LogicalVolume(treatmentRoom, airNist,
                                           "logicTreatmentRoom", 0, 0, 0);
  physicalTreatmentRoom =
      new G4PVPlacement(0, G4ThreeVector(), "physicalTreatmentRoom",
                        logicTreatmentRoom, 0, false, 0);

  // The treatment room is invisible in the Visualisation
  logicTreatmentRoom->SetVisAttributes(G4VisAttributes::GetInvisible());

  // -----------------------------
  // Applicator + phantom +Default dimensions
  //------------------------------
 

  

  Collimator = new Applicator(physicalTreatmentRoom);
  

   
 fPhantom_physical =
        ConstructPhantom(Collimator->fFinalApplicatorXPositionFlash +
	Collimator->fHightFinalApplicatorFlash+fAirGap);
 ConstructDetector();
   

  return physicalTreatmentRoom;
}



void FlashDetectorConstruction::ConstructSDandField() {
//modify this function if you want to insert a sensitive detector
}


/////MESSANGER ///

G4bool FlashDetectorConstruction::SetPhantomMaterial(G4String material)
{

    if (G4Material* pMat = G4NistManager::Instance()->FindOrBuildMaterial(material, false) )
    {
	fPhantomMaterial  = pMat;

	if (fPhantomLogicalVolume) 
	{
	    
	    fPhantomLogicalVolume ->  SetMaterial(pMat);

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


void FlashDetectorConstruction::SetPhantomSize(G4double sizeX, G4double sizeY, G4double sizeZ)
{
    if (sizeX > 0.) fPhantomSizeX = sizeX;
    if (sizeY > 0.) fPhantomSizeY = sizeY;
    if (sizeZ > 0.) fPhantomSizeZ = sizeZ;
}

void FlashDetectorConstruction::SetAirGap(G4double displ)
{
  
   fAirGap=displ;
}

G4bool FlashDetectorConstruction::SetDetectorMaterial(G4String material)
{

    if (G4Material* pMat = G4NistManager::Instance()->FindOrBuildMaterial(material, false) )
    {
	fDetectorMaterial  = pMat;

	if (fDetLogicalVolume) 
	{
	    
	    fDetLogicalVolume ->  SetMaterial(pMat);

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




void FlashDetectorConstruction::SetDetectorThickness(G4double thickness)
{
 
   fDet_thickness=thickness;
}

void FlashDetectorConstruction::SetDetectorWidth(G4double width)
{
 
   fDet_width=width;
}

void FlashDetectorConstruction::SetDetector_subThickness(G4double thickness_sub)
{
 
    fDet_sub_thickness= thickness_sub;
}


 


 void FlashDetectorConstruction::SetAirGap_water_detector(G4double spost)
{
  
  fAirGap_phantom_det=spost;
}


void FlashDetectorConstruction::SetDetectorPosition(G4double position)
{
 
   fDetectorPosition=position;
}
