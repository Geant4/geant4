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

#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"
#include "G4Colour.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4NistManager.hh"

#include "GammaKnifeDetectorMessenger.hh"
#include "GammaKnifeDetectorConstruction.hh"

#include "G4SystemOfUnits.hh"

GammaKnifeDetectorConstruction::GammaKnifeDetectorConstruction()
  : physicalTreatmentRoom(0),
    patientPhysicalVolume(0), 
    patientLogicalVolume(0),
    solidColl_helmet(0),
    helmetSize(4)
{
  // Messenger to change parameters of the geometry
  detectorMessenger = new GammaKnifeDetectorMessenger(this);
}

GammaKnifeDetectorConstruction::~GammaKnifeDetectorConstruction()
{
  delete detectorMessenger;
}

G4VPhysicalVolume* GammaKnifeDetectorConstruction::Construct()
{
  // Define the geometry components
  ConstructBeamLine();
 
  return physicalTreatmentRoom;
}

void GammaKnifeDetectorConstruction::ConstructBeamLine()
{
    // NIST Materials
    G4Material* air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    G4Material* water = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
    G4Material* cobalt = G4NistManager::Instance()->FindOrBuildMaterial("G4_Co");
    G4Material* Pb = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");
    G4Material* tungsten = G4NistManager::Instance()->FindOrBuildMaterial("G4_W");
    G4Material* Al = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
    G4Material* Fe = G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe");

    // Steel as non-NIST material
    G4Element* elFe = G4NistManager::Instance()->FindOrBuildElement("Fe");
    G4Element* elNi = G4NistManager::Instance()->FindOrBuildElement("Ni");
    G4Element* elCr = G4NistManager::Instance()->FindOrBuildElement("Cr");
    G4Material* steel = new G4Material("StainlessSteel", 7.80 * g/cm3, 3 /* components */);
    steel -> AddElement(elFe, 70 * perCent);
    steel -> AddElement(elCr, 18 * perCent);
    steel -> AddElement(elNi, 12 * perCent);

  // -----------------------------
  // Treatment room - World volume
  // -----------------------------

  // Treatment room sizes
  const G4double worldX = 400.0 *cm;
  const G4double worldY = 400.0 *cm;
  const G4double worldZ = 400.0 *cm;

  G4Box* treatmentRoom = new G4Box("TreatmentRoom",worldX,worldY,worldZ);

  G4LogicalVolume* logicTreatmentRoom = new G4LogicalVolume(treatmentRoom, 
                                                            air, 
                                                            "logicTreatmentRoom", 
							    0,0,0);

  physicalTreatmentRoom = new G4PVPlacement(0,
					    G4ThreeVector(),
					    "physicalTreatmentRoom", 
					    logicTreatmentRoom, 
					    0,false,0);


  // The treatment room is invisible in the Visualisation
  logicTreatmentRoom -> SetVisAttributes (G4VisAttributes::GetInvisible());


  // Visualisation attributes of all elements colours 
  G4VisAttributes * grayFe = new G4VisAttributes(G4Colour(0.5 ,0.5 ,0.5));
  grayFe -> SetVisibility(true);
  grayFe -> SetForceSolid(true);

  G4VisAttributes * blueCobalt = new G4VisAttributes(G4Colour(0. ,0. ,0.7));
  blueCobalt -> SetVisibility(true);
  blueCobalt -> SetForceSolid(true);

  G4VisAttributes * graySS = new G4VisAttributes(G4Colour(0.9 ,0.9 ,0.9));
  graySS -> SetVisibility(true);
  graySS -> SetForceSolid(true);

  G4VisAttributes * grayAl = new G4VisAttributes(G4Colour(0.7 ,0.7 ,0.7));
  grayAl -> SetVisibility(true);
  grayAl -> SetForceSolid(true);

  G4VisAttributes * blackLead = new G4VisAttributes(G4Colour(0.2 ,0.2 ,0.2));
  blackLead -> SetVisibility(true);
  blackLead -> SetForceSolid(true);
 

  G4VisAttributes * colorTungsten = new G4VisAttributes(G4Colour(0.3 ,0.3 ,0.3));
  colorTungsten -> SetVisibility(true);
  colorTungsten -> SetForceSolid(true);
 




  //--------------------------------------------
  // Cylinder source "Tube_source"
  //--------------------------------------------
  G4double innerRadiusOfTheTube_source = 0.;
  G4double outerRadiusOfTheTube_source = 0.5*mm;
  G4double hightOfTheTube_source = 1*cm;
  G4double startAngleOfTheTube = 0.*deg;
  G4double spanningAngleOfTheTube = 360.*deg;

  G4ThreeVector positionTube_source = G4ThreeVector(0,0,-40.1*cm);
  
  solidTube_source = new G4Tubs("solidTube_source",    
				 innerRadiusOfTheTube_source,    
				 outerRadiusOfTheTube_source,   
				 hightOfTheTube_source,       
				 startAngleOfTheTube,             
				 spanningAngleOfTheTube);          
  logicTube_source = new G4LogicalVolume(solidTube_source,cobalt,"logicTube_source",0,0,0);
  physiTube_source = new G4PVPlacement(0,           
				positionTube_source,
				logicTube_source,
				"Tube_source",
				logicTreatmentRoom,
				false,
				0);

  logicTube_source -> SetVisAttributes(blueCobalt);


  //-------------------------------------
  // Cylinder covering source "Tube"
  //-------------------------------------
  G4double innerRadiusOfTheTube = 0.5*mm;
  G4double outerRadiusOfTheTube = 4.*mm;
  G4double hightOfTheTube = 1*cm;
  

  G4ThreeVector positionTube = G4ThreeVector(0,0,-40.1*cm);
  
  solidTube = new G4Tubs("solidTube",    
				 innerRadiusOfTheTube,    
				 outerRadiusOfTheTube,  
				 hightOfTheTube,          
				 startAngleOfTheTube,    
				 spanningAngleOfTheTube);
  logicTube = new G4LogicalVolume(solidTube,steel,"logicTube",0,0,0);
  physiTube = new G4PVPlacement(0,          
				positionTube,
				logicTube,
				"Tube",
				logicTreatmentRoom,
				false,
				0);

  logicTube -> SetVisAttributes(graySS);

  //---------------------------------------
  // Cylinder covering source "Tube_Al"
  //---------------------------------------
  G4double innerRadiusOfTheTube_Al = 4.*mm;
  G4double outerRadiusOfTheTube_Al = 15.*mm;
  G4double hightOfTheTube_Al = 1*cm;
 
  G4ThreeVector positionTube_Al = G4ThreeVector(0,0,-40.1*cm);
  
  solidTube_Al = new G4Tubs("solidTube_Al",    
				 innerRadiusOfTheTube_Al,   
				 outerRadiusOfTheTube_Al,    
				 hightOfTheTube_Al,         
				 startAngleOfTheTube,     
				 spanningAngleOfTheTube); 
  logicTube_Al = new G4LogicalVolume(solidTube_Al,Al,"logicTube_Al",0,0,0);
  physiTube_Al = new G4PVPlacement(0,          
				positionTube_Al,
				logicTube_Al,
				"Tube_Al",
				logicTreatmentRoom,
				false,
				0);

  logicTube_Al -> SetVisAttributes(grayAl);

  //----------------------------------------------
  // Cylinder covering external part of the source "Tube_Fe"
  //----------------------------------------------
  G4double innerRadiusOfTheTube_Fe = 15.*mm;
  G4double outerRadiusOfTheTube_Fe = 50.*mm;
  G4double hightOfTheTube_Fe = 1*cm;
  

  G4ThreeVector positionTube_Fe = G4ThreeVector(0,0,-40.1*cm);
  
  solidTube_Fe = new G4Tubs("solidTube_Fe",    
				 innerRadiusOfTheTube_Fe,    
				 outerRadiusOfTheTube_Fe,   
				 hightOfTheTube_Fe,          
				 startAngleOfTheTube,     
				 spanningAngleOfTheTube); 
  logicTube_Fe = new G4LogicalVolume(solidTube_Fe,Fe,"logicTube_Fe",0,0,0);
  physiTube_Fe = new G4PVPlacement(0,           
				positionTube_Fe,
				logicTube_Fe,
				"Tube_Fe",
				logicTreatmentRoom,
				false,
				0);

  logicTube_Fe -> SetVisAttributes(grayFe);




 //------------------------------------------------
 // Cylinder covering posterior part "Tube_post"
 //------------------------------------------------

  G4double innerRadiusOfTheTube_post = 0;
  G4double outerRadiusOfTheTube_post = 50*mm;
  G4double hightOfTheTube_post = 1*cm;
 
  G4ThreeVector positionTube_post = G4ThreeVector(0,0,-42.2*cm);
  
  solidTube_post = new G4Tubs("solidTube_post",    
				 innerRadiusOfTheTube_post,   
				 outerRadiusOfTheTube_post,   
				 hightOfTheTube_post,          
				 startAngleOfTheTube,     
				 spanningAngleOfTheTube); 
 logicTube_post = new G4LogicalVolume(solidTube_post,Fe,"logicTube_post",0,0,0);
 physiTube_post = new G4PVPlacement(0,           
				positionTube_post,
				logicTube_post,
				"Tube_post",
			        logicTreatmentRoom,
				false,
				0);

 logicTube_post -> SetVisAttributes(grayFe);


 //------------------------------------------------
 // Fixed cylinder collimator "Tube_coll"
 //------------------------------------------------

  G4double innerRadiusOfTheTube_coll = 2.5*mm;
  G4double outerRadiusOfTheTube_coll = 15.*mm;
  G4double hightOfTheTube_coll = 3.25*cm;
 
  G4ThreeVector positionTube_coll = G4ThreeVector(0,0,-35.2*cm);
  
  solidTube_coll = new G4Tubs("solidTube_coll",    
				 innerRadiusOfTheTube_coll,    
				 outerRadiusOfTheTube_coll,    
				 hightOfTheTube_coll,         
				 startAngleOfTheTube,     
				 spanningAngleOfTheTube); 
 logicTube_coll = new G4LogicalVolume(solidTube_coll,tungsten,"logicTube_coll",0,0,0);
 physiTube_coll = new G4PVPlacement(0,           
				positionTube_coll,
				logicTube_coll,
				"Tube_coll",
				logicTreatmentRoom,
				false,
				0);

 logicTube_coll -> SetVisAttributes(colorTungsten);


 //------------------------------------------------
 // Cylinder covering fixed collimator "Tube_coll_Fe"
 //------------------------------------------------

  G4double innerRadiusOfTheTube_coll_Fe = 15.*mm;
  G4double outerRadiusOfTheTube_coll_Fe = 50.*mm;
  G4double hightOfTheTube_coll_Fe = 3.25*cm;
 
  G4ThreeVector positionTube_coll_Fe = G4ThreeVector(0,0,-35.2*cm);
  
  solidTube_coll_Fe = new G4Tubs("solidTube_coll_Fe",    
				 innerRadiusOfTheTube_coll_Fe,    
				 outerRadiusOfTheTube_coll_Fe,    
				 hightOfTheTube_coll_Fe,          
				 startAngleOfTheTube,     
				 spanningAngleOfTheTube); 
 logicTube_coll_Fe = new G4LogicalVolume(solidTube_coll_Fe,Fe,"logicTube_coll_Fe",0,0,0);
 physiTube_coll_Fe = new G4PVPlacement(0,          
				positionTube_coll_Fe,
				logicTube_coll_Fe,
				"Tube_coll_Fe",
				logicTreatmentRoom,
				false,
				0);

 logicTube_coll_Fe -> SetVisAttributes(grayFe);


 //------------------------------------------------
 // Fixed truncated cone collimator "Coll_fixed"
 //------------------------------------------------

  G4double Rmin1Coll_fixed = 2.5*mm;
  G4double Rmax1Coll_fixed = 15.*mm;
  G4double Rmin2Coll_fixed = 4.25*mm;
  G4double Rmax2Coll_fixed = 15.*mm;
  G4double hightColl_fixed = 4.625*cm;
 

  G4ThreeVector positionColl_fixed = G4ThreeVector(0,0,-27.325*cm);

  solidColl_fixed = new G4Cons("solidColl_fixed",
			       Rmin1Coll_fixed,  
			       Rmax1Coll_fixed, 
			       Rmin2Coll_fixed,  
			       Rmax2Coll_fixed, 
			       hightColl_fixed,  
			       startAngleOfTheTube,    
			       spanningAngleOfTheTube);
  logicColl_fixed = new G4LogicalVolume(solidColl_fixed,Pb,"logicColl_fixed",0,0,0);
  physiColl_fixed = new G4PVPlacement(0,
				      positionColl_fixed,
				      logicColl_fixed,
				      "Coll_fixed",
				      logicTreatmentRoom,
				      false,
				      0);

 logicColl_fixed -> SetVisAttributes(blackLead);


 //-----------------------------------------------------------
 // Cilinder covering fixed collimator "Coll_fixed_Fe"
 //-----------------------------------------------------------

  G4double Rmin1Coll_fixed_Fe = 15.*mm;
  G4double Rmax1Coll_fixed_Fe = 50.*mm;
  G4double Rmin2Coll_fixed_Fe = 15.*mm;
  G4double Rmax2Coll_fixed_Fe = 40.*mm;
  G4double hightColl_fixed_Fe = 4.625*cm;
 

  G4ThreeVector positionColl_fixed_Fe = G4ThreeVector(0,0,-27.325*cm);

  solidColl_fixed_Fe = new G4Cons("solidColl_fixed_Fe",
			       Rmin1Coll_fixed_Fe,  
			       Rmax1Coll_fixed_Fe,  
			       Rmin2Coll_fixed_Fe,  
			       Rmax2Coll_fixed_Fe,  
			       hightColl_fixed_Fe,  
			       startAngleOfTheTube,    // 
			       spanningAngleOfTheTube);
  logicColl_fixed_Fe = new G4LogicalVolume(solidColl_fixed_Fe,Fe,"logicColl_fixed_Fe",0,0,0);
  physiColl_fixed_Fe = new G4PVPlacement(0,
				      positionColl_fixed_Fe,
				      logicColl_fixed_Fe,
				      "Coll_fixed_Fe",
				      logicTreatmentRoom,
				      false,
				      0);

  logicColl_fixed_Fe -> SetVisAttributes(grayFe);


 //------------------------------------------------
 // Mobile truncate cone collimator "Coll_helmet"
 //------------------------------------------------
  G4double Rmax1Coll_helmet = 15.*mm;
  G4double Rmax2Coll_helmet = 15.*mm;
  G4double hightColl_helmet = 3.0*cm;
 

  G4ThreeVector positionColl_helmet = G4ThreeVector(0,0,-19.5*cm);

  solidColl_helmet = new G4Cons("solidColl_helmet",
                               0.0,  // will be set later
                               Rmax1Coll_helmet, 
                               0.0,  // will be set later
			       Rmax2Coll_helmet,  
			       hightColl_helmet, 
			       startAngleOfTheTube,   
                               spanningAngleOfTheTube);
  UpdateHelmet(); // Set the proper inner radii

  logicColl_helmet = new G4LogicalVolume(solidColl_helmet,tungsten,"logicColl_helmet",0,0,0);
  physiColl_helmet = new G4PVPlacement(0,
				      positionColl_helmet,
				      logicColl_helmet,
				      "Coll_helmet",
				      logicTreatmentRoom,
				      false,
				      0);

  logicColl_helmet -> SetVisAttributes(colorTungsten);

 //--------------------------------------------------------------
 // Truncated cone covering mobile collimator "Coll_helmet_Fe"
 //--------------------------------------------------------------

  G4double Rmin1Coll_helmet_Fe = 15.*mm;
  G4double Rmax1Coll_helmet_Fe = 40.*mm;
  G4double Rmin2Coll_helmet_Fe = 15.*mm;
  G4double Rmax2Coll_helmet_Fe = 30.*mm;
  G4double hightColl_helmet_Fe = 3.0*cm;

  G4ThreeVector positionColl_helmet_Fe = G4ThreeVector(0,0,-19.5*cm);

  solidColl_helmet_Fe = new G4Cons("solidColl_helmet_Fe",
			       Rmin1Coll_helmet_Fe, 
			       Rmax1Coll_helmet_Fe,  
			       Rmin2Coll_helmet_Fe,  
			       Rmax2Coll_helmet_Fe,  
			       hightColl_helmet_Fe, 
			       startAngleOfTheTube,   
			       spanningAngleOfTheTube);
  logicColl_helmet_Fe = new G4LogicalVolume(solidColl_helmet_Fe,Fe,"logicColl_helmet_Fe",0,0,0);
  physiColl_helmet_Fe = new G4PVPlacement(0,
				      positionColl_helmet_Fe,
				      logicColl_helmet_Fe,
				      "Coll_helmet_Fe",
				      logicTreatmentRoom,
				      false,
				      0);

  logicColl_helmet_Fe -> SetVisAttributes(grayFe);

  //-----------------------------------------
  // Patient --> water spherical phantom
  //-----------------------------------------

 
  G4Orb* patient = new G4Orb("patient",8.*cm);
  patientLogicalVolume = new G4LogicalVolume(patient,
							      water, 
							      "patientLog", 0, 0, 0);
  patientPhysicalVolume = new G4PVPlacement( new G4RotationMatrix(),
                                            G4ThreeVector(0., 0., 0.),
					    "patientPhys",
					    patientLogicalVolume,
					    physicalTreatmentRoom,
					    false,0);

  // Visualisation attributes of the patient 
  G4VisAttributes * redWire = new G4VisAttributes(G4Colour(0.8 ,0. ,0.));
  redWire -> SetVisibility(true);
  redWire -> SetForceWireframe(true);
  redWire -> SetForceAuxEdgeVisible(true);
  patientLogicalVolume -> SetVisAttributes(redWire); 

}

void GammaKnifeDetectorConstruction::UpdateHelmet()
{
    if (solidColl_helmet)
    {
        switch( helmetSize )
        {
        case 18:
            solidColl_helmet->SetInnerRadiusMinusZ( 4.15 * mm );
            solidColl_helmet->SetInnerRadiusPlusZ( 5.3  * mm );
            break;

        case 14:
            solidColl_helmet->SetInnerRadiusMinusZ( 3.15 * mm );
            solidColl_helmet->SetInnerRadiusPlusZ( 4.25  * mm );
            break;

        case 8:
            solidColl_helmet->SetInnerRadiusMinusZ( 1.9 * mm );
            solidColl_helmet->SetInnerRadiusPlusZ( 2.5  * mm );
            break;

        case 4:
            solidColl_helmet->SetInnerRadiusMinusZ( 1. * mm );
            solidColl_helmet->SetInnerRadiusPlusZ( 1.25  * mm );
            break;
        }
        // Inform the run manager about change in the geometry
        G4RunManager::GetRunManager()->GeometryHasBeenModified();
    }
}

void GammaKnifeDetectorConstruction::SetHelmetSize(G4int size)
{
    if (size != helmetSize) // Only if the size changes
    {
        // Allow only valid numbers
        switch( size )
        {
        case 18:
        case 14:
        case 8:
        case 4:
            helmetSize = size;
            G4cout << "Helmet size set to " << helmetSize << std::endl;
            UpdateHelmet();
            break;
        default:
      G4Exception("GammaKnifeDetectorConstruction::SetHelmetSize()", 
		  "GammaKnife001", FatalException, 
		  "Error: Invalid helmet size.");
            return;
        }
    }
}



