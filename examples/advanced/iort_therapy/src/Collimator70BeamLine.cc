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

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"  
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4NistManager.hh"
#include "G4NistElementBuilder.hh"
#include "G4SubtractionSolid.hh"   
#include "IORTDetectorConstruction.hh" 
#include "Collimator70BeamLine.hh"
#include "Collimator70BeamLineMessenger.hh"  

Collimator70BeamLine::Collimator70BeamLine():
  physicalTreatmentRoom(0),iortDetectorConstruction(0),
  
  

  solidFinalCollimatorIORT(0),
  physiFinalCollimatorIORT(0),

  solidGiunz1FinalCollIORT(0),
  physiGiunz1FinalCollIORT(0),

  solidGiunz2FinalCollIORT(0),  
  physiGiunz2FinalCollIORT(0),

  solidGiunz3FinalCollIORT(0),
  physiGiunz3FinalCollIORT(0),  

  solidGiunz3FinalCollIntIORT(0),
  physiGiunz3FinalCollIntIORT(0), 
  
  solidGiunz4FinalCollIORT(0),
  physiGiunz4FinalCollIORT(0),

  solidGiunz5FinalCollIORT(0),
  physiGiunz5FinalCollIORT(0),
  
  solidBlocco1IORT(0),
  physiBlocco1IORT(0),

  solidBlocco2IORT(0),
  physiBlocco2IORT(0),

  solidBlocco3IORT(0),
  physiBlocco3IORT(0),

  solidBlocco20mmIORT(0),
  physiBlocco20mmIORT(0),  

  solidCM1_1_2IORT(0),
  physiCM1_1_2IORT(0),
  
  solidCM1_2_2IORT(0),
  physiCM1_2_2IORT(0),
  
  solidCM2_1_2IORT(0),
  physiCM2_1_2IORT(0),

  solidCM2_2_2IORT(0),
  physiCM2_2_2IORT(0),  

  solidCCMIORT(0),
  physiCCMIORT(0),

  solidPFS1IORT(0),
  physiPFS1IORT(0),

  solidPFS2IORT(0),
  physiPFS2IORT(0),

  solidPFS3IORT(0),
  physiPFS3IORT(0),

  solidFTIORT(0),
  physiFTIORT(0)


{
  // Messenger to change parameters of the collimator70BeamLine geometry
  collimatorMessenger = new Collimator70BeamLineMessenger(this);

}
/////////////////////////////////////////////////////////////////////////////
Collimator70BeamLine::~Collimator70BeamLine()
{
  delete collimatorMessenger;
  delete iortDetectorConstruction;
}

/////////////////////////////////////////////////////////////////////////////


G4VPhysicalVolume* Collimator70BeamLine::Construct()
{ 
  // Sets default geometry and materials
  SetDefaultDimensions();
  
  // Construct the whole Collimator Beam Line 
  ConstructCollimator70BeamLine();

  
  // IORTDetectorConstruction builds ONLY the phantom and the detector with its associated ROGeometry
  iortDetectorConstruction = new IORTDetectorConstruction(physicalTreatmentRoom); 
  
  return physicalTreatmentRoom;
}

// In the following method the DEFAULTS used in the geometry of 
// collimator beam line are provided
// HERE THE USER CAN CHANGE THE GEOMETRY CHARACTERISTICS OF BEAM
// LINE ELEMENTS, ALTERNATIVELY HE/SHE CAN USE THE MACRO FILE (IF A 
// MESSENGER IS PROVIDED)
//
// DEFAULT MATERIAL ARE ALSO PROVIDED	
// and COLOURS ARE ALSO DEFINED
// ----------------------------------------------------------
/////////////////////////////////////////////////////////////////////////////
void Collimator70BeamLine::SetDefaultDimensions()
{

   // Set of coulors that can be used
  white = new G4VisAttributes( G4Colour());
  white -> SetVisibility(true);
  //white -> SetForceSolid(true);
	
  blue = new G4VisAttributes(G4Colour(0. ,0. ,1.));
  blue -> SetVisibility(true);
  //blue -> SetForceSolid(true);
	
  gray = new G4VisAttributes( G4Colour(0.5, 0.5, 0.5 ));
  gray-> SetVisibility(true);
  //gray-> SetForceSolid(true);
	
  red = new G4VisAttributes(G4Colour(1. ,0. ,0.));
  red-> SetVisibility(true);
  //red-> SetForceSolid(true);
	
  yellow = new G4VisAttributes(G4Colour(1., 1., 0. ));
  yellow-> SetVisibility(true);
  //yellow-> SetForceSolid(true);
	
  green = new G4VisAttributes( G4Colour(25/255. , 255/255. ,  25/255. ));
  green -> SetVisibility(true);
  //green -> SetForceSolid(true);
	
  darkGreen = new G4VisAttributes( G4Colour(0/255. , 100/255. ,  0/255. ));
  darkGreen -> SetVisibility(true);
  //darkGreen -> SetForceSolid(true);
		
  darkOrange3 = new G4VisAttributes( G4Colour(205/255. , 102/255. ,  000/255. ));
  darkOrange3 -> SetVisibility(true);
  //darkOrange3 -> SetForceSolid(true);
	
  skyBlue = new G4VisAttributes( G4Colour(135/255. , 206/255. ,  235/255. ));
  skyBlue -> SetVisibility(true);
  //skyBlue -> SetForceSolid(true);
  

  // Geometry FINAL COLLIMATOR DEFAULTS

  G4double defaultOuterRadiusFinalCollimatorIORT = 40. *mm;
  OuterRadiusFinalCollimatorIORT = defaultOuterRadiusFinalCollimatorIORT;

  G4double defaultinnerRadiusFinalCollimatorIORT = 35. *mm;
  innerRadiusFinalCollimatorIORT = defaultinnerRadiusFinalCollimatorIORT;

  // DEFAULT DEFINITION OF THE MATERIALS
  // All elements and compound definition follows the NIST database
 
  // ELEMENTS
  G4bool isotopes = false;
  G4Material* aluminumNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al", isotopes);
  //G4Material* tantalumNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ta", isotopes); 
  //G4Material* copperNistAsMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu", isotopes);
  G4Element* zincNist = G4NistManager::Instance()->FindOrBuildElement("Zn");
  G4Element* copperNist = G4NistManager::Instance()->FindOrBuildElement("Cu");

  // COMPOUND
  G4Material* airNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);
  //G4Material* kaptonNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON", isotopes);
  G4Material* galacticNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic", isotopes);
  G4Material* PMMANist = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLEXIGLASS", isotopes);
  //G4Material* mylarNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_MYLAR", isotopes);
  G4Material* titanioNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ti", isotopes); 
  
  
  G4double d; // Density
  G4int nComponents;// Number of components 
  G4double fractionmass; // Fraction in mass of an element in a material

  d = 8.40*g/cm3;   // brass 
  nComponents = 2;
  G4Material* brass = new G4Material("Brass", d, nComponents);  
  brass -> AddElement(zincNist, fractionmass = 30 *perCent);
  brass -> AddElement(copperNist, fractionmass = 70 *perCent);

 
  // MATERIAL ASSIGNMENT


 // Material of the FINAL COLLIMATOR IORT
  finalCollimatorMaterialIORT = PMMANist;

 // Junction 1 FINAL COLLIMATOR IORT
  Giunz1FinalCollMaterialIORT = PMMANist;

 // Junction 2 FINAL COLLIMATOR IORT
  Giunz2FinalCollMaterialIORT = PMMANist;
 
 // Junction 3 FINAL COLLIMATOR IORT
  Giunz3FinalCollMaterialIORT = PMMANist;
 
 // Junction 3 FINAL COLLIMATOR Int IORT
  Giunz3FinalCollMaterialIntIORT = airNist;

 // Junction 4 FINAL COLLIMATOR IORT
  Giunz4FinalCollMaterialIORT = PMMANist;

 // Junction 5 FINAL COLLIMATOR IORT
  Giunz5FinalCollMaterialIORT = PMMANist;

 // Block 1 Diameter 30 mm 
  Blocco1IORTMaterialIORT = PMMANist; 

 // Block 2 Diameter 30 mm 
  Blocco2IORTMaterialIORT = PMMANist; 

 // Block 3 Diameter 30 mm 
  Blocco3IORTMaterialIORT = PMMANist;

 // Block Diameter 20 mm 
  Blocco20mmIORTMaterialIORT = PMMANist;

 // First Monitor Chamber Lamina Al 1 of 2  
    CM1_1_2IORTMaterialIORT = aluminumNist;

 // First Monitor Chamber Lamina Al 2 of 2  
    CM1_2_2IORTMaterialIORT = aluminumNist;

 // Second Monitor Chamber Lamina Al 1 of 2  
    CM2_1_2IORTMaterialIORT = aluminumNist;

 // Second Monitor Chamber Lamina Al 2 of 2  
    CM2_2_2IORTMaterialIORT = aluminumNist;
    
 // Monitor Chamber Cylinder 
    CCMIORTMaterialIORT = PMMANist;

 // Superior Final Part Monitor Chambers
    PFS1IORTMaterialIORT = PMMANist;

 // Superior Final Part Monitor Chambers
    PFS2IORTMaterialIORT = PMMANist;

 // Superior Final Part Monitor Chambers
    PFS3IORTMaterialIORT = PMMANist;

 // Superior Final Part Monitor Chambers Material
    FTIORTMaterialIORT = titanioNist;

 // Vacuum Source
    VSIORTMaterialIORT = galacticNist;

}

/////////////////////////////////////////////////////////////////////////////
void Collimator70BeamLine::ConstructCollimator70BeamLine()
{ 
  // -----------------------------
  // Treatment room - World volume
  //------------------------------
  // Treatment room sizes
  const G4double worldX = 400.0 *cm;
  const G4double worldY = 400.0 *cm;
  const G4double worldZ = 400.0 *cm;
  G4bool isotopes = false;
 
  G4Material* airNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);
  G4Box* treatmentRoom = new G4Box("TreatmentRoom",worldX,worldY,worldZ);
  G4LogicalVolume* logicTreatmentRoom = new G4LogicalVolume(treatmentRoom, 
                                                            airNist, 
                                                            "logicTreatmentRoom", 
							    0,0,0);
  physicalTreatmentRoom = new G4PVPlacement(0,
					    G4ThreeVector(),
					    "physicalTreatmentRoom", 
					    logicTreatmentRoom, 
					    0,false,0);
 

  // The treatment room is invisible in the Visualisation
  logicTreatmentRoom -> SetVisAttributes (G4VisAttributes::GetInvisible());
 
  // Components of the Collimator Beam Line

  IortBeamLineVacuumSource();
  IortBeamLineTitaniumWindows();
  IortBeamLineMonitorChambers();
  IortBeamLineBlocks() ;
  IortBeamLineJunctions(); 
  IortBeamLineFinalCollimator();
    
}


void Collimator70BeamLine::IortBeamLineVacuumSource()
{
 // ---------------------------------------------------------------//
  //                    Vacuum Source                            //
  // ---------------------------------------------------------------//

  
  G4double phi1 = 90. *deg;     

            
   G4RotationMatrix rm1;               
   rm1.rotateY(phi1);

  const G4double outRadiusVSIORT = 44.75 *mm;
  const G4double innRadiusVSIORT = 0.*mm;
  const G4double hightVSIORT = 1. *mm;
  const G4double startAngleVSIORT = 0.*deg;
  const G4double spanningAngleVSIORT = 360.*deg;
  const G4double XPositionVSIORT = -862.797 *mm;
    
  solidVSIORT = new G4Tubs("VSIORT", innRadiusVSIORT, 
				    outRadiusVSIORT,
				    hightVSIORT, 
				    startAngleVSIORT, 
				    spanningAngleVSIORT);

  G4LogicalVolume* logVSIORT = new G4LogicalVolume(solidVSIORT, 
							      VSIORTMaterialIORT, "VSIORT", 0, 0, 0);

  physiVSIORT = new G4PVPlacement(G4Transform3D(rm1, G4ThreeVector((XPositionVSIORT),0.,0.)),
					   "VSIORT", logVSIORT, physicalTreatmentRoom, false, 0); 

  logVSIORT -> SetVisAttributes(green);
}

void Collimator70BeamLine::IortBeamLineTitaniumWindows()
{
// ---------------------------------------------------------------//
  //                     Titanium Window                        //
  // ---------------------------------------------------------------//

  G4double phi2 = 90. *deg;     

            
   G4RotationMatrix rm2;               
   rm2.rotateY(phi2);


  const G4double outRadiusFTIORT = 44.75 *mm;
  const G4double innRadiusFTIORT = 8.5 *mm;
  const G4double hightFTIORT = 0.006 *mm;
  const G4double startAngleFTIORT = 0.*deg;
  const G4double spanningAngleFTIORT = 360.*deg;
  const G4double XPositionFTIORT = -861.791 *mm;

  solidFTIORT = new G4Tubs("FTIORT", innRadiusFTIORT, 
				    outRadiusFTIORT,
				    hightFTIORT, 
				    startAngleFTIORT, 
				    spanningAngleFTIORT);

  G4LogicalVolume* logFTIORT = new G4LogicalVolume(solidFTIORT, 
							      FTIORTMaterialIORT, "FTIORT", 0, 0, 0);

  physiFTIORT = new G4PVPlacement(G4Transform3D(rm2, G4ThreeVector((XPositionFTIORT),0.,0.)),
					   "FTIORT", logFTIORT, physicalTreatmentRoom, false, 0); 

  logFTIORT -> SetVisAttributes(yellow);
}

void Collimator70BeamLine::IortBeamLineMonitorChambers()
{

   G4double phi3 = 90. *deg;     

     // Matrix definition for a 90 deg rotation. Also used for other volumes       
   G4RotationMatrix rm3;               
   rm3.rotateY(phi3);
///////////////////////////////////////////////////////////////////////////////

  // Monitor Chambers System

///////////////////////////////////////////////////////////////////////////////
 

  // ---------------------------------------------------------------//
  //             Superior Final Part Monitor Chambers   3      //
  // ---------------------------------------------------------------//

  const G4double outRadiusPFS3IORT = 44.75 *mm;
  const G4double innRadiusPFS3IORT = 17.5 *mm;
  const G4double hightPFS3IORT = 3.03 *mm;
  const G4double startAnglePFS3IORT = 0.*deg;
  const G4double spanningAnglePFS3IORT = 360.*deg;
  const G4double XPositionPFS3IORT = -848.755 *mm;

  solidPFS3IORT = new G4Tubs("PFS3IORT", innRadiusPFS3IORT, 
				    outRadiusPFS3IORT,
				    hightPFS3IORT, 
				    startAnglePFS3IORT, 
				    spanningAnglePFS3IORT);

  G4LogicalVolume* logPFS3IORT = new G4LogicalVolume(solidPFS3IORT, 
							      PFS3IORTMaterialIORT, "PFS3IORT", 0, 0, 0);

  physiPFS3IORT = new G4PVPlacement(G4Transform3D(rm3, G4ThreeVector((XPositionPFS3IORT),0.,0.)),
					   "PFS3IORT", logPFS3IORT, physicalTreatmentRoom, false, 0); 

  logPFS3IORT -> SetVisAttributes(white);

  
  // ---------------------------------------------------------------//
  //             Superior Final Part Monitor Chambers   2       //
  // ---------------------------------------------------------------//

  const G4double outRadiusPFS2IORT = 44.75 *mm;
  const G4double innRadiusPFS2IORT = 10. *mm;
  const G4double hightPFS2IORT = 1.47 *mm;
  const G4double startAnglePFS2IORT = 0.*deg;
  const G4double spanningAnglePFS2IORT = 360.*deg;
  const G4double XPositionPFS2IORT = -844.255 *mm;

  solidPFS2IORT = new G4Tubs("PFS2IORT", innRadiusPFS2IORT, 
				    outRadiusPFS2IORT,
				    hightPFS2IORT, 
				    startAnglePFS2IORT, 
				    spanningAnglePFS2IORT);

  G4LogicalVolume* logPFS2IORT = new G4LogicalVolume(solidPFS2IORT, 
							      PFS2IORTMaterialIORT, "PFS2IORT", 0, 0, 0);

  physiPFS2IORT = new G4PVPlacement(G4Transform3D(rm3, G4ThreeVector((XPositionPFS2IORT),0.,0.)),
					   "PFS2IORT", logPFS2IORT, physicalTreatmentRoom, false, 0); 

  logPFS2IORT -> SetVisAttributes(green);

  // ---------------------------------------------------------------//
  //             Superior Final Part Monitor Chambers   1       //
  // ---------------------------------------------------------------//

  const G4double outRadiusPFS1IORT = 35. *mm;
  const G4double innRadiusPFS1IORT = 10. *mm;
  const G4double hightPFS1IORT = 0.88 *mm;
  const G4double startAnglePFS1IORT = 0.*deg;
  const G4double spanningAnglePFS1IORT = 360.*deg;
  const G4double XPositionPFS1IORT = -841.905 *mm;

  solidPFS1IORT = new G4Tubs("PFS1IORT", innRadiusPFS1IORT, 
				    outRadiusPFS1IORT,
				    hightPFS1IORT, 
				    startAnglePFS1IORT, 
				    spanningAnglePFS1IORT);

  G4LogicalVolume* logPFS1IORT = new G4LogicalVolume(solidPFS1IORT, 
							      PFS1IORTMaterialIORT, "PFS1IORT", 0, 0, 0);

  physiPFS1IORT = new G4PVPlacement(G4Transform3D(rm3, G4ThreeVector((XPositionPFS1IORT),0.,0.)),
					   "PFS1IORT", logPFS1IORT, physicalTreatmentRoom, false, 0); 

  logPFS1IORT -> SetVisAttributes(green);

  // ------------------------------------------------//
  //           Monitor Chambers Cylinder               //
  // ------------------------------------------------//

  const G4double outRadiusCCMIORT = 35. *mm;
  const G4double innRadiusCCMIORT = 10. *mm;
  const G4double hightCCMIORT = 4.0125 *mm;
  const G4double startAngleCCMIORT = 0.*deg;
  const G4double spanningAngleCCMIORT = 360.*deg;
  const G4double XPositionCCMIORT = -837.0125 *mm;

  solidCCMIORT = new G4Tubs("CCMIORT", innRadiusCCMIORT, 
				    outRadiusCCMIORT,
				    hightCCMIORT, 
				    startAngleCCMIORT, 
				    spanningAngleCCMIORT);

  G4LogicalVolume* logCCMIORT = new G4LogicalVolume(solidCCMIORT, 
							      CCMIORTMaterialIORT, "CCMIORT", 0, 0, 0);

  physiCCMIORT = new G4PVPlacement(G4Transform3D(rm3, G4ThreeVector((XPositionCCMIORT),0.,0.)),
					   "CCMIORT", logCCMIORT, physicalTreatmentRoom, false, 0); 

  logCCMIORT -> SetVisAttributes(green);


  // ------------------------------------------------//
  //        Second Monitor Chamber Lamina Al 2 of 2  //
  // ------------------------------------------------//

  const G4double outRadiusCM2_2_2IORT = 20. *mm;
  const G4double innRadiusCM2_2_2IORT = 0. *mm;
  const G4double hightCM2_2_2IORT = 0.025 *mm;
  const G4double startAngleCM2_2_2IORT = 0.*deg;
  const G4double spanningAngleCM2_2_2IORT = 360.*deg;
  const G4double XPositionCM2_2_2IORT = -841. *mm;

  solidCM2_2_2IORT = new G4Tubs("CM2_2_2IORT", innRadiusCM2_2_2IORT, 
				    outRadiusCM2_2_2IORT,
				    hightCM2_2_2IORT, 
				    startAngleCM2_2_2IORT, 
				    spanningAngleCM2_2_2IORT);

  G4LogicalVolume* logCM2_2_2IORT = new G4LogicalVolume(solidCM2_2_2IORT, 
							      CM2_2_2IORTMaterialIORT, "CM2_2_2IORT", 0, 0, 0);

  physiCM2_2_2IORT = new G4PVPlacement(G4Transform3D(rm3, G4ThreeVector((XPositionCM2_2_2IORT),0.,0.)),
					   "CM2_2_2ORT", logCM2_2_2IORT, physicalTreatmentRoom, false, 0); 

  logCM2_2_2IORT -> SetVisAttributes(green);  


// ------------------------------------------------//
  //        Second Monitor Chamber Lamina Al 1 of 2  //
  // ------------------------------------------------//

  const G4double outRadiusCM2_1_2IORT = 20. *mm;
  const G4double innRadiusCM2_1_2IORT = 0. *mm;
  const G4double hightCM2_1_2IORT = 0.025 *mm;
  const G4double startAngleCM2_1_2IORT = 0.*deg;
  const G4double spanningAngleCM2_1_2IORT = 360.*deg;
  const G4double XPositionCM2_1_2IORT = -839. *mm;

  solidCM2_1_2IORT = new G4Tubs("CM2_1_2IORT", innRadiusCM2_1_2IORT, 
				    outRadiusCM2_1_2IORT,
				    hightCM2_1_2IORT, 
				    startAngleCM2_1_2IORT, 
				    spanningAngleCM2_1_2IORT);

  G4LogicalVolume* logCM2_1_2IORT = new G4LogicalVolume(solidCM2_1_2IORT, 
							      CM2_1_2IORTMaterialIORT, "CM2_1_2IORT", 0, 0, 0);

  physiCM2_1_2IORT = new G4PVPlacement(G4Transform3D(rm3, G4ThreeVector((XPositionCM2_1_2IORT),0.,0.)),
					   "CM2_1_2ORT", logCM2_1_2IORT, physicalTreatmentRoom, false, 0); 

  logCM2_1_2IORT -> SetVisAttributes(yellow); 

  // ------------------------------------------------//
  //        First Monitor Chamber Lamina Al 2 of 2    //
  // ------------------------------------------------//

  const G4double outRadiusCM1_2_2IORT = 20. *mm;
  const G4double innRadiusCM1_2_2IORT = 0. *mm;
  const G4double hightCM1_2_2IORT = 0.025 *mm;
  const G4double startAngleCM1_2_2IORT = 0.*deg;
  const G4double spanningAngleCM1_2_2IORT = 360.*deg;
  const G4double XPositionCM1_2_2IORT = -837. *mm;

  solidCM1_2_2IORT = new G4Tubs("CM1_2_2IORT", innRadiusCM1_2_2IORT, 
				    outRadiusCM1_2_2IORT,
				    hightCM1_2_2IORT, 
				    startAngleCM1_2_2IORT, 
				    spanningAngleCM1_2_2IORT);

  G4LogicalVolume* logCM1_2_2IORT = new G4LogicalVolume(solidCM1_2_2IORT, 
							  CM1_2_2IORTMaterialIORT, "CM1_2_2IORT", 0, 0, 0);

  physiCM1_2_2IORT = new G4PVPlacement(G4Transform3D(rm3, G4ThreeVector((XPositionCM1_2_2IORT),0.,0.)),
					   "CM1_2_2ORT", logCM1_2_2IORT, physicalTreatmentRoom, false, 0); 

  logCM1_2_2IORT -> SetVisAttributes(yellow);
  
  // ------------------------------------------------//
  //        First Monitor Chamber Lamina Al 1 of 2         //
  // ------------------------------------------------//

  const G4double outRadiusCM1_1_2IORT = 20. *mm;
  const G4double innRadiusCM1_1_2IORT = 0. *mm;
  const G4double hightCM1_1_2IORT = 0.025 *mm;
  const G4double startAngleCM1_1_2IORT = 0.*deg;
  const G4double spanningAngleCM1_1_2IORT = 360.*deg;
  const G4double XPositionCM1_1_2IORT = -835. *mm;

  solidCM1_1_2IORT = new G4Tubs("CM1_1_2IORT", innRadiusCM1_1_2IORT, 
				    outRadiusCM1_1_2IORT,
				    hightCM1_1_2IORT, 
				    startAngleCM1_1_2IORT, 
				    spanningAngleCM1_1_2IORT);

  G4LogicalVolume* logCM1_1_2IORT = new G4LogicalVolume(solidCM1_1_2IORT, 
							      CM1_1_2IORTMaterialIORT, "CM1_1_2IORT", 0, 0, 0);

  physiCM1_1_2IORT = new G4PVPlacement(G4Transform3D(rm3, G4ThreeVector((XPositionCM1_1_2IORT),0.,0.)),
					   "CM1_1_2ORT", logCM1_1_2IORT, physicalTreatmentRoom, false, 0); 

  logCM1_1_2IORT -> SetVisAttributes(yellow);
}

void Collimator70BeamLine::IortBeamLineBlocks()
{

   G4double phi4 = 90. *deg;     

           
   G4RotationMatrix rm4;               
   rm4.rotateY(phi4);

 ///////////////////////////////////////////////////////////////////////////////

  // IORT BEAM LINE BLOCKS
  
///////////////////////////////////////////////////////////////////////////////

  // ------------------------------------------------//
  //        Block 4       //
  // ------------------------------------------------//

  const G4double outRadiusBlocco20mmIORT = 36.5 *mm;
  const G4double innRadiusBlocco20mmIORT = 10. *mm;
  const G4double hightBlocco20mmIORT = 3. *mm;
  const G4double startAngleBlocco20mmIORT = 0.*deg;
  const G4double spanningAngleBlocco20mmIORT = 360.*deg;
  const G4double XPositionBlocco20mmIORT = -830. *mm;

  solidBlocco20mmIORT = new G4Tubs("Blocco20mmIORT", innRadiusBlocco20mmIORT, 
				    outRadiusBlocco20mmIORT,
				    hightBlocco20mmIORT, 
				    startAngleBlocco20mmIORT, 
				    spanningAngleBlocco20mmIORT);

  G4LogicalVolume* logBlocco20mmIORT = new G4LogicalVolume(solidBlocco20mmIORT, 
							      Blocco20mmIORTMaterialIORT, "Blocco20mmIORT", 0, 0, 0);

  physiBlocco20mmIORT = new G4PVPlacement(G4Transform3D(rm4, G4ThreeVector((XPositionBlocco20mmIORT),0.,0.)),
					   "Blocco20mmORT", logBlocco20mmIORT, physicalTreatmentRoom, false, 0); 

  logBlocco20mmIORT -> SetVisAttributes(green);


  // -----------------------//
  //        Block 3        //
  // -----------------------//

  const G4double outRadiusBlocco3IORT = 36.5 *mm;
  const G4double innRadiusBlocco3IORT = 15. *mm;
  const G4double hightBlocco3IORT = 3.5 *mm;
  const G4double startAngleBlocco3IORT = 0.*deg;
  const G4double spanningAngleBlocco3IORT = 360.*deg;
  const G4double XPositionBlocco3IORT = -823.5 *mm;

  solidBlocco3IORT = new G4Tubs("Blocco3IORT", innRadiusBlocco3IORT, 
				    outRadiusBlocco3IORT,
				    hightBlocco3IORT, 
				    startAngleBlocco3IORT, 
				    spanningAngleBlocco3IORT);

  G4LogicalVolume* logBlocco3IORT = new G4LogicalVolume(solidBlocco3IORT, 
							      Blocco3IORTMaterialIORT, "Blocco3IORT", 0, 0, 0);

  physiBlocco3IORT = new G4PVPlacement(G4Transform3D(rm4, G4ThreeVector((XPositionBlocco3IORT),0.,0.)),
					   "Blocco3ORT", logBlocco3IORT, physicalTreatmentRoom, false, 0); 

  logBlocco3IORT -> SetVisAttributes(yellow);

 // -----------------------//
  //        Block 2        //
  // -----------------------//

  const G4double outRadiusBlocco2IORT = 41.5 *mm;
  const G4double innRadiusBlocco2IORT = 15. *mm;
  const G4double hightBlocco2IORT = 8. *mm;
  const G4double startAngleBlocco2IORT = 0.*deg;
  const G4double spanningAngleBlocco2IORT = 360.*deg;
  const G4double XPositionBlocco2IORT = -812. *mm;

  solidBlocco2IORT = new G4Tubs("Blocco2IORT", innRadiusBlocco2IORT, 
				    outRadiusBlocco2IORT,
				    hightBlocco2IORT, 
				    startAngleBlocco2IORT, 
				    spanningAngleBlocco2IORT);

  G4LogicalVolume* logBlocco2IORT = new G4LogicalVolume(solidBlocco2IORT, 
							      Blocco2IORTMaterialIORT, "Blocco2IORT", 0, 0, 0);

  physiBlocco2IORT = new G4PVPlacement(G4Transform3D(rm4, G4ThreeVector((XPositionBlocco2IORT),0.,0.)),
					   "Blocco2IORT", logBlocco2IORT, physicalTreatmentRoom, false, 0); 

  logBlocco2IORT -> SetVisAttributes(red);

  // ----------------------- //
  //       Block 1          //
  // ----------------------- //

  const G4double outRadiusBlocco1IORT = 52.0 *mm;
  const G4double innRadiusBlocco1IORT = 15. *mm;
  const G4double hightBlocco1IORT = 8.5 *mm;
  const G4double startAngleBlocco1IORT = 0.*deg;
  const G4double spanningAngleBlocco1IORT = 360.*deg;
  const G4double XPositionBlocco1IORT = -795.5*mm;

  solidBlocco1IORT = new G4Tubs("Blocco1IORT", innRadiusBlocco1IORT, 
				    outRadiusBlocco1IORT,
				    hightBlocco1IORT, 
				    startAngleBlocco1IORT, 
				    spanningAngleBlocco1IORT);

  G4LogicalVolume* logBlocco1IORT = new G4LogicalVolume(solidBlocco1IORT, 
							      Blocco1IORTMaterialIORT, "Blocco1IORT", 0, 0, 0);

  physiBlocco1IORT = new G4PVPlacement(G4Transform3D(rm4, G4ThreeVector((XPositionBlocco1IORT),0.,0.)),
					   "Blocco1IORT", logBlocco1IORT, physicalTreatmentRoom, false, 0); 

  logBlocco1IORT -> SetVisAttributes(white);
}

void Collimator70BeamLine::IortBeamLineJunctions()
{
 

  G4double phi5 = 90. *deg;     

          
   G4RotationMatrix rm5;               
   rm5.rotateY(phi5);
// --------------------------------- //
  // Junction 5 FINAL COLLIMATOR IORT //
  // --------------------------------- //

  const G4double outRadiusGiunz5FinalCollIORT = 48.25 *mm;
  const G4double innRadiusGiunz5FinalCollIORT = 13.75 *mm;
  const G4double hightGiunz5FinalCollIORT = 3.5 *mm;
  const G4double startAngleGiunz5FinalCollIORT = 0.*deg;
  const G4double spanningAngleGiunz5FinalCollIORT = 360.*deg;
  const G4double Giunz5FinalCollXPositionIORT = -783.5 *mm;

  solidGiunz5FinalCollIORT = new G4Tubs("Giunz5FinalCollIORT", innRadiusGiunz5FinalCollIORT, 
				    outRadiusGiunz5FinalCollIORT,
				    hightGiunz5FinalCollIORT, 
				    startAngleGiunz5FinalCollIORT, 
				    spanningAngleGiunz5FinalCollIORT);

  G4LogicalVolume* logGiunz5FinalCollIORT = new G4LogicalVolume(solidGiunz5FinalCollIORT, 
							      Giunz5FinalCollMaterialIORT, "Giunz5FinalCollIORT", 0, 0, 0);

  physiGiunz5FinalCollIORT = new G4PVPlacement(G4Transform3D(rm5, G4ThreeVector((Giunz5FinalCollXPositionIORT),0.,0.)),
					   "Giunz5FinalCollIORT", logGiunz5FinalCollIORT, physicalTreatmentRoom, false, 0); 

  logGiunz5FinalCollIORT -> SetVisAttributes(yellow);

// --------------------------------- //
  // Junction 4 FINAL COLLIMATOR IORT //
  // --------------------------------- //

  const G4double outRadiusGiunz4FinalCollIORT = 42. *mm;
  const G4double innRadiusGiunz4FinalCollIORT = 13.75 *mm;
  const G4double hightGiunz4FinalCollIORT = 8.5 *mm;
  const G4double startAngleGiunz4FinalCollIORT = 0.*deg;
  const G4double spanningAngleGiunz4FinalCollIORT = 360.*deg;
  const G4double Giunz4FinalCollXPositionIORT = -771.5 *mm;

  solidGiunz4FinalCollIORT = new G4Tubs("Giunz4FinalCollIORT", innRadiusGiunz4FinalCollIORT, 
				    outRadiusGiunz4FinalCollIORT,
				    hightGiunz4FinalCollIORT, 
				    startAngleGiunz4FinalCollIORT, 
				    spanningAngleGiunz4FinalCollIORT);

  G4LogicalVolume* logGiunz4FinalCollIORT = new G4LogicalVolume(solidGiunz4FinalCollIORT, 
							      Giunz4FinalCollMaterialIORT, "Giunz4FinalCollIORT", 0, 0, 0);

  physiGiunz4FinalCollIORT = new G4PVPlacement(G4Transform3D(rm5, G4ThreeVector((Giunz4FinalCollXPositionIORT),0.,0.)),
					   "Giunz4FinalCollIORT", logGiunz4FinalCollIORT, physicalTreatmentRoom, false, 0); 

  logGiunz4FinalCollIORT -> SetVisAttributes(blue); 


  
 // --------------------------------- //
  // Junction 3 FINAL COLLIMATOR IORT //
  // --------------------------------- //
   
  const G4double outRadiusGiunz3FinalCollIORT = 42. *mm;
  const G4double innRadiusGiunz3FinalCollIORT = 0. *mm;
  const G4double hightGiunz3FinalCollIORT = 4.25 *mm;
  const G4double startAngleGiunz3FinalCollIORT = 0.*deg;
  const G4double spanningAngleGiunz3FinalCollIORT = 360.*deg;
  const G4double Giunz3FinalCollXPositionIORT = -758.75 *mm;

  solidGiunz3FinalCollIORT = new G4Tubs("Giunz3FinalCollIORT", innRadiusGiunz3FinalCollIORT, 
				    outRadiusGiunz3FinalCollIORT,
				    hightGiunz3FinalCollIORT, 
				    startAngleGiunz3FinalCollIORT, 
				    spanningAngleGiunz3FinalCollIORT);

  G4LogicalVolume* logicsolidGiunz3FinalCollIORT = new G4LogicalVolume(solidGiunz3FinalCollIORT, 
							      Giunz3FinalCollMaterialIORT, "Giunz3FinalCollIORT", 0, 0, 0);

  physiGiunz3FinalCollIORT = new G4PVPlacement(G4Transform3D(rm5, G4ThreeVector((Giunz3FinalCollXPositionIORT),0.,0.)),
					   "Giunz3FinalCollIORT", logicsolidGiunz3FinalCollIORT, physicalTreatmentRoom, false, 0); 

  logicsolidGiunz3FinalCollIORT -> SetVisAttributes(yellow);
  //  logicsolidGiunz3FinalCollIORT -> SetVisAttributes (G4VisAttributes::GetInvisible());



  // --------------------------------- //
  // Junction 3 FINAL COLLIMATOR IORT internal //
  // --------------------------------- //
   
 
     
  solidGiunz3FinalCollIntIORT = new G4Cons("Giunz3FinalCollIntIORT",0.*mm,13.75*mm,0.*mm,22.25*mm,4.25*mm,0.*deg,360.*deg);

  G4LogicalVolume* logicsolidGiunz3FinalCollIntIORT = new G4LogicalVolume(solidGiunz3FinalCollIntIORT, 
							      Giunz3FinalCollMaterialIntIORT, "Giunz3FinalCollIntIORT", 0, 0, 0);

  physiGiunz3FinalCollIntIORT = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),"Giunz3FinalCollIntIORT", logicsolidGiunz3FinalCollIntIORT,physiGiunz3FinalCollIORT, false, 0); 

  logicsolidGiunz3FinalCollIntIORT -> SetVisAttributes(yellow); 


// --------------------------------- //
  // Junction 2 FINAL COLLIMATOR IORT //
  // --------------------------------- //

  const G4double outRadiusGiunz2FinalCollIORT = 42. *mm;
  const G4double innRadiusGiunz2FinalCollIORT = 22.25 *mm;
  const G4double hightGiunz2FinalCollIORT = 5.75 *mm;
  const G4double startAngleGiunz2FinalCollIORT = 0.*deg;
  const G4double spanningAngleGiunz2FinalCollIORT = 360.*deg;
  const G4double Giunz2FinalCollXPositionIORT = -748.75 *mm;

  solidGiunz2FinalCollIORT = new G4Tubs("Giunz2FinalCollIORT", innRadiusGiunz2FinalCollIORT, 
				    outRadiusGiunz2FinalCollIORT,
				    hightGiunz2FinalCollIORT, 
				    startAngleGiunz2FinalCollIORT, 
				    spanningAngleGiunz2FinalCollIORT);

  G4LogicalVolume* logGiunz2FinalCollIORT = new G4LogicalVolume(solidGiunz2FinalCollIORT, 
							      Giunz2FinalCollMaterialIORT, "Giunz2FinalCollIORT", 0, 0, 0);

  physiGiunz2FinalCollIORT = new G4PVPlacement(G4Transform3D(rm5, G4ThreeVector((Giunz2FinalCollXPositionIORT),0.,0.)),
					   "Giunz2FinalCollIORT", logGiunz2FinalCollIORT, physicalTreatmentRoom, false, 0); 

  logGiunz2FinalCollIORT -> SetVisAttributes(red); 
 
// --------------------------------- //
  // Junction 1 FINAL COLLIMATOR IORT //
  // --------------------------------- //

  const G4double outRadiusGiunz1FinalCollIORT = 50. *mm;
  const G4double innRadiusGiunz1FinalCollIORT = 22.25 *mm;
  const G4double hightGiunz1FinalCollIORT = 10. *mm;
  const G4double startAngleGiunz1FinalCollIORT = 0.*deg;
  const G4double spanningAngleGiunz1FinalCollIORT = 360.*deg;
  const G4double Giunz1FinalCollXPositionIORT = -733.*mm;
   
  solidGiunz1FinalCollIORT = new G4Tubs("Giunz1FinalCollIORT", innRadiusGiunz1FinalCollIORT, 
				    outRadiusGiunz1FinalCollIORT,
				    hightGiunz1FinalCollIORT, 
				    startAngleGiunz1FinalCollIORT, 
				    spanningAngleGiunz1FinalCollIORT);

  G4LogicalVolume* logGiunz1FinalCollIORT = new G4LogicalVolume(solidGiunz1FinalCollIORT, 
							      Giunz1FinalCollMaterialIORT, "Giunz1FinalCollIORT", 0, 0, 0);

  physiGiunz1FinalCollIORT = new G4PVPlacement(G4Transform3D(rm5, G4ThreeVector((Giunz1FinalCollXPositionIORT),0.,0.)),
					   "Giunz1FinalCollIORT", logGiunz1FinalCollIORT, physicalTreatmentRoom, false, 0); 

  logGiunz1FinalCollIORT -> SetVisAttributes(gray); 

}

void Collimator70BeamLine::IortBeamLineFinalCollimator()
{
// -----------------------//
  // FINAL COLLIMATOR IORT  //
  //------------------------//

 // const G4double outRadiusFinalCollimatorIORT = 40. *mm;
 // const G4double innRadiusFinalCollimatorIORT = 35. *mm;
  const G4double hightFinalCollimatorIORT = 334. *mm;
  const G4double startAngleFinalCollimatorIORT = 0.*deg;
  const G4double spanningAngleFinalCollimatorIORT = 360.*deg;
  const G4double finalCollimatorXPositionIORT = -389.*mm;

  
 
  
  G4double phi6 = 90. *deg;     

           
   G4RotationMatrix rm6;               
   rm6.rotateY(phi6);

    
  solidFinalCollimatorIORT = new G4Tubs("FinalCollimatorIORT", innerRadiusFinalCollimatorIORT, 
				    OuterRadiusFinalCollimatorIORT,
				    hightFinalCollimatorIORT, 
				    startAngleFinalCollimatorIORT, 
				    spanningAngleFinalCollimatorIORT);

  G4LogicalVolume* logFinalCollimatorIORT = new G4LogicalVolume(solidFinalCollimatorIORT, 
							      finalCollimatorMaterialIORT, "FinalCollimatorIORT", 0, 0, 0);

  physiFinalCollimatorIORT = new G4PVPlacement(G4Transform3D(rm6, G4ThreeVector((finalCollimatorXPositionIORT),0.,0.)),
					   "FinalCollimatorIORT", logFinalCollimatorIORT, physicalTreatmentRoom, false, 0); 

  //  logFinalCollimatorIORT -> SetVisAttributes(G4VisAttributes::GetInvisible()); 
  logFinalCollimatorIORT -> SetVisAttributes(darkOrange3);
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////// MESSENGER ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


void Collimator70BeamLine::SetInnerRadiusFinalCollimatorIORT(G4double value)
{
  solidFinalCollimatorIORT -> SetInnerRadius(value);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout<<"Inner Radius of the final collimator IORT is (mm):"
	<< solidFinalCollimatorIORT -> GetInnerRadius()/mm
	<< G4endl; 
}

/////////////////////////////////////////////////////////////////////////

void Collimator70BeamLine::SetOuterRadiusFinalCollimatorIORT(G4double value)
{
  solidFinalCollimatorIORT -> SetOuterRadius(value);
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout<<"Outer Radius of the final collimator IORT is (mm):"
	<< solidFinalCollimatorIORT -> GetOuterRadius()/mm
	<< G4endl; 
}

/////////////////////////////////////////////////////////////////////////////



