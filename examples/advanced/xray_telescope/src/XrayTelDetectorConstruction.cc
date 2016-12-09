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
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 xray_telescope advanced example         *
// *                                                                    *
// * MODULE:            XrayTelDetectorConstruction.cc                  *     
// * -------                                                            *
// *                                                                    *
// * Version:           0.4                                             *
// * Date:              06/11/00                                        *
// * Author:            R Nartallo                                      *
// * Organisation:      ESA/ESTEC, Noordwijk, THe Netherlands           *
// *                                                                    *
// **********************************************************************
// 
// CHANGE HISTORY
// --------------
//
// 06.11.2000 R.Nartallo
// - First implementation of xray_telescope geometry
// - Based on Chandra and XMM models by R Nartallo, P Truscott, F Lei 
//   and P Arce
//
//
// **********************************************************************

#include "XrayTelDetectorConstruction.hh"

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

XrayTelDetectorConstruction::XrayTelDetectorConstruction()
{
  world_x = 2500.*cm;
  world_y = 2500.*cm;
  world_z = 2500.*cm;
}

XrayTelDetectorConstruction::~XrayTelDetectorConstruction()
{;}

G4VPhysicalVolume* XrayTelDetectorConstruction::Construct( )
{

  // Material: Vacuum
  G4Material* Vacuum = new G4Material("Vacuum",
				      1.0 , 1.01*g/mole, 1.0E-25*g/cm3,
				      kStateGas, 2.73*kelvin, 3.0E-18*pascal );

  // Visualization attributes
  G4VisAttributes* VisAttWorld= new G4VisAttributes( G4Colour(204/255.,255/255.,255/255.));

  // World
  G4Box * solidWorld = new G4Box( "world_S", world_x, world_y, world_z );
  G4LogicalVolume * logicalWorld = new G4LogicalVolume( solidWorld, // solid
							Vacuum,                          // material
							"world_L",                       // name 
							0,0,0);

  logicalWorld -> SetVisAttributes(VisAttWorld);

  // Physical volume
  physicalWorld= new G4PVPlacement( 0,
				    G4ThreeVector(),
				    "world_P",        // name (2nd constructor)
				    logicalWorld,     // logical volume
				    NULL,             // mother volume
				    false,            // no boolean operation
				    0);               // copy number

  // Make Invisible
  logicalWorld -> SetVisAttributes(G4VisAttributes::GetInvisible());

  // Construct geometry
  ConstructTelescope();
  ConstructFocalPlane();

  return physicalWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Construct Telescope

void XrayTelDetectorConstruction::ConstructTelescope()
{
  // Construct Mirror
  // Single shell mirror made of Nickel with thin Gold coating
  // Mirror made up of two cones approximating the parabolic section and
  // two cones approximating the hyperbolic section
  // The centre of the mirror is filled wiith a solid aluminium rod shaped as the
  // mirrors, so as to leave a constant BaffleGap distance from the mirror surface

  // Build materials
  G4Material* Ni = new G4Material("Nickel", 28., 58.6934*g/mole, 8.902*g/cm3);
  G4Material* Au = new G4Material("Gold", 79., 196.96654*g/mole, 19.300*g/cm3);
  G4Material* Al = new G4Material("Aluminium", 13., 26.98*g/mole, 2.700*g/cm3);

  // Visualization attributes
  G4VisAttributes* VisAttMirror = new G4VisAttributes(
						      G4Colour(0/255., 0/255.,255/255.));
  G4VisAttributes* VisAttAuCoating = new G4VisAttributes(
							 G4Colour(255/255., 255/255., 0/255.));
  G4VisAttributes* VisAttBaffle = new G4VisAttributes(
						      G4Colour(128/255., 128/255., 128/255.));

  // Rotation Matrix
  G4RotationMatrix *rotateMatrix = new G4RotationMatrix();
  rotateMatrix -> rotateY(90.*deg);

  // Construct cones to make  Mirror sections
  G4int i;
  G4double MirrorEnd[5] = { 34.9995975*cm, 34.8277209*cm, 34.6549918*cm,
			    34.1347834*cm, 33.6137753*cm };
  G4double MirrorPosition[4] = { 772.5*cm, 757.5*cm, 742.5*cm, 727.5*cm };
  G4double MirrorSectionLength = 15.0*cm;
  G4double MirrorNiThickness = 1.07*mm;
  G4double MirrorAuCoating = 50.0e-6*mm;
  G4double BaffleGap = 4.0*mm;

  G4Cons* MirrorSolid[4];
  G4Cons* MirrorAuCoatingSolid[4];
  G4Cons* BaffleSolid[4];

  G4LogicalVolume* MirrorLogicalVolume[4];
  G4LogicalVolume* MirrorAuCoatingLogicalVolume[4];
  G4LogicalVolume* BaffleLogicalVolume[4];

  for ( i=0; i<4; i++ ) {

    // Mirror Nickel base
    MirrorSolid[i] = new G4Cons( "Mirror_S",
				 MirrorEnd[i], MirrorEnd[i] + MirrorNiThickness,
				 MirrorEnd[i+1], MirrorEnd[i+1] + MirrorNiThickness,
				 MirrorSectionLength/2, 0*deg, 360.*deg);
    MirrorLogicalVolume[i] = new G4LogicalVolume(
						 MirrorSolid[i], Ni, "Mirror_L", 0, 0, 0 );
    MirrorLogicalVolume[i]->SetVisAttributes(VisAttMirror);

    // Gold coating on mirror
    MirrorAuCoatingSolid[i] = new G4Cons( 
					 "MirrorAuCoating_S",
					 MirrorEnd[i] - MirrorAuCoating, MirrorEnd[i],
					 MirrorEnd[i+1] - MirrorAuCoating, MirrorEnd[i+1], 
					 MirrorSectionLength/2, 0*deg, 360.*deg);
    MirrorAuCoatingLogicalVolume[i] = new G4LogicalVolume(
							  MirrorAuCoatingSolid[i],
							  Au,
							  "MirrorAuCoating_L",
							  0, 0, 0 );
    MirrorAuCoatingLogicalVolume[i]->SetVisAttributes(VisAttAuCoating);

    // Aluminium baffle inside mirror
    BaffleSolid[i] = new G4Cons( "Baffle_S",
				 0, MirrorEnd[i] - BaffleGap,
				 0, MirrorEnd[i+1] - BaffleGap,
				 MirrorSectionLength/2, 0*deg, 360.*deg);
    BaffleLogicalVolume[i] = new G4LogicalVolume(
						 BaffleSolid[i], Al, "Baffle_L", 0, 0, 0 );
    BaffleLogicalVolume[i]-> SetVisAttributes(VisAttBaffle);
  }

  // Physical volume
 
  for ( i=0; i<4; i++ ) {
     new G4PVPlacement(
		       rotateMatrix,
		       G4ThreeVector( MirrorPosition[i], 0.0*cm, 0.0*cm ),
		       "Mirror_P",
		       MirrorLogicalVolume[i],
		       physicalWorld, false, 0 );
     new G4PVPlacement(
		       rotateMatrix,
		       G4ThreeVector( MirrorPosition[i], 0.0*cm, 0.0*cm ),
		       "MirrorAuCoating_P",
		       MirrorAuCoatingLogicalVolume[i],
		       physicalWorld, false, 0 );
    new G4PVPlacement(
		      rotateMatrix,
		      G4ThreeVector( MirrorPosition[i], 0.0*cm, 0.0*cm ),
		      "Baffle_P",
		      BaffleLogicalVolume[i],
		      physicalWorld, false, 0 );
  }

  // Make Mirror Invisible

  for ( i=0; i<4; i++ ) {
    //   MirrorLogicalVolume[i] -> SetVisAttributes(G4VisAttributes::GetInvisible());
    //   MirrorAuCoatingLogicalVolume[i] -> SetVisAttributes(G4VisAttributes::GetInvisible());
    BaffleLogicalVolume[i] -> SetVisAttributes(G4VisAttributes::GetInvisible());
  }


  // Construct Optical Bench
  // Main Telescope carbon fibre tube and two aluminium end caps

  G4int nel;
  G4String symbol;

  // Elements
  G4Element* C = new G4Element("Carbon", symbol="C", 6., 12.011*g/mole);
  G4Element* H = new G4Element("Hydrogen",symbol="H", 1., 1.00794*g/mole);

  // Materials from Combination
  G4Material* Cf = new G4Material("Carbon Fibre", 2.0*g/cm3, nel=2);
  Cf->AddElement(C,1);
  Cf->AddElement(H,2);

  // Visualization attributes
  G4VisAttributes* VisAttBench = new G4VisAttributes(
						     G4Colour(0/255., 200/255., 0/255.));

  // Construct Optical bench
  G4double BenchThickness = 1.0*cm;
  G4double BenchFrontEndMinRadiusOut = MirrorEnd[4] + 
    ( MirrorEnd[3] - MirrorEnd[4] )*7.5/15
    + MirrorNiThickness;
  G4double BenchFrontEndMinRadiusIn  = MirrorEnd[4] + 
    ( MirrorEnd[3] - MirrorEnd[4] )*7.4/15
    + MirrorNiThickness;
  G4double BenchFrontEndMaxRadius = MirrorEnd[4] + MirrorNiThickness + 25.*cm;
  G4double BenchBackEndMinRadius = 0.0*cm;
  G4double BenchBackEndMaxRadius =  MirrorEnd[4] + MirrorNiThickness + 5.*cm;
  G4double BenchMainLength;

  BenchMainLength = MirrorPosition[3] - BenchThickness;

  G4Cons* BenchFrontEndSolid;
  G4Tubs* BenchBackEndSolid;
  G4Cons* BenchMainSolid;

  G4LogicalVolume* BenchFrontEndLogicalVolume;
  G4LogicalVolume* BenchBackEndLogicalVolume;
  G4LogicalVolume* BenchMainLogicalVolume;

  BenchFrontEndSolid = new G4Cons( "BenchFrontEnd_S",
				   BenchFrontEndMinRadiusOut, BenchFrontEndMaxRadius,
				   BenchFrontEndMinRadiusIn, BenchFrontEndMaxRadius,
				   BenchThickness/2, 0*deg, 360.*deg );
  BenchFrontEndLogicalVolume = new G4LogicalVolume(
						   BenchFrontEndSolid, Al, "BenchFrontEnd_L", 0, 0, 0 );
  BenchFrontEndLogicalVolume->SetVisAttributes(VisAttBench);

  BenchBackEndSolid  = new G4Tubs( "BenchBackEnd_S",
				   BenchBackEndMinRadius, BenchBackEndMaxRadius,
				   BenchThickness/2, 0*deg, 360.*deg );
  BenchBackEndLogicalVolume = new G4LogicalVolume(
						  BenchBackEndSolid, Al, "BenchBackEnd_L", 0, 0, 0 );
  BenchBackEndLogicalVolume->SetVisAttributes(VisAttBench);

  BenchMainSolid     = new G4Cons( "BenchMain_S",
				   BenchFrontEndMaxRadius - BenchThickness,
				   BenchFrontEndMaxRadius,
				   BenchBackEndMaxRadius - BenchThickness,
				   BenchBackEndMaxRadius,
				   BenchMainLength/2, 0*deg, 360.*deg);
  BenchMainLogicalVolume = new G4LogicalVolume(
					       BenchMainSolid, Cf, "BenchMain_L", 0, 0, 0 );
  BenchMainLogicalVolume -> SetVisAttributes(VisAttBench);

  // Physical volume

  new G4PVPlacement(
		    rotateMatrix,
		    G4ThreeVector( MirrorPosition[3] - BenchThickness/2,
				   0.0*cm, 0.0*cm ),
		    "BenchFrontEnd_P",
		    BenchFrontEndLogicalVolume,
		    physicalWorld, false, 0 );

  new G4PVPlacement(
		    rotateMatrix,
		    G4ThreeVector(0.0*cm - BenchThickness/2, 0.0*cm, 0.0*cm ),
		    "BenchBackEnd_P",
		    BenchBackEndLogicalVolume,
		    physicalWorld, false, 0 );
  
  new G4PVPlacement(
		    rotateMatrix,
		    G4ThreeVector( BenchMainLength/2, 0.0*cm, 0.0*cm ),
		    "BenchMain_P",
		    BenchMainLogicalVolume,
		    physicalWorld, false, 0 );
  
  //--- Make Bench Invisible

  // BenchFrontEndLogicalVolume -> SetVisAttributes(G4VisAttributes::GetInvisible())

  // BenchBackEndLogicalVolume -> SetVisAttributes(G4VisAttributes::GetInvisible());
  BenchMainLogicalVolume -> SetVisAttributes(G4VisAttributes::GetInvisible());

  return;
}

// Construct Focal Plane
// Conical Titanium baffle and silicon detector

void XrayTelDetectorConstruction::ConstructFocalPlane()
{

  // Elements
  G4Material* Ti = new G4Material("Titanium", 22., 47.867*g/mole, 4.54*g/cm3);
  G4Material* Si = new G4Material("Silicon", 14., 28.090*g/mole, 2.33*g/cm3);

  // Visualization attributes
  G4VisAttributes* VisDetectorBaffle = new G4VisAttributes(
							   G4Colour(190/255., 255/255., 0/255.) );
  G4VisAttributes* VisDetector = new G4VisAttributes(
						     G4Colour(255/255., 0/255., 0/255.) );

  // Rotation Matrix
  G4RotationMatrix *rotateMatrix = new G4RotationMatrix();
  rotateMatrix -> rotateY(90.*deg);

  // Construct Detector Baffle
  G4double DetectorBaffleLength = 57.2*cm;
  G4double DetectorBaffleOuterRadiusIn = 7.1*cm;
  G4double DetectorBaffleOuterRadiusOut = 7.35*cm;
  G4double DetectorBaffleInnerRadiusIn = 4.55*cm;
  G4double DetectorBaffleInnerRadiusOut = 5.75*cm;

  G4Cons* DetectorBaffleSolid;

  G4LogicalVolume* DetectorBaffleLogicalVolume;

  DetectorBaffleSolid = new G4Cons( "DetectorBaffle_S",              
				    DetectorBaffleOuterRadiusIn,
				    DetectorBaffleOuterRadiusOut,
				    DetectorBaffleInnerRadiusIn,
				    DetectorBaffleInnerRadiusOut,
				    DetectorBaffleLength/2, 0*deg, 360.*deg);  
  DetectorBaffleLogicalVolume = new G4LogicalVolume( 
						    DetectorBaffleSolid, Ti, "DetectorBaffle_L", 0, 0, 0 );
  DetectorBaffleLogicalVolume -> SetVisAttributes( VisDetectorBaffle );

  // Physical volume
 
  /* G4VPhysicalVolume* DetectorBafflePhysicalVolume = */
  new G4PVPlacement(
		    rotateMatrix,
		    G4ThreeVector( DetectorBaffleLength/2, 0.0*cm, 0.0*cm),
		    "DetectorBaffle_P",
		    DetectorBaffleLogicalVolume,
		    physicalWorld, false, 0 );
  
  //--- Make Invisible

  // DetectorBaffleLogicalVolume -> SetVisAttributes( G4VisAttributes::GetInvisible() );

  // Construct Detector

  G4double DetectorRadius = 32.5*mm;
  G4double DetectorThickness = 50e-6*m;

  G4Tubs* DetectorSolid;

  G4LogicalVolume* DetectorLogicalVolume;

  DetectorSolid = new G4Tubs( "Detector_S",                        
			      0, DetectorRadius,
			      DetectorThickness/2, 0*deg, 360.*deg);  
  DetectorLogicalVolume = new G4LogicalVolume( 
					      DetectorSolid, Si, "Detector_L", 0, 0, 0 );
  DetectorLogicalVolume -> SetVisAttributes( VisDetector );

  // Physical volume
  /*G4VPhysicalVolume* DetectorPhysicalVolume = */
  new G4PVPlacement( 
		    rotateMatrix,
		    G4ThreeVector( DetectorThickness/2, 0.0*cm, 0.0*cm),
		    "Detector_P",
		    DetectorLogicalVolume,
		    physicalWorld, false, 0 );
  
  //--- Make Invisible
  // DetectorLogicalVolume -> SetVisAttributes( G4VisAttributes::GetInvisible() );

  return;
}
