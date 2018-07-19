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
// --------------------------------------------------------------
//                 GEANT 4 - Brachytherapy example
// --------------------------------------------------------------
//
// Code developed by:
// S. Guatelli, D. Cutajar, J. Poder
// Centre For Medical Radiation Physics, University of Wollongong
//
//    ****************************************
//    *                                      *
//    *    BrachyDetectorConstructionFlexi.cc*
//    *                                      *
//    ****************************************
//
//
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "BrachyDetectorConstructionFlexi.hh"
#include "G4RunManager.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4TransportationManager.hh"
#include "BrachyMaterial.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

BrachyDetectorConstructionFlexi::BrachyDetectorConstructionFlexi()
  : steel_shell(0),logical_steel_shell(0),air_gap(0), logical_air_gap(0), physical_air_gap(0),
    End1_steel_shell(0),logical_End1_steel_shell(0), physical_End1_steel_shell(0),
    End2_steel_shell(0),logical_End2_steel_shell(0), physical_End2_steel_shell(0),
    cable(0),logical_cable(0),physical_cable(0),
    iridium_core(0),logical_iridium_core(0),physical_iridium_core(0),
    steelAttributes(0), endAttributes(0), simpleIridiumVisAtt(0)
{
  pMat = new BrachyMaterial();
}

BrachyDetectorConstructionFlexi::~BrachyDetectorConstructionFlexi()
{ 
  delete pMat; 
}

void BrachyDetectorConstructionFlexi::ConstructFlexi(G4VPhysicalVolume* mother)
{
  G4Material* steelMat = pMat -> GetMat("Stainless steel 304");
  G4Material* iridiumMat = pMat -> GetMat("Iridium");
  G4Material* airMat = pMat -> GetMat("Air");

 //Define dimensions of the outer Steel shell around the solid source - not including the ends 

  G4double shellr_min = 0.00 * mm;
  G4double shellr_max = 0.85 * mm;
  G4double shell_length = 3.6 * mm; 
  steel_shell = new G4Tubs("steel_shell",shellr_min, shellr_max/2, shell_length/2.,0.*deg,360.*deg);
  logical_steel_shell = new G4LogicalVolume(steel_shell, steelMat, "steel_shell_log", 0, 0, 0);
  physical_steel_shell = new G4PVPlacement(0,G4ThreeVector(0,0,0),"phys_steel_shell", logical_steel_shell, mother, false, 0, true);

//Define dimensions of the air gap between Steel shell and Iridium core
  G4double airr_min = 0.00 * mm;
  G4double airr_max = 0.67 * mm;
  G4double air_length = 3.6 * mm; 
  air_gap = new G4Tubs("air_gap", airr_min, airr_max/2, air_length/2, 0.*deg, 360.*deg);
  logical_air_gap = new G4LogicalVolume(air_gap, airMat, "air_gap_log", 0, 0, 0);
  physical_air_gap = new G4PVPlacement(0, G4ThreeVector(0,0,0), "phys_air_gap",logical_air_gap, physical_steel_shell, false, 0, true);

//Define the non-cable weld end of the Steel shell 
  G4double End1r_min = 0.0 * mm;
  G4double End1r_max = 0.85 * mm;
  G4double End1length = 0.65 * mm; 
  End1_steel_shell = new G4Tubs("End_1_steel_shell",End1r_min, End1r_max/2, End1length/2.,0.*deg,360.*deg);
  logical_End1_steel_shell = new G4LogicalVolume(End1_steel_shell, steelMat, "End1_steel_shell_log", 0, 0, 0);
  G4double end1offset_x = 0.0 * mm;
  G4double end1offset_y = 0.0 * mm;
  G4double end1offset_z = 2.125 * mm;
  physical_End1_steel_shell = new G4PVPlacement(0,G4ThreeVector(end1offset_x,end1offset_y,end1offset_z),"phys_End1_steel_shell", logical_End1_steel_shell,mother, false, 0, true);

//Define the cable weld end of the Steel shell 
  G4double End2r_min1 = 0.0 * mm;
  G4double End2r_max1 = 0.85 * mm;
  G4double End2r_min2 = 0.0 * mm;
  G4double End2r_max2 = 0.5 * mm;
  G4double End2length = 0.4 * mm;
  End2_steel_shell = new G4Cons("End_2_steel_shell",End2r_min2, End2r_max2/2, End2r_min1, End2r_max1/2, End2length/2.0, 0.0, 360.0*deg);
  logical_End2_steel_shell = new G4LogicalVolume(End2_steel_shell, steelMat, "End2_steel_shell_log", 0, 0, 0);
  G4double end2offset_x = 0.0 * mm;
  G4double end2offset_y = 0.0 * mm;
  G4double end2offset_z = -2.0 * mm; 
  physical_End2_steel_shell = new G4PVPlacement(0,G4ThreeVector(end2offset_x,end2offset_y,end2offset_z), "phys_End2_steel_shell", logical_End2_steel_shell,mother, false, 0, true);

//Define the cable 
  G4double cable_min = 0.0 * mm;
  G4double cable_max = 0.5 * mm;
  G4double cablelength = 5.0 * mm; 
  cable = new G4Tubs("cable",cable_min, cable_max/2, cablelength/2.,0.*deg,360.*deg);
  logical_cable = new G4LogicalVolume(cable, steelMat, "cable_log", 0, 0, 0);
  G4double cableoffset_x = 0.0 * mm;
  G4double cableoffset_y = 0.0 * mm;
  G4double cableoffset_z = -4.7 * mm;
  physical_cable = new G4PVPlacement(0,G4ThreeVector(cableoffset_x,cableoffset_y,cableoffset_z),"phys_cable", logical_cable, mother, false, 0, true);

// Define the Iridium core
  G4double corer_min = 0.0 * mm;	
  G4double corer_max = 0.6 * mm;
  G4double core_length = 3.5 * mm; 
  iridium_core = new G4Tubs("iridium_core",corer_min, corer_max/2,core_length/2.,0.*deg,360.*deg);
  logical_iridium_core = new G4LogicalVolume(iridium_core, iridiumMat, "iridium_core_log", 0, 0, 0);
  physical_iridium_core = new G4PVPlacement(0,G4ThreeVector(0,0,0), "phys_iridium_core", logical_iridium_core, physical_air_gap, false, 0, true);

// Visualisations

//Shell/cable attributes    
  steelAttributes = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // red
  steelAttributes -> SetVisibility(true);
  steelAttributes -> SetForceAuxEdgeVisible(true);

  endAttributes = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // red
  endAttributes -> SetVisibility(true);
  endAttributes -> SetForceAuxEdgeVisible(true);
  logical_steel_shell -> SetVisAttributes(steelAttributes);
  logical_End1_steel_shell -> SetVisAttributes(endAttributes);
  logical_End2_steel_shell -> SetVisAttributes(endAttributes);
  logical_cable -> SetVisAttributes(steelAttributes);
 
  G4Colour  magenta (1.0, 0.0, 1.0) ; 

  simpleIridiumVisAtt = new G4VisAttributes(magenta);
  simpleIridiumVisAtt -> SetVisibility(true);
  simpleIridiumVisAtt -> SetForceWireframe(true);
  logical_iridium_core -> SetVisAttributes(simpleIridiumVisAtt);
}

void BrachyDetectorConstructionFlexi::CleanFlexi()
{ 
  delete simpleIridiumVisAtt; 
  simpleIridiumVisAtt = 0;
  
  delete endAttributes; 
  endAttributes = 0;
  
  delete steelAttributes; 
  steelAttributes = 0;
  
  delete physical_iridium_core; 
  physical_iridium_core = 0 ;

  delete logical_iridium_core; 
  logical_iridium_core = 0;
  
  delete iridium_core; 
  iridium_core = 0;
  
  delete physical_cable;
  physical_cable = 0;
 
  delete logical_cable; 
  logical_cable = 0;

  delete cable; 
  cable = 0;
  
  delete physical_End2_steel_shell; 
  physical_End2_steel_shell = 0;
   
  delete logical_End2_steel_shell; 
  logical_End2_steel_shell = 0;
  
  delete End2_steel_shell; 
  End2_steel_shell = 0;
  
  delete physical_End1_steel_shell; 
  physical_End1_steel_shell = 0;
   
  delete logical_End1_steel_shell; 
  logical_End1_steel_shell = 0;
  
  delete End1_steel_shell; 
  End1_steel_shell = 0;

  delete physical_air_gap;
  physical_air_gap = 0;

  delete logical_air_gap;
  logical_air_gap = 0;

  delete air_gap;
  air_gap = 0;

  delete physical_steel_shell;
  physical_steel_shell = 0;

  delete logical_steel_shell;
  logical_steel_shell = 0;
 
  delete steel_shell;
  steel_shell = 0;
 
 G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}
