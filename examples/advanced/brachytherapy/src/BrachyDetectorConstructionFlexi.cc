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
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4NistManager.hh"

BrachyDetectorConstructionFlexi::BrachyDetectorConstructionFlexi()
  : fSteelShell(nullptr), fLogicalSteelShell(nullptr), fAirGap(nullptr), fLogicalAirGap(nullptr), fPhysicalAirGap(nullptr),
    fEnd1SteelShell(nullptr), fLogicalEnd1SteelShell(nullptr), fPhysicalEnd1SteelShell(nullptr),
    fEnd2SteelShell(nullptr), fLogicalEnd2SteelShell(nullptr), fPhysicalEnd2SteelShell(nullptr),
    fCable(nullptr), fLogicalCable(nullptr), fPhysicalCable(nullptr),
    fIridiumCore(nullptr), fLogicalIridiumCore(nullptr), fPhysicalIridiumCore(nullptr),
    fSteelAttributes(nullptr), fEndAttributes(nullptr), fSimpleIridiumVisAtt(nullptr)
{}

void BrachyDetectorConstructionFlexi::ConstructFlexi(G4VPhysicalVolume* mother)
{
  G4NistManager* nist = G4NistManager::Instance();
  
  G4Material* iridiumMat = nist -> FindOrBuildMaterial("G4_Ir");
  G4Material* airMat = nist -> FindOrBuildMaterial("G4_AIR");
  
  //Define Stainless-steel-304 - Flexi source
  G4int Z; //atomic number of the element
  G4Element* elC = nist -> FindOrBuildElement(Z=6);
  G4Element* elMn = nist -> FindOrBuildElement(Z=12);
  G4Element* elSi = nist -> FindOrBuildElement(Z=14);
  G4Element* elCr = nist -> FindOrBuildElement(Z=24);
  G4Element* elFe = nist -> FindOrBuildElement(Z=26);
  G4Element* elNi = nist -> FindOrBuildElement(Z=28);
 
  constexpr G4double d = 7.999*g/cm3;
  auto steelMat = new G4Material("Stainless steel 304",d,6);
  steelMat -> AddElement(elMn, 0.02);
  steelMat -> AddElement(elSi, 0.01);
  steelMat -> AddElement(elCr, 0.19);
  steelMat -> AddElement(elNi, 0.10);
  steelMat -> AddElement(elFe, 0.6792);
  steelMat -> AddElement(elC, 0.0008);
  
 //Define dimensions of the outer Steel shell around the solid source - not including the ends 

  G4double shellr_min = 0.00 * mm;
  G4double shellr_max = 0.85 * mm;
  G4double shell_length = 3.6 * mm; 
  fSteelShell = new G4Tubs("steel_shell",shellr_min, shellr_max/2, shell_length/2.,0.*deg,360.*deg);
  fLogicalSteelShell = new G4LogicalVolume(fSteelShell, steelMat, "steel_shell_log", nullptr, nullptr, nullptr);
  fPhysicalSteelShell = new G4PVPlacement(nullptr,G4ThreeVector(0,0,0),"phys_steel_shell", fLogicalSteelShell, mother, false, 0, true);

//Define dimensions of the air gap between Steel shell and Iridium core
  G4double airr_min = 0.00 * mm;
  G4double airr_max = 0.67 * mm;
  G4double air_length = 3.6 * mm; 
  fAirGap = new G4Tubs("air_gap", airr_min, airr_max/2, air_length/2, 0.*deg, 360.*deg);
  fLogicalAirGap = new G4LogicalVolume(fAirGap, airMat, "air_gap_log", nullptr, nullptr, nullptr);
  fPhysicalAirGap = new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), "phys_air_gap", fLogicalAirGap, fPhysicalSteelShell, false, 0, true);

//Define the non-cable weld end of the Steel shell 
  G4double end1r_min = 0.0 * mm;
  G4double end1r_max = 0.85 * mm;
  G4double end1length = 0.65 * mm; 
  fEnd1SteelShell = new G4Tubs("End_1_steel_shell", end1r_min, end1r_max/2, end1length/2.,0.*deg,360.*deg);
  fLogicalEnd1SteelShell = new G4LogicalVolume(fEnd1SteelShell, steelMat, "End1_steel_shell_log", nullptr, nullptr, nullptr);
  G4double end1offset_x = 0.0 * mm;
  G4double end1offset_y = 0.0 * mm;
  G4double end1offset_z = 2.125 * mm;
  fPhysicalEnd1SteelShell = new G4PVPlacement(nullptr,G4ThreeVector(end1offset_x,end1offset_y,end1offset_z),"phys_End1_steel_shell", fLogicalEnd1SteelShell,mother, false, 0, true);

//Define the cable weld end of the Steel shell 
  G4double end2r_min1 = 0.0 * mm;
  G4double end2r_max1 = 0.85 * mm;
  G4double end2r_min2 = 0.0 * mm;
  G4double end2r_max2 = 0.5 * mm;
  G4double end2length = 0.4 * mm;
  fEnd2SteelShell = new G4Cons("End_2_steel_shell", end2r_min2, end2r_max2/2, end2r_min1, end2r_max1/2, end2length/2.0, 0.0, 360.0*deg);
  fLogicalEnd2SteelShell = new G4LogicalVolume(fEnd2SteelShell, steelMat, "End2_steel_shell_log", nullptr, nullptr, nullptr);
  G4double end2offset_x = 0.0 * mm;
  G4double end2offset_y = 0.0 * mm;
  G4double end2offset_z = -2.0 * mm; 
  fPhysicalEnd2SteelShell = new G4PVPlacement(nullptr,G4ThreeVector(end2offset_x,end2offset_y,end2offset_z), "phys_End2_steel_shell", fLogicalEnd2SteelShell, mother, false, 0, true);

//Define the cable 
  G4double cable_min = 0.0 * mm;
  G4double cable_max = 0.5 * mm;
  G4double cablelength = 5.0 * mm; 
  fCable = new G4Tubs("cable",cable_min, cable_max/2, cablelength/2.,0.*deg,360.*deg);
  fLogicalCable = new G4LogicalVolume(fCable, steelMat, "cable_log", nullptr, nullptr, nullptr);
  G4double cableoffset_x = 0.0 * mm;
  G4double cableoffset_y = 0.0 * mm;
  G4double cableoffset_z = -4.7 * mm;
  fPhysicalCable = new G4PVPlacement(nullptr,G4ThreeVector(cableoffset_x,cableoffset_y,cableoffset_z),"phys_cable", fLogicalCable, mother, false, 0, true);

// Define the Iridium core
  G4double corer_min = 0.0 * mm;	
  G4double corer_max = 0.6 * mm;
  G4double core_length = 3.5 * mm; 
  fIridiumCore = new G4Tubs("iridium_core",corer_min, corer_max/2,core_length/2.,0.*deg,360.*deg);
  fLogicalIridiumCore = new G4LogicalVolume(fIridiumCore, iridiumMat, "iridium_core_log", nullptr, nullptr, nullptr);
  fPhysicalIridiumCore = new G4PVPlacement(nullptr,G4ThreeVector(0,0,0), "phys_iridium_core", fLogicalIridiumCore, fPhysicalAirGap, false, 0, true);

// Visualisations

//Shell/cable attributes    
  fSteelAttributes = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // red
  fSteelAttributes -> SetVisibility(true);
  fSteelAttributes -> SetForceAuxEdgeVisible(true);

  fEndAttributes = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // red
  fEndAttributes -> SetVisibility(true);
  fEndAttributes -> SetForceAuxEdgeVisible(true);
  fLogicalSteelShell -> SetVisAttributes(fSteelAttributes);
  fLogicalEnd1SteelShell -> SetVisAttributes(fEndAttributes);
  fLogicalEnd2SteelShell -> SetVisAttributes(fEndAttributes);
  fLogicalCable -> SetVisAttributes(fSteelAttributes);
 
  G4Colour  magenta (1.0, 0.0, 1.0) ; 

  fSimpleIridiumVisAtt = new G4VisAttributes(magenta);
  fSimpleIridiumVisAtt -> SetVisibility(true);
  fSimpleIridiumVisAtt -> SetForceWireframe(true);
  fLogicalIridiumCore -> SetVisAttributes(fSimpleIridiumVisAtt);
}

void BrachyDetectorConstructionFlexi::CleanFlexi()
{ 
  delete fSimpleIridiumVisAtt; 
  fSimpleIridiumVisAtt = nullptr;
  
  delete fEndAttributes; 
  fEndAttributes = nullptr;
  
  delete fSteelAttributes; 
  fSteelAttributes = nullptr;
  
  delete fPhysicalIridiumCore; 
  fPhysicalIridiumCore = nullptr ;

  delete fLogicalIridiumCore; 
  fLogicalIridiumCore = nullptr;
  
  delete fIridiumCore; 
  fIridiumCore = nullptr;
  
  delete fPhysicalCable;
  fPhysicalCable = nullptr;
 
  delete fLogicalCable; 
  fLogicalCable = nullptr;

  delete fCable; 
  fCable = nullptr;
  
  delete fPhysicalEnd2SteelShell; 
  fPhysicalEnd2SteelShell = nullptr;
   
  delete fLogicalEnd2SteelShell; 
  fLogicalEnd2SteelShell = nullptr;
  
  delete fEnd2SteelShell; 
  fEnd2SteelShell = nullptr;
  
  delete fPhysicalEnd1SteelShell; 
  fPhysicalEnd1SteelShell = nullptr;
   
  delete fLogicalEnd1SteelShell; 
  fLogicalEnd1SteelShell = nullptr;
  
  delete fEnd1SteelShell; 
  fEnd1SteelShell = nullptr;

  delete fPhysicalAirGap;
  fPhysicalAirGap = nullptr;

  delete fLogicalAirGap;
  fLogicalAirGap = nullptr;

  delete fAirGap;
  fAirGap = nullptr;

  delete fPhysicalSteelShell;
  fPhysicalSteelShell = nullptr;

  delete fLogicalSteelShell;
  fLogicalSteelShell = nullptr;
 
  delete fSteelShell;
  fSteelShell = nullptr;
 
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}
