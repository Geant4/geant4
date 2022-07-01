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
// D. Cutajar, A. Le
//
//    ********************************************
//    *                                          *
//    *  BrachyDetectorConstructionOncura6711.cc *
//    *                                          *
//    ********************************************
//
// 
//
#include "BrachyDetectorConstructionOncura6711.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Sphere.hh"
#include "G4RunManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4TransportationManager.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

BrachyDetectorConstructionOncura6711::BrachyDetectorConstructionOncura6711()
  : 
   fOncuraCapsule(nullptr), fOncuraCapsuleLog(nullptr),
   fOncuraCapsulePhys(nullptr),
   fOncuraCapsuleTip1(nullptr), fOncuraCapsuleTip1Log(nullptr),
   fOncuraCapsuleTip1Phys(nullptr),
   fOncuraCapsuleTip2(nullptr), fOncuraCapsuleTip2Log(nullptr),
   fOncuraCapsuleTip2Phys(nullptr),
   fOncuraAirGap(nullptr), fOncuraAirGapLog(nullptr),
   fOncuraAirGapPhys(nullptr),
   fOncuraSilverCore(nullptr), fOncuraSilverCoreLog(nullptr),
   fOncuraSilverCorePhys(nullptr),
   fOncuraCapsuleShellVisAtt(nullptr), fOncuraCapsuleTipVisAtt(nullptr),
   fOncuraSilverCoreVisAtt(nullptr)
{}

void BrachyDetectorConstructionOncura6711::ConstructOncura6711(G4VPhysicalVolume* mother)
{
G4Colour  red     (1.0, 0.0, 0.0) ;
G4Colour  magenta (1.0, 0.0, 1.0) ; 

G4NistManager* nist = G4NistManager::Instance();
G4Material* titanium = nist -> FindOrBuildMaterial("G4_Ti");
G4Material* air = nist -> FindOrBuildMaterial("G4_AIR");
G4Material* silver = nist -> FindOrBuildMaterial("G4_Ag");

//Capsule shell
fOncuraCapsule = new G4Tubs("OncuraCapsule",0,0.4*mm,1.875*mm,0.*deg,360.*deg);
fOncuraCapsuleLog = new G4LogicalVolume(fOncuraCapsule,titanium,"OncuraCapsuleLog", nullptr,nullptr,nullptr);
fOncuraCapsulePhys = new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), 
  								"OncuraCapsulePhys", fOncuraCapsuleLog, 
  								mother, false, 0, true);
						
//Capsule tips
fOncuraCapsuleTip1 = new G4Sphere("OncuraCapsuleTip1", 0, 0.4*mm, 0., 360*deg, 0., 90*deg);
fOncuraCapsuleTip1Log = new G4LogicalVolume(fOncuraCapsuleTip1, titanium, "OncuraCapsuleTip1Log",nullptr,nullptr,nullptr);
fOncuraCapsuleTip1Phys = new G4PVPlacement(nullptr, G4ThreeVector(0,0,1.875*mm), 
								"OncuraCapsuleTip1Phys", fOncuraCapsuleTip1Log,
								mother, false, 
								0, true);

fOncuraCapsuleTip2 = new G4Sphere("OncuraCapsuleTip2", 0, 0.4*mm, 0., 360*deg, 90*deg, 90*deg);
fOncuraCapsuleTip2Log = new G4LogicalVolume(fOncuraCapsuleTip2, titanium, "OncuraCapsuleTip2Log", nullptr,nullptr,nullptr);
fOncuraCapsuleTip2Phys = new G4PVPlacement(nullptr, G4ThreeVector(0,0,-1.875*mm), 
								"OncuraCapsuleTip2Phys", fOncuraCapsuleTip2Log,
								mother, false, 
								0, true);

//Air gap
fOncuraAirGap = new G4Tubs("OncuraAirGap",0,0.33*mm,1.825*mm,0.*deg,360.*deg);
fOncuraAirGapLog = new G4LogicalVolume(fOncuraAirGap, air, "OncuraAirGapLog");
fOncuraAirGapPhys = new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), 
                                "OncuraAirGapPhys", fOncuraAirGapLog, 
                                fOncuraCapsulePhys, false,
                                0, true);

//Silver core
fOncuraSilverCore = new G4Tubs("OncuraSilverCore",0,0.25*mm,1.4*mm,0.*deg,360.*deg);
fOncuraSilverCoreLog = new G4LogicalVolume(fOncuraSilverCore, silver, "silverCoreLog");
fOncuraSilverCorePhys = new G4PVPlacement(nullptr, G4ThreeVector(0,0,0), 
								"OncuraSilverCorePhys", fOncuraSilverCoreLog,
                                fOncuraAirGapPhys, false,
                                0, true);

fOncuraCapsuleShellVisAtt = new G4VisAttributes(red);
fOncuraCapsuleShellVisAtt -> SetVisibility(true);  
fOncuraCapsuleShellVisAtt -> SetForceWireframe(true);
fOncuraCapsuleLog -> SetVisAttributes(fOncuraCapsuleShellVisAtt);

fOncuraCapsuleTipVisAtt = new G4VisAttributes(red);
fOncuraCapsuleTipVisAtt -> SetVisibility(true);  
fOncuraCapsuleTipVisAtt -> SetForceSolid(true);
fOncuraCapsuleTip1Log -> SetVisAttributes(fOncuraCapsuleTipVisAtt);
fOncuraCapsuleTip2Log -> SetVisAttributes(fOncuraCapsuleTipVisAtt);

fOncuraSilverCoreVisAtt = new G4VisAttributes(magenta);
fOncuraSilverCoreVisAtt -> SetVisibility(true);
fOncuraSilverCoreVisAtt -> SetForceSolid(true);
fOncuraSilverCoreLog -> SetVisAttributes(fOncuraSilverCoreVisAtt);
}

void BrachyDetectorConstructionOncura6711::CleanOncura6711()
{ 
 delete fOncuraSilverCoreVisAtt;
 fOncuraSilverCoreVisAtt = nullptr;

 delete fOncuraCapsuleTipVisAtt;
 fOncuraCapsuleTipVisAtt = nullptr;

 delete fOncuraCapsuleShellVisAtt;
 fOncuraCapsuleShellVisAtt = nullptr;

 delete fOncuraSilverCorePhys;
 fOncuraSilverCorePhys = nullptr;

 delete fOncuraSilverCoreLog;
 fOncuraSilverCoreLog = nullptr;

 delete fOncuraSilverCore;
 fOncuraSilverCore = nullptr;

 delete fOncuraAirGapPhys;
 fOncuraAirGapPhys = nullptr;

 delete fOncuraAirGapLog;
 fOncuraAirGapLog = nullptr;

 delete fOncuraAirGap;
 fOncuraAirGap = nullptr;

 delete fOncuraCapsuleTip2Phys;
 fOncuraCapsuleTip2Phys = nullptr;

 delete fOncuraCapsuleTip2Log;
 fOncuraCapsuleTip2Log = nullptr;
 
 delete fOncuraCapsuleTip2;
 fOncuraCapsuleTip2 = nullptr;

 delete fOncuraCapsuleTip1Phys;
 fOncuraCapsuleTip1Phys = nullptr;

 delete fOncuraCapsuleTip1Log;
 fOncuraCapsuleTip1Log = nullptr;

 delete fOncuraCapsuleTip1;
 fOncuraCapsuleTip1 = nullptr;

 delete fOncuraCapsulePhys;
 fOncuraCapsulePhys = nullptr;
	
 delete fOncuraCapsuleLog;
 fOncuraCapsuleLog = nullptr;

 delete fOncuraCapsule;
 fOncuraCapsule = nullptr;
  
 G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}
