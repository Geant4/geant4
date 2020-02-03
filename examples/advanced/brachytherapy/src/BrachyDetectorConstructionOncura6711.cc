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
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "BrachyDetectorConstructionOncura6711.hh"
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
#include "BrachyMaterial.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

BrachyDetectorConstructionOncura6711::BrachyDetectorConstructionOncura6711()
  : 
   OncuraCapsule(0),OncuraCapsuleLog(0),
   OncuraCapsulePhys(0),
   OncuraCapsuleTip1(0),OncuraCapsuleTip1Log(0),
   OncuraCapsuleTip1Phys(0),
   OncuraCapsuleTip2(0),OncuraCapsuleTip2Log(0),
   OncuraCapsuleTip2Phys(0),
   OncuraAirGap(0),OncuraAirGapLog(0),
   OncuraAirGapPhys(0),
   OncuraSilverCore(0), OncuraSilverCoreLog(0),
   OncuraSilverCorePhys(0),
   OncuraCapsuleShellVisAtt(0),OncuraCapsuleTipVisAtt(0),
   OncuraSilverCoreVisAtt(0)
{
  pMat = new BrachyMaterial();
}

BrachyDetectorConstructionOncura6711::~BrachyDetectorConstructionOncura6711()
{ 
  delete pMat; 
}

void BrachyDetectorConstructionOncura6711::ConstructOncura6711(G4VPhysicalVolume* mother)
{
G4Colour  red     (1.0, 0.0, 0.0) ;
G4Colour  magenta (1.0, 0.0, 1.0) ; 

G4Material* titaniumMat = pMat -> GetMat("titanium");
G4Material* airMat = pMat -> GetMat("Air");
G4Material* silverMat = pMat -> GetMat("Silver");

//Capsule shell
OncuraCapsule = new G4Tubs("OncuraCapsule",0,0.4*mm,1.875*mm,0.*deg,360.*deg);
OncuraCapsuleLog = new G4LogicalVolume(OncuraCapsule,titaniumMat,"OncuraCapsuleLog", 0,0,0);
OncuraCapsulePhys = new G4PVPlacement(0, G4ThreeVector(0,0,0), 
  								"OncuraCapsulePhys", OncuraCapsuleLog, 
  								mother, false,
                                0, true);
						
//Capsule tips
OncuraCapsuleTip1 = new G4Sphere("OncuraCapsuleTip1", 0, 0.4*mm, 0., 360*deg, 0., 90*deg);
OncuraCapsuleTip1Log = new G4LogicalVolume(OncuraCapsuleTip1, titaniumMat, "OncuraCapsuleTip1Log",0,0,0);
OncuraCapsuleTip1Phys = new G4PVPlacement(0, G4ThreeVector(0,0,1.875*mm), 
								"OncuraCapsuleTip1Phys", OncuraCapsuleTip1Log,
								mother, false, 
								0, true);

OncuraCapsuleTip2 = new G4Sphere("OncuraCapsuleTip2", 0, 0.4*mm, 0., 360*deg, 90*deg, 90*deg);
OncuraCapsuleTip2Log = new G4LogicalVolume(OncuraCapsuleTip2, titaniumMat, "OncuraCapsuleTip2Log",0,0,0);
OncuraCapsuleTip2Phys = new G4PVPlacement(0, G4ThreeVector(0,0,-1.875*mm), 
								"OncuraCapsuleTip2Phys", OncuraCapsuleTip2Log,
								mother, false, 
								0, true);

//Airgap
OncuraAirGap = new G4Tubs("OncuraAirGap",0,0.33*mm,1.825*mm,0.*deg,360.*deg);
OncuraAirGapLog = new G4LogicalVolume(OncuraAirGap, airMat, "OncuraAirGapLog");
OncuraAirGapPhys = new G4PVPlacement(0, G4ThreeVector(0,0,0), 
                                "OncuraAirGapPhys", OncuraAirGapLog, 
                                OncuraCapsulePhys, false,
                                0, true);

//Silver core
OncuraSilverCore = new G4Tubs("OncuraSilverCore",0,0.25*mm,1.4*mm,0.*deg,360.*deg);
OncuraSilverCoreLog = new G4LogicalVolume(OncuraSilverCore, silverMat, "silverCoreLog");
OncuraSilverCorePhys = new G4PVPlacement(0, G4ThreeVector(0,0,0), 
								"OncuraSilverCorePhys", OncuraSilverCoreLog,
                                OncuraAirGapPhys, false,
                                0, true);

OncuraCapsuleShellVisAtt = new G4VisAttributes(red);
OncuraCapsuleShellVisAtt -> SetVisibility(true);  
OncuraCapsuleShellVisAtt -> SetForceWireframe(true);
OncuraCapsuleLog -> SetVisAttributes(OncuraCapsuleShellVisAtt);

OncuraCapsuleTipVisAtt = new G4VisAttributes(red);
OncuraCapsuleTipVisAtt -> SetVisibility(true);  
OncuraCapsuleTipVisAtt -> SetForceSolid(true);
OncuraCapsuleTip1Log -> SetVisAttributes(OncuraCapsuleTipVisAtt);
OncuraCapsuleTip2Log -> SetVisAttributes(OncuraCapsuleTipVisAtt);

OncuraSilverCoreVisAtt = new G4VisAttributes(magenta);
OncuraSilverCoreVisAtt -> SetVisibility(true);
OncuraSilverCoreVisAtt -> SetForceSolid(true);
OncuraSilverCoreLog -> SetVisAttributes(OncuraSilverCoreVisAtt);

}

void BrachyDetectorConstructionOncura6711::CleanOncura6711()
{ 
 delete OncuraSilverCoreVisAtt;
 OncuraSilverCoreVisAtt = 0;

 delete OncuraCapsuleTipVisAtt;
 OncuraCapsuleTipVisAtt = 0;

 delete OncuraCapsuleShellVisAtt;
 OncuraCapsuleShellVisAtt = 0;

 delete OncuraSilverCorePhys;
 OncuraSilverCorePhys = 0;

 delete OncuraSilverCoreLog;
 OncuraSilverCoreLog = 0;

 delete OncuraSilverCore;
 OncuraSilverCore = 0;

 delete OncuraAirGapPhys;
 OncuraAirGapPhys = 0;

 delete OncuraAirGapLog;
 OncuraAirGapLog = 0;

 delete OncuraAirGap;
 OncuraAirGap = 0;

 delete OncuraCapsuleTip2Phys;
 OncuraCapsuleTip2Phys = 0;

 delete OncuraCapsuleTip2Log;
 OncuraCapsuleTip2Log = 0;
 
 delete OncuraCapsuleTip2;
 OncuraCapsuleTip2 = 0;

 delete OncuraCapsuleTip1Phys;
 OncuraCapsuleTip1Phys = 0;

 delete OncuraCapsuleTip1Log;
 OncuraCapsuleTip1Log = 0;

 delete OncuraCapsuleTip1;
 OncuraCapsuleTip1 = 0;

 delete OncuraCapsulePhys;
 OncuraCapsulePhys = 0;
	
 delete OncuraCapsuleLog;
 OncuraCapsuleLog = 0;

 delete OncuraCapsule;
 OncuraCapsule = 0;
  
 G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}
