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
// Code developed by: S.Guatelli
//
///    ****************************************
//    *                                      *
//    *    BrachyDetectorConstructionI.cc     *
//    *                                      *
//    ****************************************
//
//
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "BrachyDetectorConstructionI.hh"
#include "G4CSGSolid.hh"
#include "G4Sphere.hh"
#include "G4MaterialPropertyVector.hh"
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

BrachyDetectorConstructionI::BrachyDetectorConstructionI():
defaultTub(0), capsule(0), capsuleTip(0),iodiumCore(0), defaultTubLog(0),
capsuleLog(0), capsuleTipLog(0), iodiumCoreLog(0),defaultTubPhys(0), 
capsulePhys(0),capsuleTipPhys1(0),capsuleTipPhys2(0), iodiumCorePhys(0),
simpleiodiumVisAtt(0), simpleCapsuleVisAtt(0), simpleCapsuleTipVisAtt(0)
{ 
  pMaterial = new BrachyMaterial();
}

BrachyDetectorConstructionI::~BrachyDetectorConstructionI()
{ 
  delete pMaterial; 
}

void BrachyDetectorConstructionI::ConstructIodium(G4VPhysicalVolume* mother)
{
  // source Bebig Isoseed I-125 ...

  //Get materials for source construction ...
  G4Material* titanium = pMaterial -> GetMat("titanium");
  G4Material* air = pMaterial -> GetMat("Air");
  G4Material* iodium = pMaterial -> GetMat("Iodine");
 
  G4Colour  red     (1.0, 0.0, 0.0) ;
  G4Colour  magenta (1.0, 0.0, 1.0) ; 
  G4Colour  lblue   (0.0, 0.0, .75);

  // Air tub
  defaultTub = new G4Tubs("DefaultTub",0.*mm, 0.40*mm, 1.84*mm, 0.*deg, 360.*deg);
  defaultTubLog = new G4LogicalVolume(defaultTub,air,"DefaultTub_Log");
  defaultTubPhys = new G4PVPlacement(0,
                                      G4ThreeVector(),
                                      "defaultTub_Phys",
                                      defaultTubLog,
                                      mother,
                                      false,
                                      0, true); 
  //  Capsule main body ...
  G4double capsuleR = 0.35*mm;
  capsule = new G4Tubs("Capsule", capsuleR,0.40*mm,1.84*mm,0.*deg,360.*deg);
  capsuleLog = new G4LogicalVolume(capsule,titanium,"CapsuleLog");
  capsulePhys = new G4PVPlacement(0,
                                  G4ThreeVector(),
                                  "CapsulePhys",
                                  capsuleLog,
                                  defaultTubPhys,
                                  false,
                                  0, true);
  // Capsule tips
  capsuleTip = new G4Sphere("CapsuleTip",
                            0.*mm,
                            0.40*mm,
                            0.*deg,
                            360.*deg,
                            0.*deg,
                            90.*deg);
  capsuleTipLog = new G4LogicalVolume(capsuleTip,titanium,"CapsuleTipLog");
  capsuleTipPhys1 = new G4PVPlacement(0,
                                      G4ThreeVector(0.,0.,1.84*mm),
                                      "IodineCapsuleTipPhys1",
                                      capsuleTipLog,
                                      mother,
                                      false,
                                      0, true);

  G4RotationMatrix* rotateMatrix = new G4RotationMatrix();
  rotateMatrix -> rotateX(180.0*deg);
  capsuleTipPhys2 = new G4PVPlacement(rotateMatrix, 
                                      G4ThreeVector(0,0,-1.84*mm),
                                      "IodineCapsuleTipPhys2",
                                      capsuleTipLog,
                                      mother,
                                      false,
                                      0, true);
 
  // Radiactive core ...
  iodiumCore = new G4Tubs("ICore",0.085*mm,0.35*mm,1.75*mm,0.*deg,360.*deg);
  iodiumCoreLog = new G4LogicalVolume(iodiumCore,iodium,"iodiumCoreLog");
  iodiumCorePhys = new G4PVPlacement(0,
                                     G4ThreeVector(0.,0.,0.),
                                     "iodiumCorePhys",
                                     iodiumCoreLog,
                                     defaultTubPhys,
                                     false,
                                     0, true);
 
  // Visual attributes ...
  
  simpleiodiumVisAtt= new G4VisAttributes(magenta);
  simpleiodiumVisAtt -> SetVisibility(true);
  simpleiodiumVisAtt -> SetForceSolid(true);
  iodiumCoreLog -> SetVisAttributes(simpleiodiumVisAtt);

  simpleCapsuleVisAtt= new G4VisAttributes(red);
  simpleCapsuleVisAtt -> SetVisibility(true);  
  simpleCapsuleVisAtt -> SetForceWireframe(true);
  capsuleLog -> SetVisAttributes( simpleCapsuleVisAtt);

  simpleCapsuleTipVisAtt= new G4VisAttributes(red);
  simpleCapsuleTipVisAtt -> SetVisibility(true); 
  simpleCapsuleTipVisAtt -> SetForceSolid(true);
  capsuleTipLog -> SetVisAttributes( simpleCapsuleTipVisAtt);
}

void BrachyDetectorConstructionI::CleanIodium()
{
delete simpleiodiumVisAtt; simpleiodiumVisAtt = 0;
delete simpleCapsuleVisAtt; simpleCapsuleVisAtt = 0;
delete simpleCapsuleTipVisAtt; simpleCapsuleTipVisAtt = 0;
delete capsuleTipPhys1; capsuleTipPhys1 = 0; 
delete capsuleTipPhys2; capsuleTipPhys2 = 0;
delete iodiumCorePhys; iodiumCorePhys = 0;
delete capsulePhys; capsulePhys = 0;
delete defaultTubPhys; defaultTubPhys = 0; 
delete defaultTubLog; defaultTubLog = 0;
delete capsuleLog; capsuleLog = 0;
delete capsuleTipLog; capsuleTipLog = 0;
delete iodiumCoreLog; iodiumCoreLog = 0;
delete defaultTub; defaultTub = 0;
delete capsule; capsule = 0;
delete capsuleTip; capsuleTip = 0;
delete iodiumCore; iodiumCore = 0;

G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}
