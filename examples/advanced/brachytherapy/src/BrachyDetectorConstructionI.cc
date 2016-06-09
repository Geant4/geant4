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
// $Id: BrachyDetectorConstructionI.cc,v 1.9 2006-06-29 15:48:11 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "globals.hh"
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
  capsulePhys(0),capsuleTipPhys1(0),capsuleTipPhys2(0), iodiumCorePhys(0),markerPhys(0)
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
  G4Material* iodium = pMaterial -> GetMat("Iodium");
  G4Material* gold = pMaterial -> GetMat("gold");

  G4Colour  red     (1.0, 0.0, 0.0) ;
  G4Colour  magenta (1.0, 0.0, 1.0) ; 
  G4Colour  lblue   (0.0, 0.0, .75);

  // Air tub
  G4Tubs* defaultTub = new G4Tubs("DefaultTub",0.*mm, 0.40*mm, 1.84*mm, 0.*deg, 360.*deg);
  G4LogicalVolume* defaultTubLog = new G4LogicalVolume(defaultTub,air,"DefaultTub_Log");
  G4VPhysicalVolume* defaultTubPhys = new G4PVPlacement(0,
                                      G4ThreeVector(),
                                      "defaultTub_Phys",
                                      defaultTubLog,
                                      mother,
                                      false,
                                      0); 
  //  Capsule main body ...
  G4double capsuleR = 0.35*mm;
  G4Tubs* capsule = new G4Tubs("Capsule", capsuleR,0.40*mm,1.84*mm,0.*deg,360.*deg);
  G4LogicalVolume* capsuleLog = new G4LogicalVolume(capsule,titanium,"CapsuleLog");
  capsulePhys = new G4PVPlacement(0,
                                  G4ThreeVector(),
                                  "CapsulePhys",
                                  capsuleLog,
                                  defaultTubPhys,
                                  false,
                                  0);
  // Capsule tips
  G4Sphere* capsuleTip = new G4Sphere("CapsuleTip",
                            0.*mm,
                            0.40*mm,
                            0.*deg,
                            360.*deg,
                            0.*deg,
                            90.*deg);
  G4LogicalVolume* capsuleTipLog = new G4LogicalVolume(capsuleTip,titanium,"CapsuleTipLog");
  capsuleTipPhys1 = new G4PVPlacement(0,
                                      G4ThreeVector(0.,0.,1.84*mm),
                                      "CapsuleTipPhys1",
                                      capsuleTipLog,
                                      mother,
                                      false,
                                      0);

  G4RotationMatrix* rotateMatrix = new G4RotationMatrix();
  rotateMatrix -> rotateX(180.0*deg);
  capsuleTipPhys2 = new G4PVPlacement(rotateMatrix, 
                                      G4ThreeVector(0,0,-1.84*mm),
                                      "CapsuleTipPhys2",
                                      capsuleTipLog,
                                      mother,
                                      false,
                                      0);
 
  // Radiactive core ...
  G4Tubs* iodiumCore = new G4Tubs("ICore",0.085*mm,0.35*mm,1.75*mm,0.*deg,360.*deg);
  G4LogicalVolume* iodiumCoreLog = new G4LogicalVolume(iodiumCore,iodium,"iodiumCoreLog");
  iodiumCorePhys = new G4PVPlacement(0,
                                     G4ThreeVector(0.,0.,0.),
                                     "iodiumCorePhys",
                                     iodiumCoreLog,
                                     defaultTubPhys,
                                     false,
                                     0);
  // Golden marker
  G4Tubs* marker = new G4Tubs("GoldenMarker",0.*mm,0.085*mm,1.75*mm,0.*deg,360.*deg);
  G4LogicalVolume* markerLog = new G4LogicalVolume(marker,gold,"MarkerLog");
  markerPhys = new G4PVPlacement(0,
                                 G4ThreeVector(0.,0.,0.),
                                 "MarkerPhys",
                                 markerLog,
                                 defaultTubPhys,
                                 false,
                                 0);
  // Visual attributes ...
  G4VisAttributes* simpleMarkerVisAtt= new G4VisAttributes(lblue);
  simpleMarkerVisAtt -> SetVisibility(true);
  simpleMarkerVisAtt -> SetForceSolid(true);
  markerLog -> SetVisAttributes( simpleMarkerVisAtt);  

  G4VisAttributes* simpleiodiumVisAtt= new G4VisAttributes(magenta);
  simpleiodiumVisAtt -> SetVisibility(true);
  simpleiodiumVisAtt -> SetForceWireframe(true);
  iodiumCoreLog -> SetVisAttributes(simpleiodiumVisAtt);

  G4VisAttributes* simpleCapsuleVisAtt= new G4VisAttributes(red);
  simpleCapsuleVisAtt -> SetVisibility(true);  
  simpleCapsuleVisAtt -> SetForceWireframe(true);
  capsuleLog -> SetVisAttributes( simpleCapsuleVisAtt);

  G4VisAttributes* simpleCapsuleTipVisAtt= new G4VisAttributes(red);
  simpleCapsuleTipVisAtt -> SetVisibility(true); 
  simpleCapsuleTipVisAtt -> SetForceSolid(true);
  capsuleTipLog -> SetVisAttributes( simpleCapsuleTipVisAtt);
}
