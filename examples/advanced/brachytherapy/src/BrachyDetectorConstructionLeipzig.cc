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
// S.Guatelli
//
//    *******************************************
//    *                                         *
//    *    BrachyDetectorConstructionLeipzig.cc *
//    *                                         *
//    *******************************************
//
//
// $Id$
//

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "BrachyDetectorConstructionLeipzig.hh"
#include "G4CSGSolid.hh"
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

// Leipzig Applicator ...

BrachyDetectorConstructionLeipzig::BrachyDetectorConstructionLeipzig()
  : 
  capsulePhys(0),
  capsuleTipPhys(0),
  iridiumCorePhys(0),
  applicator1Phys(0),
  applicator2Phys(0)
{ 
  pMaterial = new BrachyMaterial();
}

BrachyDetectorConstructionLeipzig::~BrachyDetectorConstructionLeipzig()
{ 
  delete pMaterial; 
}

void  BrachyDetectorConstructionLeipzig::ConstructLeipzig(G4VPhysicalVolume*   mother)
{
  G4Colour  red     (1.0, 0.0, 0.0) ; 
  G4Colour  lblue   (0.0, 0.0, .75);

  G4Material* capsuleMat = pMaterial -> GetMat("Stainless steel");
  G4Material* iridium = pMaterial -> GetMat("Iridium");
  G4Material* tungsten =pMaterial -> GetMat("Tungsten");

  //Iridium source ...

  G4Tubs* capsule = new G4Tubs("Capsule",0,0.55*mm,3.725*mm,0.*deg,360.*deg);
  G4LogicalVolume* capsuleLog = new G4LogicalVolume(capsule,capsuleMat,"CapsuleLog");
  capsulePhys = new G4PVPlacement(0,
                                  G4ThreeVector(0,0,-1.975*mm),
                                  "CapsulePhys",
				  capsuleLog,
                                  mother, //mother volume: phantom
                                  false,
                                  0);

  // Capsule tip
  G4Sphere* capsuleTip = new G4Sphere("CapsuleTip",0.*mm,0.55*mm,0.*deg,360.*deg,0.*deg,90.*deg);
  G4LogicalVolume* capsuleTipLog = new G4LogicalVolume(capsuleTip,capsuleMat,"CapsuleTipLog");
  capsuleTipPhys = new G4PVPlacement(0,
                                     G4ThreeVector(0.,0.,1.75*mm),
                                     "CapsuleTipPhys",
				     capsuleTipLog,
                                     mother,
                                     false,
                                     0);
  // Iridium core

  G4Tubs* iridiumCore = new G4Tubs("IrCore",0,0.30*mm,1.75*mm,0.*deg,360.*deg);
  G4LogicalVolume* iridiumCoreLog = new G4LogicalVolume(iridiumCore,
                                       iridium,
                                       "IridiumCoreLog");
  iridiumCorePhys = new G4PVPlacement(0,
                                      G4ThreeVector(0.,0.,1.975*mm),
                                      "IridiumCorePhys",
				      iridiumCoreLog,
                                      capsulePhys,
                                      false,
                                      0);

  //Leipzig Applicator is modelled with two different volumes

  G4Tubs* applicator1 = new G4Tubs("Appl1",5*mm,10.5*mm,12*mm,0.*deg,360.*deg);
  G4LogicalVolume* applicator1Log = new G4LogicalVolume(applicator1,tungsten,"Appl1Log");
  applicator1Phys = new G4PVPlacement(0,
                                      G4ThreeVector(0,0,4.0*mm),
                                      "Appl1Phys", 
                                      applicator1Log,
				      mother,
                                      false,
                                      0);

  G4Tubs* applicator2 = new G4Tubs("Appl2",0.55*mm,5.*mm,3.125*mm,0.*deg,360.*deg);
  G4LogicalVolume* applicator2Log = new G4LogicalVolume(applicator2,tungsten,"Appl2");
  applicator2Phys = new G4PVPlacement(0,
                                      G4ThreeVector(0,0,-4.875*mm),
                                      "Appl2Phys",
                                      applicator2Log,
				      mother,
                                      false,
                                      0);

  G4VisAttributes* simpleCapsuleVisAtt = new G4VisAttributes(red); 
  simpleCapsuleVisAtt -> SetVisibility(true); 
  simpleCapsuleVisAtt -> SetForceSolid(true); 
  capsuleLog -> SetVisAttributes(simpleCapsuleVisAtt); 

  G4VisAttributes* simpleCapsuleTipVisAtt = new G4VisAttributes(red); 
  simpleCapsuleTipVisAtt -> SetVisibility(true); 
  simpleCapsuleTipVisAtt -> SetForceSolid(true); 
  capsuleTipLog -> SetVisAttributes(simpleCapsuleTipVisAtt);
 
  G4VisAttributes* applicatorVisAtt = new G4VisAttributes(lblue);
  applicatorVisAtt -> SetVisibility(true);
  applicatorVisAtt -> SetForceWireframe(true);
  applicator1Log ->  SetVisAttributes(applicatorVisAtt);  
  applicator2Log ->  SetVisAttributes(applicatorVisAtt);  
}
