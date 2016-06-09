//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// --------------------------------------------------------------
//                 GEANT 4 - Brachytherapy example
// --------------------------------------------------------------
//
// Code developed by: S.Guatelli
//
//    *******************************************
//    *                                         *
//    *    BrachyDetectorConstructionLeipzig.cc *
//    *                                         *
//    *******************************************
//
//
// $Id: BrachyDetectorConstructionLeipzig.cc,v 1.2 2004/05/25 08:36:18 guatelli Exp $
// GEANT4 tag $Name: geant4-07-01 $
//

#include "globals.hh"
#include "BrachyDetectorConstructionLeipzig.hh"
#include "G4CSGSolid.hh"
#include "G4Sphere.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4SDManager.hh"
#include "G4SubtractionSolid.hh"
#include "G4RunManager.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4MaterialTable.hh"
#include "Randomize.hh"  
#include "G4RunManager.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4PVParameterised.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4UnionSolid.hh"
#include "BrachyMaterial.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

// Leipzig Applicator ...

BrachyDetectorConstructionLeipzig::BrachyDetectorConstructionLeipzig()
: capsule(0), capsuleLog(0), capsulePhys(0),
  capsuleTip(0), capsuleTipLog(0), capsuleTipPhys(0),
  iridiumCore(0), iridiumCoreLog(0), iridiumCorePhys(0),
  applicator1(0), applicator1Log(0), applicator1Phys(0),
  applicator2(0), applicator2Log(0), applicator2Phys(0)
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

  G4Material* capsuleMat = pMaterial -> GetMat("Stainless steel");
  G4Material* iridium = pMaterial -> GetMat("Iridium");
  G4Material* tungsten =pMaterial -> GetMat("Tungsten");

  //Iridium source ...

  capsule = new G4Tubs("Capsule",0,0.55*mm,3.725*mm,0.*deg,360.*deg);
  capsuleLog = new G4LogicalVolume(capsule,capsuleMat,"CapsuleLog");
  capsulePhys = new G4PVPlacement(0,
                                  G4ThreeVector(0,0,-1.975*mm),
                                  "CapsulePhys",
				  capsuleLog,
                                  mother, //mother volume: phantom
                                  false,
                                  0);

  // Capsule tip
  capsuleTip = new G4Sphere("CapsuleTip",0.*mm,0.55*mm,0.*deg,360.*deg,0.*deg,90.*deg);
  capsuleTipLog = new G4LogicalVolume(capsuleTip,capsuleMat,"CapsuleTipLog");
  capsuleTipPhys = new G4PVPlacement(0,
                                     G4ThreeVector(0.,0.,1.75*mm),
                                     "CapsuleTipPhys",
				     capsuleTipLog,
                                     mother,
                                     false,
                                     0);
  // Iridium core

  iridiumCore = new G4Tubs("IrCore",0,0.30*mm,1.75*mm,0.*deg,360.*deg);
  iridiumCoreLog = new G4LogicalVolume(iridiumCore,
                                       iridium,
                                       "IridiumCoreLog");
  iridiumCorePhys = new G4PVPlacement(0,
                                      G4ThreeVector(0.,0.,1.975*mm),
                                      "IridiumCorePhys",
				      iridiumCoreLog,
                                      capsulePhys,
                                      false,
                                      0);
  // Applicator
  //Leipzig Applicator is given by two different volumes

  applicator1 = new G4Tubs("Appl1",5*mm,10.5*mm,12*mm,0.*deg,360.*deg);
  applicator1Log = new G4LogicalVolume(applicator1,tungsten,"Appl1Log");
  applicator1Phys = new G4PVPlacement(0,
                                      G4ThreeVector(0,0,4.0*mm),
                                      "Appl1Phys", 
                                      applicator1Log,
				      mother,
                                      false,
                                      0);

  applicator2 = new G4Tubs("Appl2",0.55*mm,5.*mm,3.125*mm,0.*deg,360.*deg);
  applicator2Log = new G4LogicalVolume(applicator2,tungsten,"Appl2");
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
}
