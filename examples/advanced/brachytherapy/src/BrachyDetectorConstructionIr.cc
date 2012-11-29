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
// S. Agostinelli, F. Foppiano, S. Garelli , M. Tropeano, S.Guatelli
//
//    ****************************************
//    *                                      *
//    *    BrachyDetectorConstructionIr.cc   *
//    *                                      *
//    ****************************************
//
// $Id$
//
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "BrachyDetectorConstructionIr.hh"
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

BrachyDetectorConstructionIr::BrachyDetectorConstructionIr()
  : 
   capsule(0),capsuleLog(0),
   capsulePhys(0),
   capsuleTip(0),capsuleTipLog(0),
   capsuleTipPhys(0),
   iridiumCore(0),iridiumCoreLog(0),
   iridiumCorePhys(0),
   simpleCapsuleVisAtt(0),simpleCapsuleTipVisAtt(0),simpleIridiumVisAtt(0)
{
  pMat = new BrachyMaterial();
}

BrachyDetectorConstructionIr::~BrachyDetectorConstructionIr()
{ 
  delete pMat; 
}

void BrachyDetectorConstructionIr::ConstructIridium(G4VPhysicalVolume* mother)
{
  G4Colour  red     (1.0, 0.0, 0.0) ;
  G4Colour  magenta (1.0, 0.0, 1.0) ; 

  G4Material* capsuleMat = pMat -> GetMat("Stainless steel");
  G4Material* iridiumMat = pMat -> GetMat("Iridium");

  // Capsule main body
  capsule = new G4Tubs("Capsule",0,0.55*mm,3.725*mm,0.*deg,360.*deg);
  capsuleLog = new G4LogicalVolume(capsule,capsuleMat,"CapsuleLog");
  capsulePhys = new G4PVPlacement(0,
                                 G4ThreeVector(0,0,-1.975*mm),
                                 "CapsulePhys",
                                 capsuleLog,
                                 mother,
                                 false,
                                 0);

  // Capsule tip
  capsuleTip = new G4Sphere("CapsuleTipIridium",
                            0.*mm,
                            0.55*mm,
                            0.*deg,
                            360.*deg,
                            0.*deg,
                            90.*deg); 
  
  capsuleTipLog = new G4LogicalVolume(capsuleTip,
                                      capsuleMat,
                                      "CapsuleTipIridumLog");
  capsuleTipPhys = new G4PVPlacement(0,
                                     G4ThreeVector(0.,0.,1.75*mm),
                                     "CapsuleTipIridiumPhys",
                                     capsuleTipLog,
                                     mother,
                                     false,
                                     0);

  // Iridium core
  iridiumCore = new G4Tubs("IrCore",0,0.30*mm,1.75*mm,0.*deg,360.*deg);
  iridiumCoreLog = new G4LogicalVolume(iridiumCore,
                                       iridiumMat,
                                       "IridiumCoreLog");
  iridiumCorePhys = new G4PVPlacement(0,
                                      G4ThreeVector(),
                                      "IridiumCorePhys",
                                      iridiumCoreLog,
                                      capsulePhys,
                                      false,
                                      0);

  simpleCapsuleVisAtt = new G4VisAttributes(red);
  simpleCapsuleVisAtt -> SetVisibility(true);  
  simpleCapsuleVisAtt -> SetForceWireframe(true);
  capsuleLog -> SetVisAttributes(simpleCapsuleVisAtt);

  simpleCapsuleTipVisAtt = new G4VisAttributes(red);
  simpleCapsuleTipVisAtt -> SetVisibility(true);  
  simpleCapsuleTipVisAtt -> SetForceSolid(true);
  capsuleTipLog -> SetVisAttributes(simpleCapsuleTipVisAtt);

  simpleIridiumVisAtt = new G4VisAttributes(magenta);
  simpleIridiumVisAtt -> SetVisibility(true);
  simpleIridiumVisAtt -> SetForceWireframe(true);
  iridiumCoreLog -> SetVisAttributes(simpleIridiumVisAtt);
}

void BrachyDetectorConstructionIr::CleanIridium()
{ 
  
  delete simpleIridiumVisAtt; 
  simpleIridiumVisAtt = 0;
  
  delete iridiumCorePhys; 
  iridiumCorePhys = 0;
  
  delete iridiumCore; 
  iridiumCore = 0;
  
  delete iridiumCoreLog; 
  iridiumCoreLog = 0 ;

  delete simpleCapsuleTipVisAtt; 
  simpleCapsuleTipVisAtt = 0;
  
  delete capsuleTipPhys; 
  capsuleTipPhys = 0;
  
  delete capsuleTip;
  capsuleTip = 0;
  delete capsuleTipLog; 
  capsuleTipLog = 0;

  delete simpleCapsuleVisAtt; 
  simpleCapsuleVisAtt = 0;
  
  delete capsulePhys; 
  capsulePhys = 0;
  
  delete capsule; 
  capsule = 0;
  
  delete capsuleLog; 
  capsuleLog = 0;
}
