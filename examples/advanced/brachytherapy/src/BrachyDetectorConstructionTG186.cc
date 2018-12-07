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
// D. Cutajar
//
//    ****************************************
//    *                                      *
//    *    BrachyDetectorConstructionTG186.cc   *
//    *                                      *
//    ****************************************
//
//
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "BrachyDetectorConstructionTG186.hh"
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

BrachyDetectorConstructionTG186::BrachyDetectorConstructionTG186()
  : 
   TG186capsule(0),TG186capsuleLog(0),
   TG186capsulePhys(0),
   TG186capsuleTip(0),TG186capsuleTipLog(0),
   TG186capsuleTipPhys(0),
   TG186iridiumCore(0),TG186iridiumCoreLog(0),
   TG186iridiumCorePhys(0),
   TG186cable(0),TG186cableLog(0),
   TG186cablePhys(0),
   TG186simpleCapsuleVisAtt(0),TG186simpleCapsuleTipVisAtt(0),TG186simpleIridiumVisAtt(0),
   TG186simpleCableVisAtt(0)
{
  pMat = new BrachyMaterial();
}

BrachyDetectorConstructionTG186::~BrachyDetectorConstructionTG186()
{ 
  delete pMat; 
}

void BrachyDetectorConstructionTG186::ConstructTG186(G4VPhysicalVolume* mother)
{
  G4Colour  red     (1.0, 0.0, 0.0) ;
  G4Colour  magenta (1.0, 0.0, 1.0) ; 

  G4Material* capsuleMat = pMat -> GetMat("Stainless steel");
  G4Material* iridiumMat = pMat -> GetMat("Iridium");

  // Capsule main body
  TG186capsule = new G4Tubs("TG186-Capsule",0,0.5*mm,2.25*mm,0.*deg,360.*deg);
  TG186capsuleLog = new G4LogicalVolume(TG186capsule,capsuleMat,"TG186-CapsuleLog");
  TG186capsulePhys = new G4PVPlacement(0,
                                 G4ThreeVector(0,0,-0.4*mm),
                                 "TG186-IridiumCapsulePhys",
                                 TG186capsuleLog,
                                 mother,
                                 false,
                                 0, true);

  // Capsule tip
  TG186capsuleTip = new G4Sphere("Tg186-CapsuleTipIridium",
                            0.*mm,
                            0.5*mm,
                            0.*deg,
                            360.*deg,
                            0.*deg,
                            90.*deg); 
  
  TG186capsuleTipLog = new G4LogicalVolume(TG186capsuleTip,
                                      capsuleMat,
                                      "CapsuleTipIridumLog");
  TG186capsuleTipPhys = new G4PVPlacement(0,
                                     G4ThreeVector(0.,0.,1.85*mm),
                                     "TG186-CapsuleTipIridiumPhys",
                                     TG186capsuleTipLog,
                                     mother,
                                     false,
                                     0, true);
									 
  TG186cable = new G4Tubs("TG186-cable",
                            0.*mm,
                            0.5*mm,
			    1.0*mm,
                            0.*deg,
                            360.*deg); 
  
  TG186cableLog = new G4LogicalVolume(TG186cable,
                                      capsuleMat,
                                      "TG186-cableLog");
  TG186cablePhys = new G4PVPlacement(0,
                                     G4ThreeVector(0.,0.,-3.65*mm),
                                     "TG186-CablePhys",
                                     TG186cableLog,
                                     mother,
                                     false,
                                     0, true);								 

  // Iridium core
  TG186iridiumCore = new G4Tubs("TG186-IrCore",0,0.30*mm,1.75*mm,0.*deg,360.*deg);
  TG186iridiumCoreLog = new G4LogicalVolume(TG186iridiumCore,
                                       iridiumMat,
                                       "TG186-IridiumCoreLog");
  TG186iridiumCorePhys = new G4PVPlacement(0,
                                      G4ThreeVector(0,0,0.4*mm),
                                      "TG186-IridiumCorePhys",
                                      TG186iridiumCoreLog,
                                      TG186capsulePhys,
                                      false,
                                      0, true);

  TG186simpleCapsuleVisAtt = new G4VisAttributes(red);
  TG186simpleCapsuleVisAtt -> SetVisibility(true);  
  TG186simpleCapsuleVisAtt -> SetForceWireframe(true);
  TG186capsuleLog -> SetVisAttributes(TG186simpleCapsuleVisAtt);

  TG186simpleCapsuleTipVisAtt = new G4VisAttributes(red);
  TG186simpleCapsuleTipVisAtt -> SetVisibility(true);  
  TG186simpleCapsuleTipVisAtt -> SetForceSolid(true);
  TG186capsuleTipLog -> SetVisAttributes(TG186simpleCapsuleTipVisAtt);

  TG186simpleIridiumVisAtt = new G4VisAttributes(magenta);
  TG186simpleIridiumVisAtt -> SetVisibility(true);
  TG186simpleIridiumVisAtt -> SetForceWireframe(true);
  TG186iridiumCoreLog -> SetVisAttributes(TG186simpleIridiumVisAtt);
  
  TG186simpleCableVisAtt = new G4VisAttributes(red);
  TG186simpleCableVisAtt -> SetVisibility(true);  
  TG186simpleCableVisAtt -> SetForceSolid(true);
  TG186cableLog -> SetVisAttributes(TG186simpleCableVisAtt);
}

void BrachyDetectorConstructionTG186::CleanTG186()
{ 
  
  delete TG186simpleIridiumVisAtt; 
  TG186simpleIridiumVisAtt = 0;
  
  delete TG186iridiumCorePhys; 
  TG186iridiumCorePhys = 0;
  
  delete TG186iridiumCore; 
  TG186iridiumCore = 0;
  
  delete TG186iridiumCoreLog; 
  TG186iridiumCoreLog = 0 ;

  delete TG186simpleCapsuleTipVisAtt; 
  TG186simpleCapsuleTipVisAtt = 0;
  
  delete TG186capsuleTipPhys; 
  TG186capsuleTipPhys = 0;
  
  delete TG186capsuleTip;
  TG186capsuleTip = 0;
  
  delete TG186capsuleTipLog; 
  TG186capsuleTipLog = 0;

  delete TG186simpleCapsuleVisAtt; 
  TG186simpleCapsuleVisAtt = 0;
  
  delete TG186capsulePhys; 
  TG186capsulePhys = 0;
  
  delete TG186capsule; 
  TG186capsule = 0;
  
  delete TG186capsuleLog; 
  TG186capsuleLog = 0;
  
  delete TG186cable;
  TG186cable = 0;
  
  delete TG186cableLog;
  TG186cableLog = 0; 
   
  delete TG186cablePhys;
  TG186cablePhys = 0;
  
  delete TG186simpleCableVisAtt; 
  TG186simpleCableVisAtt = 0;

  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}
