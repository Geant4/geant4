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
// Code developed by:
// S. Agostinelli, F. Foppiano, S. Garelli , M. Tropeano, S.Guatelli
//
//    ****************************************
//    *                                      *
//    *    BrachyDetectorConstructionIr.cc   *
//    *                                      *
//    ****************************************
//
// $Id: BrachyDetectorConstructionIr.cc,v 1.8 2004/03/11 15:38:42 guatelli Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
#include "globals.hh"
#include "BrachyDetectorConstructionIr.hh"
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

BrachyDetectorConstructionIr::BrachyDetectorConstructionIr()
 : capsule(0),capsuleLog(0),capsulePhys(0),
   capsuleTip(0),capsuleTipLog(0),capsuleTipPhys(0),
   iridiumCore(0),iridiumCoreLog(0),iridiumCorePhys(0),
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

  G4Material* capsuleMat = pMat->GetMat("Stainless steel");
  G4Material* iridiumMat = pMat->GetMat("Iridium");

  // Capsule main body
  capsule = new G4Tubs("Capsule",0,0.55*mm,3.725*mm,0.*deg,360.*deg);
  capsuleLog = new G4LogicalVolume(capsule,capsuleMat,"CapsuleLog");
  capsulePhys= new G4PVPlacement(0,
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
  simpleCapsuleVisAtt->SetVisibility(true);  
  simpleCapsuleVisAtt->SetForceWireframe(true);
  capsuleLog->SetVisAttributes(simpleCapsuleVisAtt);

  simpleCapsuleTipVisAtt = new G4VisAttributes(red);
  simpleCapsuleTipVisAtt->SetVisibility(true);  
  simpleCapsuleTipVisAtt->SetForceSolid(true);
  capsuleTipLog->SetVisAttributes(simpleCapsuleTipVisAtt);

  simpleIridiumVisAtt = new G4VisAttributes(magenta);
  simpleIridiumVisAtt->SetVisibility(true);
  simpleIridiumVisAtt->SetForceWireframe(true);
  iridiumCoreLog->SetVisAttributes(simpleIridiumVisAtt);
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
