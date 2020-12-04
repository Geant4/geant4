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
//
// Code by S. Guatelli
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

BrachyDetectorConstructionLeipzig::BrachyDetectorConstructionLeipzig(): 
fCapsule(nullptr), fCapsuleTip(nullptr), fIridiumCore(nullptr), fApplicator1(nullptr), fApplicator2(nullptr),
fCapsuleLog(nullptr), fCapsuleTipLog(nullptr), fIridiumCoreLog(nullptr), fApplicator1Log(nullptr),
fApplicator2Log(nullptr), fCapsulePhys(nullptr), fCapsuleTipPhys(nullptr), fIridiumCorePhys(nullptr), 
fApplicator1Phys(nullptr), fApplicator2Phys(nullptr)
{ 
  fMaterial = new BrachyMaterial();
}

BrachyDetectorConstructionLeipzig::~BrachyDetectorConstructionLeipzig()
{ 
  delete fMaterial; 
}

void  BrachyDetectorConstructionLeipzig::ConstructLeipzig(G4VPhysicalVolume* mother)
{
  G4Colour  red     (1.0, 0.0, 0.0) ; 
  G4Colour  lblue   (0.0, 0.0, .75);

  G4Material* capsuleMat = fMaterial -> GetMat("Stainless steel");
  G4Material* iridium = fMaterial -> GetMat("Iridium");
  G4Material* tungsten = fMaterial -> GetMat("Tungsten");

  //Iridium source ...

  fCapsule = new G4Tubs("Capsule",0,0.55*mm,3.725*mm,0.*deg,360.*deg);
  fCapsuleLog = new G4LogicalVolume(fCapsule,capsuleMat,"CapsuleLog");
  fCapsulePhys = new G4PVPlacement(nullptr, G4ThreeVector(0,0,-1.975*mm),"CapsulePhys",
                                  fCapsuleLog,mother, //mother volume: phantom
                                  false,0, true);

  // Capsule tip
  fCapsuleTip = new G4Sphere("CapsuleTip",0.*mm,0.55*mm,0.*deg,360.*deg,0.*deg,90.*deg);
  fCapsuleTipLog = new G4LogicalVolume(fCapsuleTip,capsuleMat,"CapsuleTipLog");
  fCapsuleTipPhys = new G4PVPlacement(nullptr,G4ThreeVector(0.,0.,1.75*mm),"CapsuleTipPhys",
                                      fCapsuleTipLog,mother,false,0, true);
  // Iridium core
  fIridiumCore = new G4Tubs("IrCore",0,0.30*mm,1.75*mm,0.*deg,360.*deg);
  fIridiumCoreLog = new G4LogicalVolume(fIridiumCore, iridium, "IridiumCoreLog");
  fIridiumCorePhys = new G4PVPlacement(nullptr,G4ThreeVector(0.,0.,1.975*mm),"IridiumCorePhys",
                                      fIridiumCoreLog, fCapsulePhys,false,0, true);

  //Leipzig Applicator is modelled with two different volumes
  fApplicator1 = new G4Tubs("Appl1",5*mm,10.5*mm,12*mm,0.*deg,360.*deg);
  fApplicator1Log = new G4LogicalVolume(fApplicator1,tungsten,"Appl1Log");
  fApplicator1Phys = new G4PVPlacement(nullptr,G4ThreeVector(0,0,4.0*mm),"Appl1Phys",fApplicator1Log,
                                      mother,false,0, true);

  fApplicator2 = new G4Tubs("Appl2",0.55*mm,5.*mm,3.125*mm,0.*deg,360.*deg);
  fApplicator2Log = new G4LogicalVolume(fApplicator2,tungsten,"Appl2");
  fApplicator2Phys = new G4PVPlacement(nullptr,G4ThreeVector(0,0,-4.875*mm),
                                      "Appl2Phys", fApplicator2Log,mother,false,0, true);

  fSimpleCapsuleVisAtt = new G4VisAttributes(red); 
  fSimpleCapsuleVisAtt -> SetVisibility(true); 
  fSimpleCapsuleVisAtt -> SetForceWireframe(true); 
  fCapsuleLog -> SetVisAttributes(fSimpleCapsuleVisAtt); 

  fSimpleCapsuleTipVisAtt = new G4VisAttributes(red); 
  fSimpleCapsuleTipVisAtt -> SetVisibility(true); 
  fSimpleCapsuleTipVisAtt -> SetForceSolid(true); 
  fCapsuleTipLog -> SetVisAttributes(fSimpleCapsuleTipVisAtt);
  fIridiumCoreLog -> SetVisAttributes(fSimpleCapsuleTipVisAtt);

  fApplicatorVisAtt = new G4VisAttributes(lblue);
  fApplicatorVisAtt -> SetVisibility(true);
  fApplicatorVisAtt -> SetForceWireframe(true);
  fApplicator1Log -> SetVisAttributes(fApplicatorVisAtt);
  fApplicator2Log -> SetVisAttributes(fApplicatorVisAtt);
}
void BrachyDetectorConstructionLeipzig::CleanLeipzigApplicator()
{
delete fApplicatorVisAtt; fApplicatorVisAtt = nullptr;
delete fSimpleCapsuleTipVisAtt; fSimpleCapsuleTipVisAtt = nullptr;
delete fSimpleCapsuleVisAtt; fSimpleCapsuleVisAtt = nullptr;
delete fApplicator2Phys; fApplicator2Phys = nullptr;
delete fApplicator1Phys; fApplicator1Phys = nullptr;
delete fIridiumCorePhys; fIridiumCorePhys = nullptr;
delete fCapsuleTipPhys; fCapsuleTipPhys = nullptr;
delete fCapsulePhys; fCapsulePhys = nullptr;
delete fApplicator2Log; fApplicator2Log = nullptr;
delete fApplicator1Log; fApplicator1Log = nullptr;
delete fIridiumCoreLog; fIridiumCoreLog = nullptr;
delete fCapsuleTipLog; fCapsuleTipLog = nullptr;
delete fCapsuleLog; fCapsuleLog = nullptr;
}
