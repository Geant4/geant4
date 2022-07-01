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
#include "BrachyDetectorConstructionI.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
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
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4NistManager.hh"

BrachyDetectorConstructionI::BrachyDetectorConstructionI():
fDefaultTub(nullptr), fCapsule(nullptr), fCapsuleTip(nullptr), fIodineCore(nullptr), fDefaultTubLog(nullptr),
fCapsuleLog(nullptr), fCapsuleTipLog(nullptr), fIodineCoreLog(nullptr), fDefaultTubPhys(nullptr), 
fCapsulePhys(nullptr),fCapsuleTipPhys1(nullptr),fCapsuleTipPhys2(nullptr), fIodineCorePhys(nullptr),
fSimpleIodineVisAtt(nullptr), fSimpleCapsuleVisAtt(nullptr), fSimpleCapsuleTipVisAtt(nullptr)
{}

void BrachyDetectorConstructionI::ConstructIodine(G4VPhysicalVolume* mother)
{
  // Model of the Bebig Isoseed I-125 brachy source
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* titanium = nist -> FindOrBuildMaterial("G4_Ti");
  G4Material* iodine = nist -> FindOrBuildMaterial("G4_I");
  G4Material* air = nist -> FindOrBuildMaterial("G4_AIR");
  
  G4Colour  red     (1.0, 0.0, 0.0) ;
  G4Colour  magenta (1.0, 0.0, 1.0) ; 
  G4Colour  lblue   (0.0, 0.0, .75);

  // Air tub
  fDefaultTub = new G4Tubs("DefaultTub",0.*mm, 0.40*mm, 1.84*mm, 0.*deg, 360.*deg);
  fDefaultTubLog = new G4LogicalVolume(fDefaultTub,air,"DefaultTub_Log");
  fDefaultTubPhys = new G4PVPlacement(nullptr,
                                      G4ThreeVector(),
                                      "defaultTub_Phys",
                                      fDefaultTubLog,
                                      mother,
                                      false,
                                      0, true); 
  // Capsule main body ...
  fCapsule = new G4Tubs("fCapsule", 0.35*mm,0.40*mm,1.84*mm,0.*deg,360.*deg);
  fCapsuleLog = new G4LogicalVolume(fCapsule,titanium,"fCapsuleLog");
  fCapsulePhys = new G4PVPlacement(nullptr,
                                  G4ThreeVector(),
                                  "fCapsulePhys",
                                  fCapsuleLog,
                                  fDefaultTubPhys,
                                  false,
                                  0, true);
  // fCapsule tips
  fCapsuleTip = new G4Sphere("fCapsuleTip",
                            0.*mm,
                            0.40*mm,
                            0.*deg,
                            360.*deg,
                            0.*deg,
                            90.*deg);
                            
  fCapsuleTipLog = new G4LogicalVolume(fCapsuleTip,titanium,"fCapsuleTipLog");
  fCapsuleTipPhys1 = new G4PVPlacement(nullptr,
                                      G4ThreeVector(0.,0.,1.84*mm),
                                      "IodinefCapsuleTipPhys1",
                                      fCapsuleTipLog,
                                      mother,
                                      false,
                                      0, true);

  auto rotateMatrix = new G4RotationMatrix();
  rotateMatrix -> rotateX(180.0*deg);
  fCapsuleTipPhys2 = new G4PVPlacement(rotateMatrix, 
                                      G4ThreeVector(0,0,-1.84*mm),
                                      "IodinefCapsuleTipPhys2",
                                      fCapsuleTipLog,
                                      mother,
                                      false,
                                      0, true);
 
  // Radiactive core ...
  fIodineCore = new G4Tubs("ICore",0.085*mm,0.35*mm,1.75*mm,0.*deg,360.*deg);
  fIodineCoreLog = new G4LogicalVolume(fIodineCore,iodine,"IodineCoreLog");
  fIodineCorePhys = new G4PVPlacement(nullptr,
                                     G4ThreeVector(0.,0.,0.),
                                     "IodineCorePhys",
                                     fIodineCoreLog,
                                     fDefaultTubPhys,
                                     false,
                                     0, true);
 
  // Visual attributes ...
  
  fSimpleIodineVisAtt = new G4VisAttributes(magenta);
  fSimpleIodineVisAtt -> SetVisibility(true);
  fSimpleIodineVisAtt -> SetForceSolid(true);
  fIodineCoreLog -> SetVisAttributes(fSimpleIodineVisAtt);

  fSimpleCapsuleVisAtt = new G4VisAttributes(red);
  fSimpleCapsuleVisAtt -> SetVisibility(true);  
  fSimpleCapsuleVisAtt -> SetForceWireframe(true);
  fCapsuleLog -> SetVisAttributes( fSimpleCapsuleVisAtt);

  fSimpleCapsuleTipVisAtt = new G4VisAttributes(red);
  fSimpleCapsuleTipVisAtt -> SetVisibility(true); 
  fSimpleCapsuleTipVisAtt -> SetForceSolid(true);
  fCapsuleTipLog -> SetVisAttributes(fSimpleCapsuleTipVisAtt);
}

void BrachyDetectorConstructionI::CleanIodine()
{
delete fSimpleIodineVisAtt; fSimpleIodineVisAtt = nullptr;
delete fSimpleCapsuleVisAtt; fSimpleCapsuleVisAtt = nullptr;
delete fSimpleCapsuleTipVisAtt; fSimpleCapsuleTipVisAtt = nullptr;
delete fCapsuleTipPhys1; fCapsuleTipPhys1 = nullptr; 
delete fCapsuleTipPhys2; fCapsuleTipPhys2 = nullptr;
delete fIodineCorePhys; fIodineCorePhys = nullptr;
delete fCapsulePhys; fCapsulePhys = nullptr;
delete fDefaultTubPhys; fDefaultTubPhys = nullptr; 
delete fDefaultTubLog; fDefaultTubLog = nullptr;
delete fCapsuleLog; fCapsuleLog = nullptr;
delete fCapsuleTipLog; fCapsuleTipLog = nullptr;
delete fIodineCoreLog; fIodineCoreLog = nullptr;
delete fDefaultTub; fDefaultTub = nullptr;
delete fCapsule; fCapsule = nullptr;
delete fCapsuleTip; fCapsuleTip = nullptr;
delete fIodineCore; fIodineCore = nullptr;

G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}
