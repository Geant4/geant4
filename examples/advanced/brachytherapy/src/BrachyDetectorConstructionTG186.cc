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
#include "BrachyDetectorConstructionTG186.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
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
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4NistManager.hh"

BrachyDetectorConstructionTG186::BrachyDetectorConstructionTG186(): 
   fTG186capsule(nullptr), fTG186capsuleLog(nullptr),
   fTG186capsulePhys(nullptr),
   fTG186capsuleTip(nullptr), fTG186capsuleTipLog(nullptr),
   fTG186capsuleTipPhys(nullptr),
   fTG186iridiumCore(nullptr), fTG186iridiumCoreLog(nullptr),
   fTG186iridiumCorePhys(nullptr),
   fTG186cable(nullptr), fTG186cableLog(nullptr),
   fTG186cablePhys(nullptr),
   fTG186simpleCapsuleVisAtt(nullptr), fTG186simpleCapsuleTipVisAtt(nullptr), fTG186simpleIridiumVisAtt(nullptr),
   fTG186simpleCableVisAtt(nullptr)
{}

void BrachyDetectorConstructionTG186::ConstructTG186(G4VPhysicalVolume* mother)
{
  G4Colour  red     (1.0, 0.0, 0.0) ;
  G4Colour  magenta (1.0, 0.0, 1.0) ; 

  G4NistManager* nist = G4NistManager::Instance();
  auto iridium = nist -> FindOrBuildMaterial("G4_Ir");
 
  // Stainless steel (Medical Physics, Vol 25, No 10, Oct 1998)
  constexpr G4double d = 8.02*g/cm3;
  G4int Z; //atomic number of the element
  G4Element* elMn = nist -> FindOrBuildElement(Z=12);
  G4Element* elSi = nist -> FindOrBuildElement(Z=14);
  G4Element* elCr = nist -> FindOrBuildElement(Z=24);
  G4Element* elFe = nist -> FindOrBuildElement(Z=26);
  G4Element* elNi = nist -> FindOrBuildElement(Z=28);
  auto capsuleMat = new G4Material("Stainless steel",d,5);
  capsuleMat -> AddElement(elMn, 0.02);
  capsuleMat -> AddElement(elSi, 0.01);
  capsuleMat -> AddElement(elCr, 0.19);
  capsuleMat -> AddElement(elNi, 0.10);
  capsuleMat -> AddElement(elFe, 0.68);
  

  // Capsule main body
  fTG186capsule = new G4Tubs("TG186-Capsule",0,0.5*mm,2.25*mm,0.*deg,360.*deg);
  fTG186capsuleLog = new G4LogicalVolume(fTG186capsule,capsuleMat,"TG186-CapsuleLog");
  fTG186capsulePhys = new G4PVPlacement(nullptr,
                                 G4ThreeVector(0,0,-0.4*mm),
                                 "TG186-IridiumCapsulePhys",
                                 fTG186capsuleLog,
                                 mother,
                                 false,
                                 0, true);

  // Capsule tip
  fTG186capsuleTip = new G4Sphere("Tg186-CapsuleTipIridium",
                            0.*mm,
                            0.5*mm,
                            0.*deg,
                            360.*deg,
                            0.*deg,
                            90.*deg); 
  
  fTG186capsuleTipLog = new G4LogicalVolume(fTG186capsuleTip,
                                      capsuleMat,
                                      "CapsuleTipIridumLog");
                                      
  fTG186capsuleTipPhys = new G4PVPlacement(nullptr,
                                     G4ThreeVector(0.,0.,1.85*mm),
                                     "TG186-CapsuleTipIridiumPhys",
                                     fTG186capsuleTipLog,
                                     mother,
                                     false,
                                     0, true);
									 
  fTG186cable = new G4Tubs("TG186-cable",
                            0.*mm,
                            0.5*mm,
	                   		    1.0*mm,
                            0.*deg,
                            360.*deg); 
  
  fTG186cableLog = new G4LogicalVolume(fTG186cable,
                                      capsuleMat,
                                      "TG186-cableLog");
                                      
  fTG186cablePhys = new G4PVPlacement(nullptr,
                                     G4ThreeVector(0.,0.,-3.65*mm),
                                     "TG186-CablePhys",
                                     fTG186cableLog,
                                     mother,
                                     false,
                                     0, true);								 

  // Iridium core
  fTG186iridiumCore = new G4Tubs("TG186-IrCore",0,0.30*mm,1.75*mm,0.*deg,360.*deg);
  
  fTG186iridiumCoreLog = new G4LogicalVolume(fTG186iridiumCore,
                                       iridium,
                                       "TG186-IridiumCoreLog");
                                       
  fTG186iridiumCorePhys = new G4PVPlacement(nullptr,
                                      G4ThreeVector(0,0,0.4*mm),
                                      "TG186-IridiumCorePhys",
                                      fTG186iridiumCoreLog,
                                      fTG186capsulePhys,
                                      false,
                                      0, true);

  fTG186simpleCapsuleVisAtt = new G4VisAttributes(red);
  fTG186simpleCapsuleVisAtt -> SetVisibility(true);  
  fTG186simpleCapsuleVisAtt -> SetForceWireframe(true);
  fTG186capsuleLog -> SetVisAttributes(fTG186simpleCapsuleVisAtt);

  fTG186simpleCapsuleTipVisAtt = new G4VisAttributes(red);
  fTG186simpleCapsuleTipVisAtt -> SetVisibility(true);  
  fTG186simpleCapsuleTipVisAtt -> SetForceSolid(true);
  fTG186capsuleTipLog -> SetVisAttributes(fTG186simpleCapsuleTipVisAtt);

  fTG186simpleIridiumVisAtt = new G4VisAttributes(magenta);
  fTG186simpleIridiumVisAtt -> SetVisibility(true);
  fTG186simpleIridiumVisAtt -> SetForceWireframe(true);
  fTG186iridiumCoreLog -> SetVisAttributes(fTG186simpleIridiumVisAtt);
  
  fTG186simpleCableVisAtt = new G4VisAttributes(red);
  fTG186simpleCableVisAtt -> SetVisibility(true);  
  fTG186simpleCableVisAtt -> SetForceSolid(true);
  fTG186cableLog -> SetVisAttributes(fTG186simpleCableVisAtt);
}

void BrachyDetectorConstructionTG186::CleanTG186()
{ 
  
  delete fTG186simpleIridiumVisAtt; 
  fTG186simpleIridiumVisAtt = nullptr;
  
  delete fTG186iridiumCorePhys; 
  fTG186iridiumCorePhys = nullptr;
  
  delete fTG186iridiumCore; 
  fTG186iridiumCore = nullptr;
  
  delete fTG186iridiumCoreLog; 
  fTG186iridiumCoreLog = nullptr ;

  delete fTG186simpleCapsuleTipVisAtt; 
  fTG186simpleCapsuleTipVisAtt = nullptr;
  
  delete fTG186capsuleTipPhys; 
  fTG186capsuleTipPhys = nullptr;
  
  delete fTG186capsuleTip;
  fTG186capsuleTip = nullptr;
  
  delete fTG186capsuleTipLog; 
  fTG186capsuleTipLog = nullptr;

  delete fTG186simpleCapsuleVisAtt; 
  fTG186simpleCapsuleVisAtt = nullptr;
  
  delete fTG186capsulePhys; 
  fTG186capsulePhys = nullptr;
  
  delete fTG186capsule; 
  fTG186capsule = nullptr;
  
  delete fTG186capsuleLog; 
  fTG186capsuleLog = nullptr;
  
  delete fTG186cable;
  fTG186cable = nullptr;
  
  delete fTG186cableLog;
  fTG186cableLog = nullptr; 
   
  delete fTG186cablePhys;
  fTG186cablePhys = nullptr;
  
  delete fTG186simpleCableVisAtt; 
  fTG186simpleCableVisAtt = nullptr;

  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}
