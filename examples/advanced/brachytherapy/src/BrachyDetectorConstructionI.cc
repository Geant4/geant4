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
// $Id: BrachyDetectorConstructionI.cc,v 1.4 2003-05-09 17:36:41 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "globals.hh"
#include "BrachyDetectorConstructionI.hh"
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

BrachyDetectorConstructionI::BrachyDetectorConstructionI()
: DefaultTub(0), DefaultTub_log(0),DefaultTub_Phys(0),
  Capsule(0),CapsuleLog(0),CapsulePhys(0),
  CapsuleTip(0),CapsuleTipLog(0),CapsuleTipPhys1(0),CapsuleTipPhys2(0),
  IodiumCore(0),ICoreLog(0),ICorePhys(0),
  Marker(0),MarkerLog(0),MarkerPhys(0)
{ 
  pMat= new BrachyMaterial();
}

//....
BrachyDetectorConstructionI::~BrachyDetectorConstructionI()
{ 
  delete pMat; 
}

//....
void BrachyDetectorConstructionI::ConstructIodium(G4VPhysicalVolume* mother)
{
  // Volumes

  G4cout << "Creation of Iodium ..." << G4endl;
  G4cout << "# of daughters of world - "
         << mother->GetMother()->GetLogicalVolume()->GetNoDaughters() << G4endl;
  G4cout << "# of daughters of phantom - "
         << mother->GetLogicalVolume()->GetNoDaughters() << G4endl;

  G4Material* titanium=pMat->GetMat("titanium");
  G4Material* air =pMat->GetMat("Air");
  G4Material* Iodio=pMat->GetMat("Iodium");
  G4Material* Gold=pMat->GetMat("gold");

  G4Colour  red     (1.0, 0.0, 0.0) ;
  G4Colour  magenta (1.0, 0.0, 1.0) ; 
  G4Colour  lblue   (0.0, 0.0, .75);

  DefaultTub = new G4Tubs("DefaultTub",0.*mm,0.40*mm,1.84*mm,0.*deg,360.*deg);
  DefaultTub_log = new G4LogicalVolume(DefaultTub,air,"DefaultTub_Log");
  DefaultTub_Phys = new G4PVPlacement(0,G4ThreeVector(),"DefaultTub_Phys",
                                      DefaultTub_log,mother,false,0); 

  G4cout << "Created DefaultTub_Phys" << G4endl;
  G4cout << "# of daughters - " << mother->GetLogicalVolume()->GetNoDaughters() << G4endl;

  //  Capsule main body

  G4double CapsuleR=0.35*mm;
  Capsule = new G4Tubs("Capsule", CapsuleR,0.40*mm,1.84*mm,0.*deg,360.*deg);
  CapsuleLog = new G4LogicalVolume(Capsule,titanium,"CapsuleLog");
  CapsulePhys = new G4PVPlacement(0,G4ThreeVector(),"CapsulePhys",
                                  CapsuleLog,DefaultTub_Phys,false,0);

  G4cout << "Created CapsulePhys" << G4endl;
  G4cout << "# of daughters - " << mother->GetLogicalVolume()->GetNoDaughters() << G4endl;

  // Capsule tips

  CapsuleTip = new G4Sphere("CapsuleTip",0.*mm,0.40*mm,0.*deg,360.*deg,0.*deg,90.*deg);
  CapsuleTipLog = new G4LogicalVolume(CapsuleTip,titanium,"CapsuleTipLog");
  CapsuleTipPhys1 = new G4PVPlacement(0,G4ThreeVector(0.,0.,1.84*mm),
                                      "CapsuleTipPhys1",CapsuleTipLog,mother,false,0);

  G4cout << "Created CapsuleTipPhys1" << G4endl;
  G4cout << "# of daughters - " << mother->GetLogicalVolume()->GetNoDaughters() << G4endl;

  // Rotated tip

  G4RotationMatrix* rotateMatrix = new G4RotationMatrix();
  rotateMatrix->rotateX(180.0*deg);
  CapsuleTipPhys2 = new G4PVPlacement(rotateMatrix, G4ThreeVector(0,0,-1.84*mm),
                                      "CapsuleTipPhys2",CapsuleTipLog,mother,false,0);

  G4cout << "Created CapsuleTipPhys2" << G4endl;
  G4cout << "# of daughters - " << mother->GetLogicalVolume()->GetNoDaughters() << G4endl;

  IodiumCore = new G4Tubs("ICore",0.085*mm,0.35*mm,1.75*mm,0.*deg,360.*deg);
  ICoreLog = new G4LogicalVolume(IodiumCore,Iodio,"IodiumCoreLog");
  ICorePhys = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),
                                "IodiumCorePhys",ICoreLog,DefaultTub_Phys,false,0);

  G4cout << "Created ICorePhys" << G4endl;
  G4cout << "# of daughters - " << mother->GetLogicalVolume()->GetNoDaughters() << G4endl;

  Marker=new G4Tubs("GoldenMarker",0.*mm,0.085*mm,1.75*mm,0.*deg,360.*deg);
  MarkerLog = new G4LogicalVolume(Marker,Gold,"MarkerLog");
  MarkerPhys=new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),
                               "MarkerPhys",MarkerLog,DefaultTub_Phys,false,0);

  G4cout << "Created MarkerPhys" << G4endl;
  G4cout << "# of daughters - " << mother->GetLogicalVolume()->GetNoDaughters() << G4endl;

  G4VisAttributes* simpleMarkerVisAtt= new G4VisAttributes(lblue);
  simpleMarkerVisAtt->SetVisibility(true);
  simpleMarkerVisAtt->SetForceSolid(true);

  G4VisAttributes* simpleIodiumVisAtt= new G4VisAttributes(magenta);
  simpleIodiumVisAtt->SetVisibility(true);
  simpleIodiumVisAtt->SetForceWireframe(true);
  ICoreLog->SetVisAttributes(simpleIodiumVisAtt);

  G4VisAttributes* simpleCapsuleVisAtt= new G4VisAttributes(red);
  simpleCapsuleVisAtt->SetVisibility(true);  
  simpleCapsuleVisAtt->SetForceWireframe(true);
  CapsuleLog->SetVisAttributes( simpleCapsuleVisAtt);

  G4VisAttributes* simpleCapsuleTipVisAtt= new G4VisAttributes(red);
  simpleCapsuleTipVisAtt->SetVisibility(true); 
  simpleCapsuleTipVisAtt->SetForceSolid(true);
  CapsuleTipLog->SetVisAttributes( simpleCapsuleTipVisAtt);
}
