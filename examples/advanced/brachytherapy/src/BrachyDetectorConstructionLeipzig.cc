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
// ********************************************************************

//
// --------------------------------------------------------------
//                 GEANT 4 - Brachytherapy example
// --------------------------------------------------------------
//
// Code developed by:
// S. Agostinelli, F. Foppiano, S. Garelli , M. Tropeano, S.Guatelli
//
//    *******************************************
//    *                                         *
//    *    BrachyDetectorConstructionLeipzig.cc *
//    *                                         *
//    *******************************************
//
//
// $Id: BrachyDetectorConstructionLeipzig.cc,v 1.3 2002-12-12 19:16:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
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
#include "globals.hh"
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
BrachyDetectorConstructionLeipzig::BrachyDetectorConstructionLeipzig():
Capsule(0),CapsuleLog(0),CapsulePhys(0),
CapsuleTip(0),CapsuleTipLog(0),CapsuleTipPhys(0),
IridiumCore(0),IridiumCoreLog(0),IridiumCorePhys(0),
Appl1(0),Appl1Log(0),Appl1Phys(0),
Appl2(0),Appl2Log(0),Appl2Phys(0)

{ 
pMat= new BrachyMaterial() ;

}


//....

BrachyDetectorConstructionLeipzig::~BrachyDetectorConstructionLeipzig()
{ 

delete pMat; 
}
void  BrachyDetectorConstructionLeipzig::ConstructLeipzig(G4VPhysicalVolume*   WaterBoxPhys)
{

G4Colour  white   (1.0, 1.0, 1.0) ;
G4Colour  grey    (0.5, 0.5, 0.5) ;
G4Colour  lgrey   (.75, .75, .75) ;
G4Colour  red     (1.0, 0.0, 0.0) ;
G4Colour  blue    (0.0, 0.0, 1.0) ;
G4Colour  cyan    (0.0, 1.0, 1.0) ;
G4Colour  magenta (1.0, 0.0, 1.0) ; 
G4Colour  yellow  (1.0, 1.0, 0.0) ;
G4Colour  lblue   (0.0, 0.0, .75);

G4Material* CapsuleMat=pMat->GetMat("Stainless steel");
G4Material*IridiumMat=pMat->GetMat("Iridium");
G4Material* matTung=pMat->GetMat("Tungsten");

Capsule = new G4Tubs("Capsule",0,0.55*mm,3.725*mm,0.*deg,360.*deg);
CapsuleLog = new G4LogicalVolume(Capsule,CapsuleMat,"CapsuleLog");
CapsulePhys = new G4PVPlacement(0,G4ThreeVector(0,0,-1.975*mm),"CapsulePhys",
				    CapsuleLog,WaterBoxPhys,false,0);

// Capsule tip

CapsuleTip = new G4Sphere("CapsuleTip",
				 0.*mm,0.55*mm,0.*deg,360.*deg,0.*deg,90.*deg);
CapsuleTipLog = new G4LogicalVolume(CapsuleTip,CapsuleMat,"CapsuleTipLog");
CapsuleTipPhys = new G4PVPlacement(0,G4ThreeVector(0.,0.,1.75*mm),"CapsuleTipPhys",
				    CapsuleTipLog,WaterBoxPhys,false,0);

// Iridium core

IridiumCore = new G4Tubs("IrCore",0,0.30*mm,1.75*mm,0.*deg,360.*deg);
IridiumCoreLog = new G4LogicalVolume(IridiumCore,IridiumMat,"IridiumCoreLog");
IridiumCorePhys = new G4PVPlacement(0,G4ThreeVector(0.,0.,1.975*mm),"IridiumCorePhys",
				    IridiumCoreLog,CapsulePhys,false,0);
// Applicator

Appl1 = new G4Tubs("Appl1",5*mm,10.5*mm,12*mm,0.*deg,360.*deg);
Appl1Log = new G4LogicalVolume(Appl1,matTung,"Appl1");
Appl1Phys = new G4PVPlacement (0,G4ThreeVector(0,0,4.0*mm),"Appl1Phys", Appl1Log,
				     WaterBoxPhys,false,0);

Appl2 = new G4Tubs("Appl2",0.55*mm,5.*mm,3.125*mm,0.*deg,360.*deg);
Appl2Log = new G4LogicalVolume(Appl2,matTung,"Appl2");
Appl2Phys = new G4PVPlacement(0,G4ThreeVector(0,0,-4.875*mm),"Appl2Phys",Appl2Log,
				     WaterBoxPhys,false,0);


simpleCapsuleVisAtt= new G4VisAttributes(red); 
simpleCapsuleVisAtt->SetVisibility(true); 
simpleCapsuleVisAtt->SetForceSolid(true); 
CapsuleLog->SetVisAttributes(simpleCapsuleVisAtt); 


//CapsuleTipLog->SetVisAttributes (G4VisAttributes::Invisible); 
simpleCapsuleTipVisAtt= new G4VisAttributes(red); 
simpleCapsuleTipVisAtt->SetVisibility(true); 
simpleCapsuleTipVisAtt->SetForceSolid(true); 
CapsuleTipLog->SetVisAttributes(simpleCapsuleTipVisAtt); 



}




