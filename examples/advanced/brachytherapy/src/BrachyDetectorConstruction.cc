
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
// Code developed by:
// S. Agostinelli, F. Foppiano, S. Garelli , M. Tropeano, S.Guatelli
//
//
//    ****************************************
//    *                                      *
//    *    BrachyDetectorConstruction.cc     *
//    *                                      *
//    ****************************************
//
// $Id: BrachyDetectorConstruction.cc,v 1.12 2002-11-18 15:18:37 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "BrachyPhantomROGeometry.hh"
#include "BrachyPhantomSD.hh"
#include "BrachyDetectorMessenger.hh"
#include "BrachyDetectorConstruction.hh"
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
#include "G4SDManager.hh"
#include "G4Colour.hh"
#include "G4UserLimits.hh"
#include "G4UnionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "BrachyMaterial.hh"
#include "BrachyFactoryLeipzig.hh"
#include "BrachyFactoryIr.hh"
#include "BrachyFactoryI.hh"

BrachyDetectorConstruction::BrachyDetectorConstruction(G4String &SDName):
NumVoxelX(0),NumVoxelZ(0),m_BoxDimX(0), m_BoxDimY(0), m_BoxDimZ(0),
ExpHall(0),ExpHallLog(0),ExpHallPhys(0),
Phantom(0), PhantomLog  (0), PhantomPhys  (0),
detectorChoice(0)

{
NumVoxelX=300;
NumVoxelZ=300;

m_BoxDimX=30*cm ;
m_BoxDimY=30*cm;
m_BoxDimZ=30*cm;

ComputeDimVoxel();

ExpHall_x = 4.0*m;
ExpHall_y = 4.0*m;
ExpHall_z = 4.0*m;

m_SDName = SDName;//pointer to sensitive detector

detectorMessenger = new BrachyDetectorMessenger(this);


G4Colour  white   (1.0, 1.0, 1.0) ;
G4Colour  grey    (0.5, 0.5, 0.5) ;
G4Colour  lgrey   (.75, .75, .75) ;
G4Colour  red     (1.0, 0.0, 0.0) ;
G4Colour  blue    (0.0, 0.0, 1.0) ;
G4Colour  cyan    (0.0, 1.0, 1.0) ;
G4Colour  magenta (1.0, 0.0, 1.0) ; 
G4Colour  yellow  (1.0, 1.0, 0.0) ;
G4Colour  lblue   (0.0, 0.0, .75);

factoryIr=new BrachyFactoryIr();

pMat= new BrachyMaterial();

}


BrachyDetectorConstruction::~BrachyDetectorConstruction()
{ 

delete  detectorMessenger;
delete pMat;

if (factoryLeipzig) delete factoryLeipzig;
if (factoryI) delete factoryI;
delete factoryIr;
}

//....
G4VPhysicalVolume* BrachyDetectorConstruction::Construct()
{ pMat-> DefineMaterials();

if (ExpHallPhys==0)
{
  ConstructPhantom(); }

switch(detectorChoice)
{ 
  case 1:
	 factoryI=new BrachyFactoryI();
	 factoryI->CreateSource(PhantomPhys);

   break;


  case 2:
         factoryLeipzig=new BrachyFactoryLeipzig();
         factoryLeipzig->CreateSource(PhantomPhys);
   break;

case 3:
       factoryIr->CreateSource(PhantomPhys);
   break;

default:

  factoryIr->CreateSource(PhantomPhys);
}


 ConstructSensitiveDetector();

return ExpHallPhys;

}
void BrachyDetectorConstruction::SwitchDetector()
{
factoryIr->CleanSource();
switch(detectorChoice)
{ 
case 1:
  factoryI=new BrachyFactoryI();
  factoryI->CreateSource(PhantomPhys);
   break;

case 2:
  factoryLeipzig=new BrachyFactoryLeipzig();
  factoryLeipzig->CreateSource(PhantomPhys);

  break;

case 3:
  factoryIr->CreateSource(PhantomPhys);

  break;   
default:
  factoryIr->CreateSource(PhantomPhys);
}

// Notify run manager that the new geometry has been built

G4RunManager::GetRunManager()->DefineWorldVolume( ExpHallPhys );
G4RunManager::GetRunManager()->GeometryHasBeenModified();

// Let the navigator to know about the new top of the new geometry
G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->SetWorldVolume(ExpHallPhys);
// We have to tweak the navigator's state. By the following dummy call we ensure that navigator's
// state is reset properly
G4ThreeVector center(0,0,0);
G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->LocateGlobalPointAndSetup(center,0,false);


}

void BrachyDetectorConstruction::SelectDetector(G4String val)
{ if(val=="Iodium") 
{ detectorChoice = 1; }
   else
     {if(val=="Leipzig") detectorChoice = 2;
     else { if(val=="Iridium")detectorChoice = 3;}}

   G4cout << "Now Detector is " << val << G4endl;
}
void BrachyDetectorConstruction::ConstructPhantom()

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

G4Material* air=pMat->GetMat("Air") ;
G4Material* soft=pMat->GetMat("tissue");

ComputeDimVoxel();

//World volume
ExpHall = new G4Box("ExpHall",ExpHall_x,ExpHall_y,ExpHall_z);
ExpHallLog = new G4LogicalVolume(ExpHall,air,"ExpHallLog",0,0,0);
ExpHallPhys = new G4PVPlacement(0,G4ThreeVector(),"ExpHallPhys",ExpHallLog,NULL,false,0);

// Water Box
Phantom = new G4Box("Phantom",m_BoxDimX/2,m_BoxDimY/2,m_BoxDimZ/2);
PhantomLog = new G4LogicalVolume(Phantom,soft,"PhantomLog",0,0,0);

PhantomPhys = new G4PVPlacement(0,G4ThreeVector(),"PhantomPhys",PhantomLog,ExpHallPhys,false,0); 

ExpHallLog->SetVisAttributes (G4VisAttributes::Invisible);

G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(lblue);
simpleBoxVisAtt->SetVisibility(true);
simpleBoxVisAtt->SetForceWireframe(true);

PhantomLog->SetVisAttributes(simpleBoxVisAtt);

}

void  BrachyDetectorConstruction::ConstructSensitiveDetector()
// Sensitive Detector and ReadOut geometry definition
{ 
G4SDManager* pSDManager = G4SDManager::GetSDMpointer();

if( pPhantomSD ==0)
 {
pPhantomSD = new BrachyPhantomSD(m_SDName,NumVoxelX,NumVoxelZ);
if(pPhantomSD)
{
      G4String ROGeometryName = "PhantomROGeometry";
        pPhantomROGeometry = new BrachyPhantomROGeometry(ROGeometryName,m_BoxDimX,m_BoxDimZ,NumVoxelX,NumVoxelZ);
      

      pPhantomROGeometry->BuildROGeometry();
      pPhantomSD->SetROgeometry(pPhantomROGeometry);
      pSDManager->AddNewDetector(pPhantomSD);

      PhantomLog->SetSensitiveDetector(pPhantomSD);


}
 }


}

void BrachyDetectorConstruction::PrintDetectorParameters()
{
G4cout << "-----------------------------------------------------------------------"
     << G4endl
    <<"the detector is a  box whose size is: " 
    <<G4endl
    <<m_BoxDimX/cm
    << " cm * "
    <<m_BoxDimY/cm
    << " cm * "
    <<m_BoxDimZ/cm
    << " cm"
    << G4endl
    <<"numVoxel: "
    <<NumVoxelX
    <<G4endl
    <<"dim voxel: "
    <<dimVoxel/mm
    <<"mm"
    <<G4endl 
    <<"material of the box : "
    <<AbsorberMaterial->GetName() 
    <<G4endl
    <<"the source is at the center  of the detector"
     <<G4endl


<<"-------------------------------------------------------------------------"
     << G4endl;
}


void BrachyDetectorConstruction::SetAbsorberMaterial(G4String materialChoice)

{
// search the material by its name   
G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
if (pttoMaterial)
 {AbsorberMaterial = pttoMaterial;
  PhantomLog->SetMaterial(pttoMaterial); 
  PrintDetectorParameters();
 } 
else G4cout<<"that's not avaiable!"<<G4endl;            
}

















