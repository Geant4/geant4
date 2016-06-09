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
// $Id: BrachyDetectorConstruction.cc,v 1.24 2004/03/11 16:05:02 guatelli Exp $
// GEANT4 tag $Name: geant4-06-01 $
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

BrachyDetectorConstruction::BrachyDetectorConstruction(G4String &SDName)
: detectorChoice(0), phantomSD(0), phantomROGeometry(0), factory(0),
  World(0), WorldLog(0), WorldPhys(0),
  Phantom(0), PhantomLog(0), PhantomPhys(0),
  phantomAbsorberMaterial(0)
{
  phantomDimensionX=15.*cm ;
  phantomDimensionY=15.*cm;
  phantomDimensionZ=15.*cm;

  numberOfVoxelsAlongX=300;
  numberOfVoxelsAlongZ=300;

  ComputeDimVoxel();

  Worldx = 4.0*m;
  Worldy = 4.0*m;
  Worldz = 4.0*m;

  sensitiveDetectorName = SDName;

  detectorMessenger = new BrachyDetectorMessenger(this);

  factory = new BrachyFactoryIr();

  pMaterial= new BrachyMaterial();
}


BrachyDetectorConstruction::~BrachyDetectorConstruction()
{ 
  delete pMaterial;
  delete factory;
  delete detectorMessenger;
  if (phantomROGeometry) delete phantomROGeometry;
}

G4VPhysicalVolume* BrachyDetectorConstruction::Construct()
{
  pMaterial-> DefineMaterials();
  ConstructPhantom();
  factory->CreateSource(PhantomPhys);
  ConstructSensitiveDetector();

  return WorldPhys;
}

void BrachyDetectorConstruction::SwitchBrachytherapicSeed()
{
  factory->CleanSource();
  delete factory;
  switch(detectorChoice)
  { 
    case 1:
      factory = new BrachyFactoryI();
      break;
    case 2:
      factory = new BrachyFactoryLeipzig();
      break;
    case 3:
      factory = new BrachyFactoryIr();
      break;   
    default:
      factory = new BrachyFactoryIr();
      break;
  }
  factory->CreateSource(PhantomPhys);

  // Notify run manager that the new geometry has been built

  G4RunManager::GetRunManager()->DefineWorldVolume( WorldPhys );
}

void BrachyDetectorConstruction::SelectBrachytherapicSeed(G4String val)
{
  if(val=="Iodium") 
  {
   detectorChoice = 1;
  }
  else
  {
    if(val=="Leipzig")
    {
      detectorChoice = 2;
    }
    else
    {
      if(val=="Iridium")
      {
        detectorChoice = 3;
      }
    }
  }
  G4cout << "Now Detector is " << val << G4endl;
}

void BrachyDetectorConstruction::ConstructPhantom()
{
  G4Colour  lblue   (0.0, 0.0, .75);

  G4Material* air=pMaterial->GetMat("Air") ;
  G4Material* water=pMaterial->GetMat("Water");

  ComputeDimVoxel();

  // World volume
  World = new G4Box("World",Worldx,Worldy,Worldz);
  WorldLog = new G4LogicalVolume(World,air,"WorldLog",0,0,0);
  WorldPhys = new G4PVPlacement(0,G4ThreeVector(),"WorldPhys",WorldLog,0,false,0);

  // Water Box
  Phantom = new G4Box("Phantom",phantomDimensionX,phantomDimensionY,phantomDimensionZ);
  PhantomLog = new G4LogicalVolume(Phantom,water,"PhantomLog",0,0,0);
  PhantomPhys = new G4PVPlacement(0,G4ThreeVector(),"PhantomPhys",PhantomLog,WorldPhys,false,0); 

  WorldLog->SetVisAttributes (G4VisAttributes::Invisible);

  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(lblue);
  simpleBoxVisAtt->SetVisibility(true);
  simpleBoxVisAtt->SetForceWireframe(true);

  PhantomLog->SetVisAttributes(simpleBoxVisAtt);
}

void  BrachyDetectorConstruction::ConstructSensitiveDetector()
// Sensitive Detector and ReadOut geometry definition
{ 
  G4SDManager* pSDManager = G4SDManager::GetSDMpointer();

  if(!phantomSD)
  {
    phantomSD = new BrachyPhantomSD(sensitiveDetectorName);
    G4String ROGeometryName = "PhantomROGeometry";
    phantomROGeometry = new BrachyPhantomROGeometry(ROGeometryName,phantomDimensionX,phantomDimensionZ,numberOfVoxelsAlongX,numberOfVoxelsAlongZ);
    phantomROGeometry->BuildROGeometry();
    phantomSD->SetROgeometry(phantomROGeometry);
    pSDManager->AddNewDetector(phantomSD);
    PhantomLog->SetSensitiveDetector(phantomSD);
  }
}

void BrachyDetectorConstruction::PrintDetectorParameters()
{
  G4cout << "-----------------------------------------------------------------------"
         << G4endl
         << "the detector is a  box whose size is: " << G4endl
         << phantomDimensionX/cm
         << " cm * "
         << phantomDimensionY/cm
         << " cm * "
         << phantomDimensionZ/cm
         << " cm" << G4endl
         << "numVoxel: "
         << numberOfVoxelsAlongX <<G4endl
         << "dim voxel: "
         << dimVoxel/mm
         << "mm" << G4endl 
         << "material of the box : "
         << phantomAbsorberMaterial->GetName() <<G4endl
         << "the source is at the center  of the detector" << G4endl
         << "-------------------------------------------------------------------------"
         << G4endl;
}


void BrachyDetectorConstruction::SetPhantomMaterial(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial)
  {
    phantomAbsorberMaterial = pttoMaterial;
    PhantomLog->SetMaterial(pttoMaterial); 
    PrintDetectorParameters();
  } 
  else
    G4cout << "WARNING: material '" << materialChoice
           << "' not available!" << G4endl;            
}
