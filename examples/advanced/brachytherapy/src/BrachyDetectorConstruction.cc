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
#include "G4CSGSolid.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4MaterialTable.hh"
#include "G4TransportationManager.hh"
#include "G4Colour.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "BrachyMaterial.hh"
#include "BrachyFactoryLeipzig.hh"
#include "BrachyFactoryIr.hh"
#include "BrachyFactoryI.hh"
#include "BrachyPhantomROGeometry.hh"
#include "BrachyPhantomSD.hh"
#include "BrachyDetectorMessenger.hh"
#include "BrachyDetectorConstruction.hh"

BrachyDetectorConstruction::BrachyDetectorConstruction(G4String &SDName)
: detectorChoice(0), phantomSD(0), phantomROGeometry(0), factory(0),
  World(0), WorldLog(0), WorldPhys(0),
  Phantom(0), PhantomLog(0), PhantomPhys(0),
  phantomAbsorberMaterial(0)
{
  // Define half size of the phantom along the x, y, z axis
  phantomSizeX = 15.*cm ;
  phantomSizeY = 15.*cm;
  phantomSizeZ = 15.*cm;

  // Define the number of voxels the phantom is subdivided along the
  // three axis 
  // Each voxel is 1 mm wide
  numberOfVoxelsAlongX = 300;
  numberOfVoxelsAlongY = 300;
  numberOfVoxelsAlongZ = 300;

  ComputeDimVoxel();
  
  // Define the sizes of the World volume contaning the phantom
  worldSizeX = 4.0*m;
  worldSizeY = 4.0*m;
  worldSizeZ = 4.0*m;

  sensitiveDetectorName = SDName;

  // Define the messenger of the Detector component
  // It is possible to modify geometrical parameters through UI
  detectorMessenger = new BrachyDetectorMessenger(this);

  // Define the Iridium source as default source modelled in the geometry
  factory = new BrachyFactoryIr();

  // BrachyMaterial defined the all the materials necessary
  // for the experimental set-up 
  pMaterial = new BrachyMaterial();
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
  pMaterial -> DefineMaterials();
 
  // Model the phantom (water box)  
  ConstructPhantom();

  // Model the source in the phantom
  factory -> CreateSource(PhantomPhys); //Build the source inside the phantom
 
  // Define the sensitive volume: phantom
  ConstructSensitiveDetector();

  return WorldPhys;
}

void BrachyDetectorConstruction::SwitchBrachytherapicSeed()
{
  // Change the source in the water phantom
  factory -> CleanSource();
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

  factory -> CreateSource(PhantomPhys);

  // Notify run manager that the new geometry has been built
  G4RunManager::GetRunManager() -> DefineWorldVolume( WorldPhys );
}

void BrachyDetectorConstruction::SelectBrachytherapicSeed(G4String val)
{
  if(val == "Iodium") 
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

  G4cout << "Now the source is " << val << G4endl;
}

void BrachyDetectorConstruction::ConstructPhantom()
{
  // Model the water phantom 
  
  // Define the light blue color
  G4Colour  lblue   (0.0, 0.0, .75);

  G4Material* air = pMaterial -> GetMat("Air") ;
  G4Material* water = pMaterial -> GetMat("Water");

  ComputeDimVoxel();

  // World volume
  World = new G4Box("World",worldSizeX,worldSizeY,worldSizeZ);
  WorldLog = new G4LogicalVolume(World,air,"WorldLog",0,0,0);
  WorldPhys = new G4PVPlacement(0,G4ThreeVector(),
                                "WorldPhys",WorldLog,0,false,0);

  // Water Box
  Phantom = new G4Box("Phantom",phantomSizeX,phantomSizeY,
                       phantomSizeZ);

  // Logical volume
  PhantomLog = new G4LogicalVolume(Phantom,water,"PhantomLog",0,0,0);

  // Physical volume
  PhantomPhys = new G4PVPlacement(0,G4ThreeVector(), // Position: rotation and translation
                                  "PhantomPhys", // Name
				  PhantomLog, // Associated logical volume
                                  WorldPhys, // Mother volume
				  false,0); 

  WorldLog -> SetVisAttributes (G4VisAttributes::Invisible);

  // Visualization attributes of the phantom
  G4VisAttributes* simpleBoxVisAtt = new G4VisAttributes(lblue);
  simpleBoxVisAtt -> SetVisibility(true);
  simpleBoxVisAtt -> SetForceWireframe(true);
  PhantomLog -> SetVisAttributes(simpleBoxVisAtt);
}

void  BrachyDetectorConstruction::ConstructSensitiveDetector()
// Sensitive Detector and ReadOut geometry definition
{ 
  G4SDManager* pSDManager = G4SDManager::GetSDMpointer();

  if(!phantomSD)
  {
    phantomSD = new BrachyPhantomSD(sensitiveDetectorName);
    G4String ROGeometryName = "PhantomROGeometry";
    phantomROGeometry = new BrachyPhantomROGeometry(ROGeometryName,
                                phantomSizeX,  
				phantomSizeY, 
				phantomSizeZ,
                                numberOfVoxelsAlongX,
				numberOfVoxelsAlongY,		    
				numberOfVoxelsAlongZ);
    phantomROGeometry -> BuildROGeometry();
    phantomSD -> SetROgeometry(phantomROGeometry);
    pSDManager -> AddNewDetector(phantomSD);
   
    PhantomLog -> SetSensitiveDetector(phantomSD);
  }
}

void BrachyDetectorConstruction::PrintDetectorParameters()
{
  G4cout << "-----------------------------------------------------------------------"
         << G4endl
         << "the phantom is a water box whose size is: " << G4endl
         << phantomSizeX *2./cm
         << " cm * "
         << phantomSizeY *2./cm
         << " cm * "
         << phantomSizeZ *2./cm
         << " cm" << G4endl
         << "number of Voxel: "
         << numberOfVoxelsAlongX <<G4endl
         << "Voxel size: "
         << dimVoxel * 2/mm
         << "mm" << G4endl 
         << "The phantom is made of "
         << phantomAbsorberMaterial -> GetName() <<G4endl
         << "the source is at the center of the phantom" << G4endl
         << "-------------------------------------------------------------------------"
         << G4endl;
}

void BrachyDetectorConstruction::SetPhantomMaterial(G4String materialChoice)
{
  // It is possible to change the material of the phantom
  // interactively

  // Search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial)
  {
    phantomAbsorberMaterial = pttoMaterial;
    PhantomLog -> SetMaterial(pttoMaterial); 
    PrintDetectorParameters();
  } 
  else
    G4cout << "WARNING: material '" << materialChoice
           << "' not available!" << G4endl;            
}
