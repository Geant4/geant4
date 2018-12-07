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
// Code developed by: S. Guatelli, D. Cutajar, A. Le
// Past developers: S. Agostinelli, F. Foppiano, S. Garelli , M. Tropeano
//
//
//    ****************************************
//    *                                      *
//    *    BrachyDetectorConstruction.cc     *
//    *                                      *
//    ****************************************
//
#include "G4SystemOfUnits.hh"
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
#include "BrachyFactoryTG186.hh"
#include "BrachyFactoryI.hh"
#include "BrachyFactoryFlexi.hh"
#include "BrachyFactoryOncura6711.hh"
#include "BrachyDetectorMessenger.hh"
#include "BrachyDetectorConstruction.hh"

BrachyDetectorConstruction::BrachyDetectorConstruction(): 
  detectorChoice(0), factory(0),
  World(0), WorldLog(0), WorldPhys(0),
  Phantom(0), PhantomLog(0), PhantomPhys(0),
  phantomAbsorberMaterial(0)
{
 // Define half size of the phantom along the x, y, z axis
 phantomSizeX = 15.*cm;
 phantomSizeY = 15.*cm;
 phantomSizeZ = 15.*cm;
 
 // Define the sizes of the World volume containing the phantom
 worldSizeX = 4.0*m;
 worldSizeY = 4.0*m;
 worldSizeZ = 4.0*m;

 // Define the messenger of the Detector component
 // It is possible to modify geometrical parameters through UI
 detectorMessenger = new BrachyDetectorMessenger(this);

  // Define the Flexi source as default source modelled in the geometry
 factory = new BrachyFactoryFlexi();

  // BrachyMaterial defined the all the materials necessary
  // for the experimental set-up 
 pMaterial = new BrachyMaterial();
}

BrachyDetectorConstruction::~BrachyDetectorConstruction()
{ 
  delete pMaterial;
  delete factory;
  delete detectorMessenger;
}

G4VPhysicalVolume* BrachyDetectorConstruction::Construct()
{
 pMaterial -> DefineMaterials();

 // Model the phantom (water box)
 ConstructPhantom();

 // Model the source in the phantom
 factory -> CreateSource(PhantomPhys);

 return WorldPhys;
}

void BrachyDetectorConstruction::SwitchBrachytherapicSeed()
{
  // Change the source in the water phantom
  factory -> CleanSource();
  G4cout << "Old Source is deleted ..." << G4endl;
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
      factory = new BrachyFactoryTG186();
      break;
   case 4:
      factory = new BrachyFactoryFlexi();
      break;
   case 5:
   	  factory = new BrachyFactoryOncura6711();
   	  break;
   default:
      factory = new BrachyFactoryFlexi();
      break;
   }

  factory -> CreateSource(PhantomPhys);
  G4cout << "... New source is created ..." << G4endl;

  // Notify run manager that the new geometry has been built
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
  G4cout << "... Geometry is notified .... THAT'S IT!!!!!" << G4endl;
}

void BrachyDetectorConstruction::SelectBrachytherapicSeed(G4String val)
{
	if (val == "Iodine") 
		detectorChoice = 1;
  		else
  		{
       		if(val=="Leipzig") 
       			detectorChoice = 2;
       			else
       			{
	            	if(val=="TG186") 
	            		detectorChoice = 3;
	             		else
             			{ 
		  					if(val=="Flexi") 
		  						detectorChoice = 4;
		  						else
		  						{
		  							if(val=="Oncura") 
		  								detectorChoice = 5;
		  		                  		else 
		  		                  			G4cout << val <<  "is not available!!!!"  << G4endl;  
		  						}            
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

  // World volume
  World = new G4Box("World",worldSizeX,worldSizeY,worldSizeZ);
  WorldLog = new G4LogicalVolume(World,air,"WorldLog",0,0,0);
  WorldPhys = new G4PVPlacement(0,G4ThreeVector(),"WorldPhys",WorldLog,0,false,0);

  // Water Box
  Phantom = new G4Box("Phantom",phantomSizeX,phantomSizeY,phantomSizeZ);

  // Logical volume
  PhantomLog = new G4LogicalVolume(Phantom,water,"PhantomLog",0,0,0);

  // Physical volume
  PhantomPhys = new G4PVPlacement(0,G4ThreeVector(), // Position: rotation and translation
                                  "PhantomPhys", // Name
                                  PhantomLog, // Associated logical volume
                                  WorldPhys, // Mother volume
                                  false,0);
  WorldLog -> SetVisAttributes (G4VisAttributes::GetInvisible());

  // Visualization attributes of the phantom
  G4VisAttributes* simpleBoxVisAtt = new G4VisAttributes(lblue);
  simpleBoxVisAtt -> SetVisibility(true);
  simpleBoxVisAtt -> SetForceWireframe(true);
  PhantomLog -> SetVisAttributes(simpleBoxVisAtt);
}

void BrachyDetectorConstruction::PrintDetectorParameters()
{
  G4cout << "----------------" << G4endl
         << "the phantom is a water box whose size is: " << G4endl
         << phantomSizeX *2./cm
         << " cm * "
         << phantomSizeY *2./cm
         << " cm * "
         << phantomSizeZ *2./cm
         << " cm" << G4endl
         << "The phantom is made of "
         << phantomAbsorberMaterial -> GetName() <<G4endl
         << "the source is at the center of the phantom" << G4endl
         << "----------------"
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
  } else 
    { G4cout << "WARNING: material '" << materialChoice << "' not available!" << G4endl;}
}
