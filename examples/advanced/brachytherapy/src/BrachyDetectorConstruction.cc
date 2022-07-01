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
#include "BrachyFactoryLeipzig.hh"
#include "BrachyFactoryTG186.hh"
#include "BrachyFactoryI.hh"
#include "BrachyFactoryFlexi.hh"
#include "BrachyFactoryOncura6711.hh"
#include "BrachyDetectorMessenger.hh"
#include "BrachyDetectorConstruction.hh"
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
#include "G4NistManager.hh"

BrachyDetectorConstruction::BrachyDetectorConstruction(): 
  fFactory(nullptr), fWorld(nullptr), fWorldLog(nullptr), fWorldPhys(nullptr),
  fPhantom(nullptr), fPhantomLog(nullptr), fPhantomPhys(nullptr),
  fDetectorChoice(0)
{
 // Define the messenger of the Detector component
 fDetectorMessenger = new BrachyDetectorMessenger(this);

 // Define half size of the phantom along the x, y, z axis
 fPhantomSizeX = 15.*cm;
 fPhantomSizeY = 15.*cm;
 fPhantomSizeZ = 15.*cm;
 
 // Define the sizes of the World volume containing the phantom
 fWorldSizeX = 4.0*m;
 fWorldSizeY = 4.0*m;
 fWorldSizeZ = 4.0*m;

  // Define the Flexi source as default source modelled in the geometry
 fFactory = new BrachyFactoryFlexi();
}

BrachyDetectorConstruction::~BrachyDetectorConstruction()
{ 
  delete fDetectorMessenger; 
  delete fFactory;
}

G4VPhysicalVolume* BrachyDetectorConstruction::Construct()
{
 // Model the phantom (water box)
 ConstructPhantom();

 // Model the source in the phantom
 fFactory -> CreateSource(fPhantomPhys);

 return fWorldPhys;
}

void BrachyDetectorConstruction::SwitchBrachytherapicSeed()
{
  // Change the source in the water phantom
  fFactory -> CleanSource();
  G4cout << "Old brachy source is deleted ..." << G4endl;
  delete fFactory;

  switch(fDetectorChoice)
  { 
   case 1:
      fFactory = new BrachyFactoryI();
      break;
   case 2:
      fFactory = new BrachyFactoryLeipzig();
      break;
   case 3:
      fFactory = new BrachyFactoryTG186();
      break;
   case 4:
      fFactory = new BrachyFactoryFlexi();
      break;
   case 5:
   	  fFactory = new BrachyFactoryOncura6711();
   	  break;
   default:
      fFactory = new BrachyFactoryFlexi();
      break;
   }

  fFactory -> CreateSource(fPhantomPhys);
  G4cout << "New brachy source is created ..." << G4endl;

  // Notify run manager that the new geometry has been built
  G4RunManager::GetRunManager() -> GeometryHasBeenModified();
}

void BrachyDetectorConstruction::SelectBrachytherapicSeed(G4String val)
{
	if (val == "Iodine") 
		fDetectorChoice = 1;
  		else
  		{
       		if(val=="Leipzig") 
       			fDetectorChoice = 2;
       			else
       			{
	            	if(val=="TG186") 
	            		fDetectorChoice = 3;
	             		else
             			{ 
		  					if(val=="Flexi") 
		  						fDetectorChoice = 4;
		  						else
		  						{
		  							if(val=="Oncura") 
		  								 fDetectorChoice = 5;
		  		                  		else 
		  		                  			G4cout << val <<  "is not available!!!!"  << G4endl;  
		  						}            
            			}
       	}		
  }
 G4cout << "Now the brachy source is " << val << G4endl;
}

void BrachyDetectorConstruction::ConstructPhantom()
{
  // Model the water phantom 

  // Define the light blue color
  G4Colour  lblue   (0.0, 0.0, .75);
  
   //Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* air = nist -> FindOrBuildMaterial("G4_AIR");
  G4Material* water = nist -> FindOrBuildMaterial("G4_WATER");

  // World volume
  fWorld = new G4Box("World", fWorldSizeX, fWorldSizeY, fWorldSizeZ);
  fWorldLog = new G4LogicalVolume(fWorld,air,"WorldLog",nullptr,nullptr,nullptr);
  fWorldPhys = new G4PVPlacement(nullptr,G4ThreeVector(),"WorldPhys",fWorldLog,nullptr,false,0);

  // Water Box
  fPhantom = new G4Box("Phantom", fPhantomSizeX, fPhantomSizeY, fPhantomSizeZ);

  // Logical volume
  fPhantomLog = new G4LogicalVolume(fPhantom,water,"PhantomLog",nullptr,nullptr,nullptr);

  // Physical volume
  fPhantomPhys = new G4PVPlacement(nullptr,G4ThreeVector(), // Position: rotation and translation
                                  "PhantomPhys", // Name
                                  fPhantomLog, // Associated logical volume
                                  fWorldPhys, // Mother volume
                                  false,0);

  fWorldLog -> SetVisAttributes (G4VisAttributes::GetInvisible());

  // Visualization attributes of the phantom
  auto simpleBoxVisAtt = new G4VisAttributes(lblue);
  simpleBoxVisAtt -> SetVisibility(true);
  simpleBoxVisAtt -> SetForceWireframe(true);
  fPhantomLog -> SetVisAttributes(simpleBoxVisAtt);
}

void BrachyDetectorConstruction::PrintDetectorParameters()
{
  G4cout << "----------------" << G4endl
         << "the phantom is a water box whose size is: " << G4endl
         << fPhantomSizeX *2./cm
         << " cm * "
         << fPhantomSizeY *2./cm
         << " cm * "
         << fPhantomSizeZ *2./cm
         << " cm" << G4endl
         << "The phantom is made of "
         << fPhantomLog -> GetMaterial() -> GetName() <<G4endl
         << "the source is at the center of the phantom" << G4endl
         << "----------------"
         << G4endl;
}

void BrachyDetectorConstruction::SetPhantomMaterial(G4String materialChoice)
{
  // Method to change the material of the phantom
 
  G4NistManager* nist = G4NistManager::Instance();

  // Search the material by its name
  G4Material* pttoMaterial =  nist -> FindOrBuildMaterial(materialChoice);

  if (pttoMaterial)
  {
    fPhantomLog -> SetMaterial(pttoMaterial);
    PrintDetectorParameters();
  } else  G4cout << "WARNING: material '" << materialChoice << "' not available!" << G4endl;
}
