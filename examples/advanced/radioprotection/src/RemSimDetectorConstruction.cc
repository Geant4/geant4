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
//    **************************************
//    *                                    *
//    *    RemSimDetectorConstruction.cc   *
//    *                                    *          
//    **************************************
//
// $Id: RemSimDetectorConstruction.cc,v 1.17 2006-06-29 16:23:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author:Susanna Guatelli, guatelli@ge.infn.it 
//
#include "RemSimDetectorConstruction.hh"
#include "RemSimMaterial.hh"
#include "RemSimDetectorMessenger.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4RunManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4UserLimits.hh"
#include "RemSimVGeometryComponent.hh"
#include "RemSimVehicle1.hh"
#include "RemSimMoonHabitat.hh"
#include "G4VisAttributes.hh"
#include "RemSimSensitiveDetector.hh"
#include "G4SDManager.hh"
#include "RemSimROGeometry.hh"
#include "RemSimDecorator.hh"
#include "RemSimShieldingDecorator.hh"
#include "RemSimShelterSPEDecorator.hh"
#include "RemSimAstronautDecorator.hh"
#include "RemSimRoofDecorator.hh"

RemSimDetectorConstruction::RemSimDetectorConstruction()
  :  experimentalHall_log(0), experimentalHall_phys(0),
     phantomPhys(0), detectorPhys(0), geometry(0)
{
 pMaterial = new RemSimMaterial();

 messenger = new RemSimDetectorMessenger(this);
 decoratorValue ="Nothing";
 astronautValue ="On";
 decorator = 0;  //pointing to shielding
 decoratorSPE = 0; //pointing to SPE shelter
 decorator1 = 0; //pointing to the phantom/astronaut
 decoratorRoof = 0; //pointing to the roof of the Moon Habitat
 moon = false;  // moon = true: moon configuration selected
                // moon = false: vehicle configuration selected
 flag = false;  // flag = true : the geometry configuration (vehicle/moon) 
                //                has been chosen 
                // flag = false : the geometry configuration (vehicle/moon) 
                //                has not been chosen}
}

RemSimDetectorConstruction::~RemSimDetectorConstruction()
{  
  delete decoratorRoof;
  delete decorator1;
  delete decoratorSPE;
  delete decorator;
  delete messenger;
  delete geometry;
  delete pMaterial;
}

G4VPhysicalVolume* RemSimDetectorConstruction::Construct()
{ 
  pMaterial -> DefineMaterials();
  G4Material* vacuum = pMaterial -> GetMaterial("Galactic");
 
  // World definition

  G4double expHall_x = 30.0*m;
  G4double expHall_y = 30.0*m;
  G4double expHall_z = 30.0*m;
  
  G4Box* experimentalHall_box
    = new G4Box("world",expHall_x,expHall_y,expHall_z);
  experimentalHall_log = new G4LogicalVolume(experimentalHall_box,
                                             vacuum,"world",0,0,0);
  
  experimentalHall_phys = new G4PVPlacement(0,G4ThreeVector(),"world",
                                      experimentalHall_log,0,
                                      false,0);
  return experimentalHall_phys;
}

void RemSimDetectorConstruction::ConstructVolume()
{ 
  if (moon == true)
    { 

      geometry = new RemSimMoonHabitat();
      geometry -> ConstructComponent(experimentalHall_phys);
      decorator1 = new RemSimAstronautDecorator(geometry, moon);
      decorator1 -> ConstructComponent(experimentalHall_phys);
    }
  else 
    {
      geometry = new RemSimVehicle1();
      geometry -> ConstructComponent(experimentalHall_phys); 
      decorator1 = new RemSimAstronautDecorator(geometry, moon); 
      decorator1 -> ConstructComponent(experimentalHall_phys);
    }
}

void RemSimDetectorConstruction::AddShielding(G4String value)
{ 
  if (moon == false)
    {
  decoratorValue = value;
   
  if (decoratorValue == "On")
	{
	  if (decorator == 0)
	    { 
	      decorator = new RemSimShieldingDecorator(geometry);
	      decorator -> ConstructComponent(experimentalHall_phys); 
	      G4RunManager::GetRunManager() -> DefineWorldVolume(experimentalHall_phys);
	    }
	  else  G4cout<<" The shielding alread exists!"<<G4endl;
	}

      if (decoratorValue == "Off")
	{
	  if (decorator != 0)
	    {
	  decorator -> DestroyComponent(); 
	  G4RunManager::GetRunManager() -> DefineWorldVolume(experimentalHall_phys);
          decorator = 0;
	   }
         
          else  G4cout<<" The shielding does not exist!"<<G4endl;
	}
    }
  else G4cout << " The shielding is not available in moon habitat configuration"<<G4endl;
}

void RemSimDetectorConstruction::AddShelterSPE(G4String value)
{ 
  if (moon == false)
    { 
      if (value == "On")
	{
	if (decoratorSPE == 0)
	  { 
	    decoratorSPE = new RemSimShelterSPEDecorator(geometry);
	    decoratorSPE -> ConstructComponent(experimentalHall_phys); 
	    G4RunManager::GetRunManager() -> DefineWorldVolume(experimentalHall_phys);
	  }
	else  G4cout<<" The SPE shelter alread exists!"<<G4endl;
	}

  if (value == "Off")
    {
      if (decoratorSPE != 0)
	{
	  decoratorSPE -> DestroyComponent(); 
	  G4RunManager::GetRunManager() -> DefineWorldVolume(experimentalHall_phys);
          decoratorSPE = 0;
	}
         
      else  G4cout<<" The SPE shelter does not exist!"<<G4endl;
    }
    }
  else  G4cout<< " It is not possible to select SPE shelter in moon habitat configuration" << G4endl;
}

void RemSimDetectorConstruction::AddHabitatRoof(G4String value)
{ 
  if (moon == true)
    { 
      if (value == "On")
      	{
	if (decoratorRoof == 0)
	  { 
	    decoratorRoof = new RemSimRoofDecorator(geometry);
	    decoratorRoof -> ConstructComponent(experimentalHall_phys); 
	    G4RunManager::GetRunManager() -> DefineWorldVolume(experimentalHall_phys);
	  }
		else  G4cout<<" The roof alread exists!"<<G4endl;
		}

  if (value == "Off")
    {
      if (decoratorRoof != 0)
	{
	  decoratorRoof -> DestroyComponent(); 
	  G4RunManager::GetRunManager() -> DefineWorldVolume(experimentalHall_phys);
          decoratorRoof = 0;
	}
         
      else  G4cout<<" The roof does not exist!"<<G4endl;
    }
    }
  else  G4cout<< " It is not possible to select the roof in the vehicle configuration" << G4endl;
}

void RemSimDetectorConstruction::ChangeShieldingThickness(G4double thick)
{  
  if (moon == false)
    {
      if (decorator != 0)
	{ 
	  decorator -> ChangeThickness(thick);
	  G4RunManager::GetRunManager() -> DefineWorldVolume(experimentalHall_phys);
	}
      else G4cout<<" Define the shielding before!"<< G4endl;
    }
}

void RemSimDetectorConstruction::ChangeRoofThickness(G4double thick)
{  
  if (moon == true)
    {
      if (decoratorRoof != 0)
	{ 
	  decoratorRoof -> ChangeThickness(thick);
	  G4RunManager::GetRunManager() -> DefineWorldVolume(experimentalHall_phys);
	}
      else G4cout<<" Define the roof before!"<< G4endl;
    }

  else G4cout<<"The roof can be selected in moon habitat configuration"<<G4endl;
}

void RemSimDetectorConstruction:: ChooseConfiguration(G4String value)
{
  if (flag == false)
    { 
      if (value == "vehicle") ConstructVolume();
  
      else if (value == "moon")
      	{
	  moon = true;
	  ConstructVolume();
	  	}
      G4RunManager::GetRunManager() -> DefineWorldVolume(experimentalHall_phys);
       flag = true;
      }
       else 
       G4cout<< "The configurations vehicle/moon can not be switched" << G4endl;
}

