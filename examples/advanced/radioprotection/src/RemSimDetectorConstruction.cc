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
//
// $Id: RemSimDetectorConstruction.cc,v 1.2 2004-02-03 09:16:46 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "RemSimDetectorConstruction.hh"
#include "RemSimDetectorMessenger.hh"
#include "RemSimMaterial.hh"
#include "RemSimVGeometryComponent.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4RunManager.hh"
#include "RemSimVehicle1.hh"
#include "RemSimVehicle2.hh"
#include "RemSimDecorator.hh"
#include "RemSimAstronautDecorator.hh"
RemSimDetectorConstruction::RemSimDetectorConstruction()
  :  experimentalHall_log(0),experimentalHall_phys(0),detectorChoice(0)
{
 pMaterial = new RemSimMaterial();
 pVehicle = new RemSimVehicle2();
 decorator = new RemSimAstronautDecorator(pVehicle);
 messenger = new RemSimDetectorMessenger(this);
}

RemSimDetectorConstruction::~RemSimDetectorConstruction()
{
  delete messenger;
  delete pVehicle;
  delete pMaterial;
}

G4VPhysicalVolume* RemSimDetectorConstruction::Construct()
{ 
  pMaterial -> DefineMaterials();
  G4Material* air = pMaterial->GetMaterial("Air") ;

  G4double expHall_x = 1000.0*m;
  G4double expHall_y = 1000.0*m;
  G4double expHall_z = 1000.0*m;
  G4Box* experimentalHall_box
    = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);
  experimentalHall_log = new G4LogicalVolume(experimentalHall_box,
                                             air,"expHall_log",0,0,0);
  experimentalHall_phys = new G4PVPlacement(0,G4ThreeVector(),
                                      experimentalHall_log,"expHall",0,false,0);
  ConstructVolume();

  return experimentalHall_phys;
}
void RemSimDetectorConstruction::ConstructVolume()
{
   //------------------------------ experimental hall (world volume)
  
  pVehicle ->ConstructComponent( experimentalHall_phys);
  
  //  decorator ->ConstructComponent( experimentalHall_phys);
}
void RemSimDetectorConstruction::SwitchVehicle(G4String value)
{
  
  pVehicle -> DestroyComponent();
  delete pVehicle;

  if(value=="Vehicle1") 
  {
   detectorChoice = 1;
  }
  else
  {
    if(value=="Vehicle2")
    {
      detectorChoice = 2;
    }
  }
   switch(detectorChoice)
  { 
    case 1:
      pVehicle = new RemSimVehicle1();
      break;
    case 2:
      pVehicle = new RemSimVehicle2();
      break;
    default:
      pVehicle = new RemSimVehicle2();
      break;
  }
  

  pVehicle ->ConstructComponent( experimentalHall_phys);
  
  if (detectorChoice == 2)
            decorator ->ConstructComponent( experimentalHall_phys);
  
  G4RunManager::GetRunManager()->DefineWorldVolume( experimentalHall_phys);
}  
G4double RemSimDetectorConstruction::GetDensity()
{
  G4double density = pVehicle -> GetDensity();
  return density;
}
void  RemSimDetectorConstruction::ChangeMaterial(G4String newMaterial)
{
  pVehicle -> ChangeMaterial(newMaterial);
}
