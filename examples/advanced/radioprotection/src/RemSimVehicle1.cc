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
// $Id$
//
// Code developed by: S.Guatelli, guatelli@ge.infn.it
//
#include "RemSimVGeometryComponent.hh"
#include "RemSimMaterial.hh"

#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "RemSimVehicle1.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4UserLimits.hh"

RemSimVehicle1::RemSimVehicle1():
  layervacuumPhys(0), layer1Phys(0), layer2Phys(0), layer3Phys(0),
  layer4Phys(0), layer5Phys(0), layer6Phys(0), layer7Phys(0), layer8Phys(0),
  layer9Phys(0), layer10Phys(0), layer11Phys(0), layer12Phys(0),
  layer13Phys(0), layer14Phys(0), layer15Phys(0), layerPhys(0),
  layer16Phys(0), layer17Phys(0), layer18Phys(0), layer19Phys(0),
  layer20Phys(0), layer21Phys(0), layer22Phys(0), layer23Phys(0), 
  layer24Phys(0), layer25Phys(0), layer26Phys(0), layer27Phys(0)
{
 pMaterial = new RemSimMaterial();
}
RemSimVehicle1::~RemSimVehicle1()
{
  delete pMaterial;
}
void RemSimVehicle1::ConstructComponent(G4VPhysicalVolume* motherVolume)
{
  
  G4double sizeX = 5.*m;
  G4double sizeY = 5.*m;
  
  G4Material* betacloth = pMaterial -> GetMaterial("betacloth");  

  //layer of betacloth
  G4double translation = - 5.*m;

  G4double thick = 0.03 *cm;  

  G4Box* layer1 = new G4Box("layer1",sizeX/2.,sizeY/2.,thick/2.);

  G4LogicalVolume* layer1Log = new G4LogicalVolume(layer1,
						   betacloth,
						   "layer1Log",
						   0,0,0);
  
  layer1Phys = new G4PVPlacement(0,
                                 G4ThreeVector(0.,0.,translation+thick/2.),
                                 "layer1Phys", 
				 layer1Log,
				 motherVolume,
				 false,
				  0);

  translation = translation + thick;

  thick = 0.01*cm;

  // Mylar layer
  G4Material* mylar = pMaterial -> GetMaterial("mylar"); 
    
  G4Box* layer2 = new G4Box("layer2",sizeX/2.,sizeY/2.,thick/2.);

  G4LogicalVolume* layer2Log = new G4LogicalVolume(layer2,
                                                   mylar,
                                                   "layer2Log",
                                                    0,0,0);
  
  layer2Phys = new G4PVPlacement(0,
                                 G4ThreeVector(0.,0.,thick/2.+ translation),
                                 "layer2Phys", 
                                 layer2Log,
                                 motherVolume,
				 false,0);
  // Nextel layer
  G4Material* nextel = pMaterial -> GetMaterial("Nextel312AF62"); 
  translation = translation + thick;    
  thick = 0.1*cm;
  
  G4Box* layer3 = new G4Box("layer3",sizeX/2.,sizeY/2.,thick/2.);

  G4LogicalVolume* layer3Log = new G4LogicalVolume(layer3,
                                                   nextel,
                                                   "layer3Log",
                                                    0,0,0);
  
  layer3Phys = new G4PVPlacement(0,
                               G4ThreeVector(0.,0.,thick/2.+ translation),
                                 "layer3Phys", 
                                 layer3Log,
                                 motherVolume,
				 false,0);
  
  translation = translation + thick;
  thick = 10.06*cm;
  G4Material* vacuum = pMaterial -> GetMaterial("Galactic") ;
  G4Box* layer4 = new G4Box("layer4",sizeX/2.,sizeY/2.,thick/2.);

  G4LogicalVolume* layer4Log = new G4LogicalVolume(layer4,
                                                   vacuum,
                                                   "layer4Log",
                                                    0,0,0);
  layer4Phys = new G4PVPlacement(0,
                                 G4ThreeVector(0.,0.,thick/2.+ translation),
                                 "layer4Phys", 
                                 layer4Log,
                                 motherVolume,
				 false,0);
  // Nextel layer
  thick = 0.1*cm;
  translation = translation + 10.06*cm;
  layer5Phys = new G4PVPlacement(0,
                                 G4ThreeVector(0.,0.,thick/2.+ translation),
                                 "layer5Phys", 
                                 layer3Log,
                                 motherVolume,
				 false,0);
  //Vacuum layer 
  translation = translation + thick;
  thick = 10.06 * cm;
  
  layer6Phys = new G4PVPlacement(0,
                                 G4ThreeVector(0.,0.,thick/2.+ translation),
                                 "layer6Phys", 
                                 layer4Log,
                                 motherVolume,
				 false,0);
 // Nextel layer 
  translation = translation + thick;
  thick = 0.2*cm;
  
  layer7Phys = new G4PVPlacement(0,
                                 G4ThreeVector(0.,0.,thick/2.+ translation),
                                 "layer7Phys", 
                                 layer3Log,
                                 motherVolume,
				 false,0);
  // mylar layer
  translation = translation + thick;
  thick = 0.20 *cm;

  G4Box* layer8 = new G4Box("layer8",sizeX/2.,sizeY/2.,thick/2.);

  G4LogicalVolume* layer8Log = new G4LogicalVolume(layer8,
                                                   mylar,
                                                   "layer8Log",
                                                    0,0,0);
  
  layer8Phys = new G4PVPlacement(0,
                                G4ThreeVector(0.,0.,thick/2.+translation),
                                "layer8Phys", 
                                layer8Log,
                                motherVolume,
				false,0);

  //Vacuum layer
  translation = translation + thick;
  thick = 10.06 * cm;

  layer9Phys = new G4PVPlacement(0,
                                 G4ThreeVector(0.,0.,thick/2.+translation),
                                 "layer6Phys", 
                                 layer4Log,
                                 motherVolume,
				 false,0);
  // mylar layer
  translation = translation + thick;
  thick = 0.2 *cm;

  G4Box* layer10 = new G4Box("layer10",sizeX/2.,sizeY/2.,thick/2.);

  G4LogicalVolume* layer10Log = new G4LogicalVolume(layer10,
                                                   mylar,
                                                   "layer10Log",
                                                    0,0,0);
  
  layer10Phys = new G4PVPlacement(0,
                                  G4ThreeVector(0.,0.,thick/2.+translation),
                                  "layer10Phys", 
                                  layer10Log,
                                  motherVolume,
				  false,0);
  //Redundant Bladder
  //Polyethylene
  G4Material* polyethilene = pMaterial -> GetMaterial("polyethylene"); 
  translation = translation + thick;
  thick = 0.0028*cm;

  G4Box* layer11 = new G4Box("layer11", sizeX/2.,sizeY/2.,thick/2.);

  G4LogicalVolume* layer11Log = new G4LogicalVolume(layer11,
                                                    polyethilene,
                                                   "layer11Log",
                                                    0,0,0);
  
  layer11Phys = new G4PVPlacement(0,
                                  G4ThreeVector(0.,0.,thick/2.+translation),
                                  "layer11Phys", 
                                  layer11Log,
                                  motherVolume,
  				  false,0);
  //ployacrylate
  translation = translation + thick;
  thick = 0.0028*cm;
 
  G4Material* polyacrylate = pMaterial ->GetMaterial("polyacrylate"); 
  G4LogicalVolume* layer12Log = new G4LogicalVolume(layer11,
                                                    polyacrylate,
                                                   "layer12Log",
                                                    0,0,0);
  
  layer12Phys = new G4PVPlacement(0,
                                  G4ThreeVector(0.,0.,thick/2.+translation),
                                  "layer12Phys", 
                                  layer12Log,
                                  motherVolume,
  				  false,0);

   //evoh layer
  translation = translation + thick;
  thick = 0.0028*cm;
  G4Material* evoh = pMaterial ->GetMaterial("evoh"); 
  G4LogicalVolume* layer13Log = new G4LogicalVolume(layer11,
                                                    evoh,
                                                   "layer13Log",
                                                    0,0,0);
  
  layer13Phys = new G4PVPlacement(0,
                                  G4ThreeVector(0.,0.,thick/2.+translation),
                                  "layer13Phys", 
                                  layer13Log,
                                  motherVolume,
  				  false,0);
  //polyacrylate
  translation = translation + thick;
  thick = 0.0028*cm;
  
  layer14Phys = new G4PVPlacement(0,
                                  G4ThreeVector(0.,0.,thick/2.+translation),
                                  "layer14Phys", 
                                  layer12Log,
                                  motherVolume,
  				  false,0);
  //polyethilene
  translation = translation + thick;
  thick = 0.0028*cm;
  
  layer15Phys = new G4PVPlacement(0,
                                  G4ThreeVector(0.,0.,thick/2.+translation),
                                  "layer15Phys", 
                                  layer11Log,
                                  motherVolume,
  				  false,0);
  //kevlar layer
  translation = translation + thick;
  thick = 0.026*cm;
  
  G4Material* kevlarAir = pMaterial -> GetMaterial("kevlarAir"); 
  
  G4Box* layer = new G4Box("layer",sizeX/2.,sizeY/2.,thick/2.);

  G4LogicalVolume* layerLog = new G4LogicalVolume(layer,
                                                    kevlarAir,
                                                   "layerLog",
                                                    0,0,0);
  
  layerPhys = new G4PVPlacement(0,
                                G4ThreeVector(0.,0.,thick/2.+translation),
                                "layerPhys", 
                                layerLog,
                                motherVolume,
  			        false,0);
  //polyethilene layer
  translation = translation + thick;
  thick = 0.0028*cm;
  
  layer16Phys = new G4PVPlacement(0,
                                  G4ThreeVector(0.,0.,thick/2.+translation),
                                                    "layer16Phys", 
                                                     layer11Log,
                                                     motherVolume,
  						     false,0);
  //polyacrylate layer
  translation = translation + thick;
  thick = 0.0028*cm;
  
  layer17Phys = new G4PVPlacement(0,
                                  G4ThreeVector(0.,0.,thick/2.+translation),
                                                    "layer17Phys", 
                                                     layer12Log,
                                                     motherVolume,
  						     false,0);
  //evoh
  translation = translation + thick;
  thick = 0.0028*cm;
  
  layer18Phys = new G4PVPlacement(0,
                                  G4ThreeVector(0.,0.,thick/2.+translation),
                                                    "layer18Phys", 
                                                     layer13Log,
                                                     motherVolume,
  						     false,0);
 //polyacrylate layer
  translation = translation + thick;
  thick = 0.0028*cm;
  layer19Phys = new G4PVPlacement(0,
                                  G4ThreeVector(0.,0.,thick/2.+translation),
                                                    "layer19Phys", 
                                                     layer12Log,
                                                     motherVolume,
  						     false,0);
 //polyethilene layer
  translation = translation + thick;
  thick = 0.0028*cm;

  layer20Phys = new G4PVPlacement(0,
                                  G4ThreeVector(0.,0.,thick/2.+translation),
                                                    "layer20Phys", 
                                                     layer11Log,
                                                     motherVolume,
  						     false,0);
  //kevlar layer
  translation = translation + thick;
  thick = 0.026*cm;
    
  layer21Phys = new G4PVPlacement(0,
                                  G4ThreeVector(0.,0.,thick/2.+translation),
                                                    "layer21Phys", 
                                                     layerLog,
                                                     motherVolume,
                                                     false,0);
 //polyethilene layer
  translation = translation + thick;
  thick = 0.0028*cm;

  layer22Phys = new G4PVPlacement(0,
                                  G4ThreeVector(0.,0.,thick/2.+translation),
                                                    "layer22Phys", 
                                                     layer11Log,
                                                     motherVolume,
  						     false,0);
 //polyacrylate layer
  translation = translation + thick;
  thick = 0.0028*cm;
  
  layer23Phys = new G4PVPlacement(0,
                                  G4ThreeVector(0.,0.,thick/2.+translation),
                                                    "layer23Phys", 
                                                     layer12Log,
                                                     motherVolume,
  						     false,0); 						  
 //evoh layer
  translation = translation + thick;
  thick = 0.0028*cm;
  
  layer24Phys = new G4PVPlacement(0,
                                  G4ThreeVector(0.,0.,thick/2.+translation),
                                                    "layer24Phys", 
                                                     layer13Log,
                                                     motherVolume,
  						     false,0);
 //polyacrylate layer
  translation = translation + thick;
  thick = 0.0028*cm;
  
  layer25Phys = new G4PVPlacement(0,
                                  G4ThreeVector(0.,0.,thick/2.+translation),
                                                    "layer25Phys", 
                                                     layer12Log,
                                                     motherVolume,
  						     false,0);
 	
//polyethilene layer
  translation = translation + thick;
  thick = 0.0028*cm;
  
  layer26Phys = new G4PVPlacement(0,
                                  G4ThreeVector(0.,0.,thick/2.+translation),
                                                    "layer26Phys", 
                                                     layer11Log,
                                                     motherVolume,
  						     false,0);
						     
//nomex layer
  translation = translation + thick;  
  thick = 0.03*cm;
 
  G4Material* nomex= pMaterial ->GetMaterial("nomexAir"); 
    
  G4Box* layer27 = new G4Box("layer27",sizeX/2.,sizeY/2.,thick/2.);

  G4LogicalVolume* layer27Log = new G4LogicalVolume(layer27,
                                                    nomex,
                                                   "layer27Log",
                                                    0,0,0);
  
  layer27Phys = new G4PVPlacement(0,
                                  G4ThreeVector(0.,0.,thick/2.+translation),
                                                    "layer27Phys", 
                                                     layer27Log,
                                                     motherVolume,
						    false,0);
  translation = translation + thick;  

  G4Colour  magenta (1.0, 0.0, 1.0) ; 
  G4Colour  green    (0.0,1.0,0.0);
 
  G4VisAttributes* layer1VisAtt = new G4VisAttributes(green);
  layer1VisAtt -> SetVisibility(true);
  layer1VisAtt -> SetForceSolid(true);
  layer1Log -> SetVisAttributes(layer1VisAtt); 
  layer3Log -> SetVisAttributes(layer1VisAtt);   
  layer8Log -> SetVisAttributes(layer1VisAtt);  
  layer11Log -> SetVisAttributes(layer1VisAtt); 
  layer13Log -> SetVisAttributes(layer1VisAtt);
 
  G4VisAttributes* layer2VisAtt = new G4VisAttributes(magenta);
  layer2VisAtt -> SetVisibility(true);
  layer2VisAtt -> SetForceSolid(true);
  layer2Log -> SetVisAttributes(layer2VisAtt);  
  layer4Log -> SetVisAttributes(layer2VisAtt);  
  layer10Log -> SetVisAttributes(layer2VisAtt);  
  layer12Log -> SetVisAttributes(layer2VisAtt);
  layer27Log -> SetVisAttributes(layer2VisAtt);
}

void RemSimVehicle1::DestroyComponent()
{
}
