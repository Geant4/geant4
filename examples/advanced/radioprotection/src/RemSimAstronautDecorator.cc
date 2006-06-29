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
//    **********************************
//    *                                *
//    *    RemSimAstronautDecorator.cc *
//    *                                *
//    **********************************
//
// $Id: RemSimAstronautDecorator.cc,v 1.7 2006-06-29 16:23:33 gunter Exp $
//
// Author:Susanna Guatelli, guatelli@ge.infn.it 

#include "RemSimVGeometryComponent.hh"
#include "RemSimMaterial.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "RemSimAstronautDecorator.hh"
#include "RemSimDecorator.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "RemSimSensitiveDetector.hh"
#include "RemSimROGeometry.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4UserLimits.hh"
#include "RemSimMoonHabitat.hh"

RemSimAstronautDecorator::RemSimAstronautDecorator(RemSimVGeometryComponent* comp, G4bool moon)
  : RemSimDecorator(comp), phantom(0), phantomLog(0), phantomPhys(0)
{ 
 flag = moon;
 sensitiveDetector = 0;
}
RemSimAstronautDecorator::~RemSimAstronautDecorator()
{
}
void RemSimAstronautDecorator::ConstructComponent(G4VPhysicalVolume* motherVolume)
{
  RemSimDecorator::ConstructComponent(motherVolume);
  ConstructAstronaut(motherVolume);
}

void RemSimAstronautDecorator::DestroyComponent()
{

}

void RemSimAstronautDecorator::ConstructAstronaut(G4VPhysicalVolume* motherVolume)
{
  // Astronaut definition: Box of water
  // The Astronaut is the sensitive detector
 
  RemSimMaterial*  pMaterial = new RemSimMaterial();

  G4Material* water = pMaterial -> GetMaterial("Water");
  delete pMaterial;

  // Phantom sizes
  G4double phantomX = 3. *m;
  G4double phantomY = 3. *m;
  G4double phantomZ = 30. *cm;
 

  phantom = new G4Box("phantom",phantomX/2.,phantomY/2.,phantomZ/2.);

  phantomLog = new G4LogicalVolume(phantom,
                                   water,
                                   "phantom",
                                   0,0,0);
  G4double translation1 = 0.; 
 if (flag == true) 
    {
      G4double thickShelter = 4.5 *m;
      translation1 = 0.5*m + thickShelter/2.;
   }
 
  phantomPhys = new G4PVPlacement(0,
                                  G4ThreeVector(0.,0.,translation1),
                                  "phantom",phantomLog, 
                                   motherVolume,false,0);
  
  // Visualisation attributes
  G4Colour  lblue   (0.0, 0.0,.75); 
  G4VisAttributes* phantomVisAtt = new G4VisAttributes(lblue);
  phantomVisAtt -> SetVisibility(true);
  phantomVisAtt -> SetForceSolid(true);
  phantomLog -> SetVisAttributes(phantomVisAtt); 
 
 //Sensitive Detector  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4String sensitiveDetectorName = "AstronautSD";
  
  
  sensitiveDetector = new  RemSimSensitiveDetector(sensitiveDetectorName);
  G4int VoxelNbAlongZ = 30;
 
  G4double  translation = 0;
    
  if (flag == true) 
    {
      G4double thickShelter = 4.5 *m;
      translation = 0.5*m + thickShelter/2.;
   }
 
  RemSimROGeometry* ROGeometry = new RemSimROGeometry(phantomX,
                                                      phantomY,
                                                      phantomZ,
						      VoxelNbAlongZ, 
                                                      translation);
  ROGeometry -> BuildROGeometry();
  sensitiveDetector -> SetROgeometry(ROGeometry);
  SDman->AddNewDetector(sensitiveDetector);
  phantomLog -> SetSensitiveDetector(sensitiveDetector);

  // Set max step allowd to particles in the phantom
  phantomLog -> SetUserLimits(new G4UserLimits(0.1*cm));
}

void RemSimAstronautDecorator::ChangeThickness(G4double)
{
  G4cout << "It is not possible to change the sizes of the Astronaut"
         << G4endl;
}

void RemSimAstronautDecorator::PrintDetectorParameters()
{
  G4cout << "-----------------------------------------------------------------------"
         << G4endl
         << "the astronaut is a box whose thickness is: " << G4endl
         << ((phantom -> GetZHalfLength())*2.)/cm
         << " cm along the Z axis"
         << G4endl
         << "material of the astronaut: "
         << phantomLog -> GetMaterial() -> GetName() <<G4endl
         << G4endl;
}

