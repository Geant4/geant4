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
// $Id: RemSimShieldingDecorator.cc,v 1.6 2006-06-29 16:24:23 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Code developed by: S.Guatelli, guatelli@ge.infn.it
//
#include "RemSimVGeometryComponent.hh"
#include "RemSimMaterial.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "RemSimShieldingDecorator.hh"
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
#include "G4VisAttributes.hh"

RemSimShieldingDecorator::RemSimShieldingDecorator(RemSimVGeometryComponent* comp)
  : RemSimDecorator(comp)
{
 shieldingX = 5.*m;
 shieldingY = 5.*m;
 shieldingZ = 10.*cm;  
 translation = -4.4* m;
 pMaterial = new RemSimMaterial();

 shieldingVisAtt = 0;
}
RemSimShieldingDecorator::~RemSimShieldingDecorator()
{
  delete pMaterial;
}
void RemSimShieldingDecorator::ConstructComponent(G4VPhysicalVolume* motherVolume)
{
  pMaterial -> DefineMaterials();
  RemSimDecorator::ConstructComponent(motherVolume);
  ConstructShielding(motherVolume);
}

void RemSimShieldingDecorator::DestroyComponent()
{
 delete shieldingVisAtt;
 shieldingVisAtt = 0;
 
 delete shieldingPhys; 
 shieldingPhys = 0;

 delete shieldingLog;
 shieldingLog = 0;

 delete shielding;
 shielding = 0;
}
void RemSimShieldingDecorator::ConstructShielding(G4VPhysicalVolume* motherVolume)
{
  // Geometry definition
  pMaterial -> DefineMaterials();

  G4Material* water = pMaterial -> GetMaterial("Water");
  
  shielding = new G4Box("shielding",shieldingX/2.,shieldingY/2.,shieldingZ/2.);

  shieldingLog = new G4LogicalVolume(shielding, water,
                                     "shieldingLog",0,0,0);
  
  shieldingPhys = new G4PVPlacement(0,
             G4ThreeVector(0.,0.,translation + shieldingZ/2.),
            "shieldingPhys", shieldingLog, motherVolume,false,0); 

  //Visualisation attributes
  G4Colour  red      (1.0,0.0,0.0);
  shieldingVisAtt = new G4VisAttributes(red);
  shieldingVisAtt -> SetVisibility(true);
  shieldingVisAtt -> SetForceSolid(true);
  shieldingLog -> SetVisAttributes(shieldingVisAtt);  
  PrintDetectorParameters();
}

void RemSimShieldingDecorator::ChangeThickness(G4double thick)
{
  PrintDetectorParameters();
  shielding -> SetZHalfLength(thick/2.);
  shieldingPhys -> SetTranslation(G4ThreeVector 
                            (0.,0.,translation + thick/2.));
}

void RemSimShieldingDecorator::PrintDetectorParameters()
{
  G4cout << "-----------------------------------------------------------------------"
         << G4endl
         << "the shielding is a box whose thickness is: " << G4endl
         << ((shielding -> GetZHalfLength())*2.)/cm
         << " cm along the Z axis"
         << G4endl
         << "material of the shielding: "
         << shieldingLog -> GetMaterial() -> GetName() <<G4endl
         << G4endl;
}
