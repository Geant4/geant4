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
//    **************************************
//    *                                    *
//    *    RemSimMoonHabitat.cc            *
//    *                                    *          
//    **************************************
//
// $Id: RemSimMoonHabitat.cc,v 1.4 2004/05/22 12:57:06 guatelli Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
// Author:Susanna Guatelli, guatelli@ge.infn.it 
//
#include "RemSimVGeometryComponent.hh"
#include "RemSimMaterial.hh"
#include "G4Material.hh"
#include "RemSimMoonHabitat.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4UserLimits.hh"
#include "RemSimDecorator.hh"
#include "RemSimAstronautDecorator.hh"

RemSimMoonHabitat::RemSimMoonHabitat()
{
  pMaterial = new RemSimMaterial(); 
  shelterPhys = 0;
}
RemSimMoonHabitat::~RemSimMoonHabitat()
{
  delete pMaterial;
}
void RemSimMoonHabitat::ConstructComponent(G4VPhysicalVolume* motherVolume)
{
  pMaterial -> DefineMaterials();

  G4Material* moonMaterial = pMaterial->GetMaterial("moon");
  G4Material* air = pMaterial -> GetMaterial("Air");

  // Moon Surface definition
  G4double sizeX = 30.*m;
  G4double sizeY = 30.*m;
  G4double sizeZ = 15.*m;

  G4Box* moon = new G4Box("moon",sizeX,sizeY,sizeZ);
     
  G4LogicalVolume* moonLogical = new G4LogicalVolume(moon, 
						     moonMaterial,
						     "moon",0,0,0);

  G4VPhysicalVolume* moonPhys = new G4PVPlacement(0,
						  G4ThreeVector(0.,0.,sizeZ),                                                     "moon",moonLogical,
						  motherVolume,
						  false,0);

  G4double thick = 4.5 *m;
  G4Box* shelter = new G4Box("shelter",5.*m,5.*m,thick/2.);
  G4LogicalVolume* shelterLog = new G4LogicalVolume(shelter,
                                                    air,
                                                    "shelter",
                                                    0,0,0);
 
  G4double translation = - sizeZ + thick/2. + 0.5 *m;
 
  shelterPhys = new G4PVPlacement(0,
                                  G4ThreeVector(0.,0.,translation),
                                  "shelter",shelterLog, 
                                  moonPhys,false,0);
  
  G4Colour  red      (1.0,0.0,0.0);
  G4Colour  blue   (0.,0.,1.);

  G4VisAttributes* moonVisAtt = new G4VisAttributes(red);
  moonVisAtt -> SetVisibility(true);
  moonVisAtt -> SetForceWireframe(true);
  moonLogical -> SetVisAttributes(moonVisAtt); 
   
  G4VisAttributes* shelterVisAtt = new G4VisAttributes(blue);
  shelterVisAtt -> SetVisibility(true);
  shelterVisAtt -> SetForceWireframe(true);
  shelterLog -> SetVisAttributes(shelterVisAtt);  
}

void RemSimMoonHabitat::DestroyComponent()
{}
 
G4VPhysicalVolume* RemSimMoonHabitat::GetShelter()
{
  return shelterPhys;
}
