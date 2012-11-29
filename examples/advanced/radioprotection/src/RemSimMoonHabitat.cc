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
//    *    RemSimMoonHabitat.cc            *
//    *                                    *          
//    **************************************
//
// $Id$
//
// Author:Susanna Guatelli, guatelli@ge.infn.it 
//
#include "RemSimMoonHabitat.hh"
#include "RemSimMaterial.hh"

#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4UserLimits.hh"
#include "RemSimDecorator.hh"
#include "RemSimAstronautDecorator.hh"
#include "G4SubtractionSolid.hh"

RemSimMoonHabitat::RemSimMoonHabitat()
{
  pMaterial = new RemSimMaterial(); 
  moonPhys = 0;
}
RemSimMoonHabitat::~RemSimMoonHabitat()
{
  delete pMaterial;
}
void RemSimMoonHabitat::ConstructComponent(G4VPhysicalVolume* motherVolume)
{
  pMaterial -> DefineMaterials();

  G4Material* moonMaterial = pMaterial->GetMaterial("moon");

  // Moon Surface definition
  G4double sizeX = 30.*m;
  G4double sizeY = 30.*m;
  G4double sizeZ = 15.*m;

  G4Box* moon = new G4Box("moon",sizeX,sizeY,sizeZ);
     
  G4double thick = 4.5 *m;
  G4Box* shelter = new G4Box("shelter",5.*m,5.*m,thick/2.);

 G4ThreeVector  translation(0, 0, -sizeZ + 0.5 * m + thick/2.);
 G4SubtractionSolid* solid  = new G4SubtractionSolid ("moon-shelter", moon,
  shelter, 0, translation);
 
 G4LogicalVolume* moonLogical = new G4LogicalVolume(solid, 
  						    moonMaterial,
  					            "moon",0,0,0);

 moonPhys = new G4PVPlacement(0,
  			      G4ThreeVector(0.,0.,sizeZ),                                                     "moon",moonLogical,
  			      motherVolume,
  			      false,0);

  G4Colour  red      (1.0,0.0,0.0);
     
  G4VisAttributes* moonVisAtt = new G4VisAttributes(red);
  moonVisAtt -> SetVisibility(true);
  moonVisAtt -> SetForceWireframe(true);
  moonLogical -> SetVisAttributes(moonVisAtt); 
}

void RemSimMoonHabitat::DestroyComponent()
{}
 
