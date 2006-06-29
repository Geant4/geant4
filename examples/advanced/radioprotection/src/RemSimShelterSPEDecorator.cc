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
// $Id: RemSimShelterSPEDecorator.cc,v 1.3 2006-06-29 16:24:21 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Code developed by: S.Guatelli, guatelli@ge.infn.it
//
#include "RemSimVGeometryComponent.hh"
#include "RemSimMaterial.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "RemSimShelterSPEDecorator.hh"
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

RemSimShelterSPEDecorator::RemSimShelterSPEDecorator(RemSimVGeometryComponent* comp)
  : RemSimDecorator(comp)
{
 shelterSPEX = 5.*m;
 shelterSPEY = 5.*m;
 shelterSPEZ = 75.*cm;  
 translation = - 2.* m;
 pMaterial = new RemSimMaterial();

 shelterSPEVisAtt = 0;
}
RemSimShelterSPEDecorator::~RemSimShelterSPEDecorator()
{
  delete pMaterial;
}
void RemSimShelterSPEDecorator::ConstructComponent(G4VPhysicalVolume* motherVolume)
{
  pMaterial -> DefineMaterials();
  RemSimDecorator::ConstructComponent(motherVolume);
  ConstructShelterSPE(motherVolume);
}

void RemSimShelterSPEDecorator::DestroyComponent()
{
 delete shelterSPEVisAtt;
 shelterSPEVisAtt = 0;
 
 delete shelterSPEPhys; 
 shelterSPEPhys = 0;

 delete shelterSPELog;
 shelterSPELog = 0;

 delete shelterSPE;
 shelterSPE = 0;
}
void RemSimShelterSPEDecorator::ConstructShelterSPE(G4VPhysicalVolume* motherVolume)
{
  // Geometry definition
  pMaterial -> DefineMaterials();

  G4Material* water = pMaterial -> GetMaterial("Water");
  
  shelterSPE = new G4Box("shelterSPE",shelterSPEX/2.,shelterSPEY/2.,shelterSPEZ/2.);

  shelterSPELog = new G4LogicalVolume(shelterSPE, water,
                                     "shelterSPELog",0,0,0);
  
  shelterSPEPhys = new G4PVPlacement(0,
             G4ThreeVector(0.,0.,translation + shelterSPEZ/2.),
            "shelterSPEPhys", shelterSPELog, motherVolume,false,0); 

  //Visualisation attributes
  G4Colour  red      (1.0,0.0,0.0);
  shelterSPEVisAtt = new G4VisAttributes(red);
  shelterSPEVisAtt -> SetVisibility(true);
  shelterSPEVisAtt -> SetForceSolid(true);
  shelterSPELog -> SetVisAttributes(shelterSPEVisAtt);  
  PrintDetectorParameters();
}

void RemSimShelterSPEDecorator::ChangeThickness(G4double)
{
  ;
}

void RemSimShelterSPEDecorator::PrintDetectorParameters()
{
  G4cout << "-----------------------------------------------------------------------"
         << G4endl
         << "the shelterSPE is a box whose thickness is: " << G4endl
         << ((shelterSPE -> GetZHalfLength())*2.)/cm
         << " cm along the Z axis"
         << G4endl
         << "material of the shelterSPE: "
         << shelterSPELog -> GetMaterial() -> GetName() <<G4endl
         << G4endl;
}
