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
//    ************************************
//    *                                  *
//    *      RemSimRoofDecorator.cc      *
//    *                                  *
//    ************************************
//
// Code developed by: S.Guatelli, guatelli@ge.infn.it
//
//
// $Id: RemSimRoofDecorator.cc,v 1.7 2006-06-29 16:24:15 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "RemSimVGeometryComponent.hh"
#include "RemSimMaterial.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "RemSimRoofDecorator.hh"
#include "RemSimDecorator.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "RemSimSensitiveDetector.hh"
#include "RemSimROGeometry.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4VisAttributes.hh"

RemSimRoofDecorator::RemSimRoofDecorator(RemSimVGeometryComponent* comp)
  : RemSimDecorator(comp)
{
  pMaterial = new RemSimMaterial();
  roofVisAtt = 0;
  roof = 0;
  roofLog = 0;
  roofPhys = 0;
}
RemSimRoofDecorator::~RemSimRoofDecorator()
{
  delete pMaterial;
}
void RemSimRoofDecorator::ConstructComponent(G4VPhysicalVolume* motherVolume)
{
  pMaterial -> DefineMaterials();
  RemSimDecorator::ConstructComponent(motherVolume);
  ConstructRoof(motherVolume);
}

void RemSimRoofDecorator::DestroyComponent()
{
  delete roofVisAtt;
  roofVisAtt = 0;
 
  delete roofPhys; 
  roofPhys = 0;

  delete roofLog;
  roofLog = 0;

  delete roof;
  roof = 0;
}
void RemSimRoofDecorator::ConstructRoof(G4VPhysicalVolume* motherVolume)
{
  // Geometry definition
  pMaterial -> DefineMaterials();

  G4Material* moonMaterial = pMaterial -> GetMaterial("moon");
  
 
  roof = new G4Trd("roof",2.25*m,7.*m,5.*m, 8.*m, 0.5 *m);

  roofLog = new G4LogicalVolume(roof,
				moonMaterial,
				"roof",
				0,0,0);

  roofPhys = new G4PVPlacement(0,
			       G4ThreeVector(0.,0., - 0.5 *m),
			       "roof",roofLog, 
			       motherVolume,false,0);

  //Visualisation attributes
  G4Colour  red      (1.0,0.0,0.0);
  roofVisAtt = new G4VisAttributes(red);
  roofVisAtt -> SetVisibility(true);
  roofVisAtt -> SetForceSolid(true);
  roofLog -> SetVisAttributes(roofVisAtt);  
  PrintDetectorParameters();
}

void RemSimRoofDecorator::ChangeThickness(G4double thick)
{
  PrintDetectorParameters();
  roof -> SetZHalfLength(thick/2.);
  roofPhys -> SetTranslation(G4ThreeVector 
			     (0.,0., - thick/2.));
}

void RemSimRoofDecorator::PrintDetectorParameters()
{
  ;
}
