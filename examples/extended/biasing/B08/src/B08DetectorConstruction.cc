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
// $Id: B08DetectorConstruction.cc,v 1.1 2002/06/04 11:14:52 dressel Exp $
// GEANT4 tag $Name: geant4-04-01 $
//

#include "g4std/strstream"
#include "globals.hh"

#include "B08DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "PhysicalConstants.h"

B08DetectorConstruction::B08DetectorConstruction()
{;}

B08DetectorConstruction::~B08DetectorConstruction()
{;}

G4VPhysicalVolume* B08DetectorConstruction::Construct()
{
  G4double pos_x;
  G4double pos_y;
  G4double pos_z; 

  G4double density, pressure, temperature;
  G4double A;

  G4String name, symbol;
  G4double z;

  G4Element *pElementH = new G4Element("Hydro","H",1.,1.0079*g/mole);
  G4Element *pElementO = new G4Element("Oxy","O",8.,15.999*g/mole);
  G4Element *pElementFe = new G4Element("Ferrum","Fe",26.,55.845*g/mole);
  G4Element *pElementNa = new G4Element("Natri","Na",11.,22.98977*g/mole);
  G4Element *pElementMg = new G4Element("Magn","Mg",12.,24.305*g/mole);
  G4Element *pElementAl = new G4Element("Alu","Al",13.,26.981538*g/mole);
  G4Element *pElementSi = new G4Element("Sili","Si",14.,28.0854*g/mole);
  G4Element *pElementK = new G4Element("Pota","K",19.,39.0983*g/mole);
  G4Element *pElementCa = new G4Element("Calc","Ca",20.,40.078*g/mole);


  G4Material *pMaterialConcrete = new G4Material("Concrete",2.31*g/cm3,9);
  pMaterialConcrete->AddElement(pElementH,/*0.18957226*/0.010853422);
  pMaterialConcrete->AddElement(pElementO, /*0.52999241*/0.48165594);
  pMaterialConcrete->AddElement(pElementNa,/*0.01556568*/0.020327181);
  pMaterialConcrete->AddElement(pElementMg,/*0.0078461149*/0.010832401);
  pMaterialConcrete->AddElement(pElementAl,/*0.039483675*/0.060514396);
  pMaterialConcrete->AddElement(pElementSi,/*0.14047077*/0.22409956);
  pMaterialConcrete->AddElement(pElementK, /*0.0048089091*/0.010680188);
  pMaterialConcrete->AddElement(pElementCa,/*0.054416603*/0.12388306);
  pMaterialConcrete->AddElement(pElementFe,/*0.017843584*/0.056603178);


  G4Material* Concrete = pMaterialConcrete;

  density     = universe_mean_density;            //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  G4Material *Galactic = 
    new G4Material(name="Galactic", z=1., A=1.01*g/mole, density,
                   kStateGas,temperature,pressure);
   
  G4Material *WorldMaterial = Galactic; // tracks entering this are 
                                        // killed immediately

  /////////////////////////////
  // world cylinder volume
  ////////////////////////////

  // world solid

  G4double innerRadiusCylinder = 0*cm;
  G4double outerRadiusCylinder = 100*cm;
  G4double hightCylinder       = 100*cm;
  G4double startAngleCylinder  = 0*deg;
  G4double spanningAngleCylinder    = 360*cm;

  G4Tubs *worldCylinder = new G4Tubs("worldCylinder",
				     innerRadiusCylinder,
				     outerRadiusCylinder,
				     hightCylinder,
				     startAngleCylinder,
				     spanningAngleCylinder);

  // logical world

  G4LogicalVolume *worldCylinder_log = 
    new G4LogicalVolume(worldCylinder, WorldMaterial, "worldCylinder_log");

  G4VisAttributes * WorldVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  WorldVisAtt->SetVisibility(true);
  worldCylinder_log->SetVisAttributes(WorldVisAtt);
  
  // physical world

  name = "worldCylinder_phys";
  G4VPhysicalVolume* worldCylinder_phys =
    new G4PVPlacement(0, G4ThreeVector(0,0,0), worldCylinder_log,
		      name, 0, false, 0);

  ///////////////////////////////////////////////
  // the concreate block
  ////////////////////////////////////////////////

  G4double innerRadiusShield = 0*cm;
  G4double outerRadiusShield = 100*cm;
  G4double hightShield       = 90*cm;
  G4double startAngleShield  = 0*deg;
  G4double spanningAngleShield    = 360*cm;

  G4Tubs *aShield = new G4Tubs("aShield",
			       innerRadiusShield,
			       outerRadiusShield,
			       hightShield,
			       startAngleShield,
			       spanningAngleShield);
  
  // logical shield

  G4LogicalVolume *aShield_log = 
    new G4LogicalVolume(aShield, Concrete, "aShield_log");

  G4VisAttributes * shieldVisAtt
    = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  shieldVisAtt->SetVisibility(true);
  shieldVisAtt->SetForceSolid(false);
  aShield_log->SetVisAttributes(shieldVisAtt);

  // physical shields

  name = "ConcreateBlock";
  pos_x = 0*cm;
  pos_y = 0*cm;
  pos_z = 0*cm;
  new G4PVPlacement(0, G4ThreeVector(pos_x, pos_y, pos_z),
		    aShield_log, name, worldCylinder_log, false, 0);
  

  return worldCylinder_phys;
}
