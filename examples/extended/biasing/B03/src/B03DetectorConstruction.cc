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
// $Id: B03DetectorConstruction.cc,v 1.8 2002-08-29 15:37:25 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "g4std/strstream"
#include "globals.hh"

#include "B03DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "PhysicalConstants.h"

// for importance biasing
#include "G4IStore.hh"

B03DetectorConstruction::B03DetectorConstruction()
{;}

B03DetectorConstruction::~B03DetectorConstruction()
{;}

G4IStore* B03DetectorConstruction::GetIStore()
{
  if (!fIStore) G4Exception("B03DetectorConstruction::fIStore empty!");
  return fIStore;
}

G4VPhysicalVolume* B03DetectorConstruction::Construct()
{
  char line[255];
  
  G4double pos_x;
  G4double pos_y;
  G4double pos_z; 

  G4double density, pressure, temperature;
  G4double A;
  G4int Z;

  G4String name, symbol;
  G4double z;
  G4double fractionmass;

  A = 1.01*g/mole;
  G4Element* elH  = new G4Element(name="Hydrogen",symbol="H" , Z= 1, A);

  A = 12.01*g/mole;
  G4Element* elC  = new G4Element(name="Carbon"  ,symbol="C" , Z = 6, A);

  A = 16.00*g/mole;
  G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , Z= 8, A);

  A = 22.99*g/mole; 
  G4Element* elNa  = new G4Element(name="Natrium"  ,symbol="Na" , Z=11 , A);

  A = 200.59*g/mole; 
  G4Element* elHg  = new G4Element(name="Hg"  ,symbol="Hg" , Z=80, A);

  A = 26.98*g/mole; 
  G4Element* elAl  = new G4Element(name="Aluminium"  ,symbol="Al" , Z=13, A);

  A = 28.09*g/mole;
  G4Element* elSi  = new G4Element(name="Silicon", symbol="Si", Z=14, A);


  A = 39.1*g/mole; 
  G4Element* elK  = new G4Element(name="K"  ,symbol="K" , Z=19 , A);

  A = 69.72*g/mole; 
  G4Element* elCa  = new G4Element(name="Calzium"  ,symbol="Ca" , Z=31 , A);

  A = 55.85*g/mole;
  G4Element* elFe = new G4Element(name="Iron"    ,symbol="Fe", Z=26, A);

  density     = universe_mean_density;            //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  G4Material *Galactic = 
    new G4Material(name="Galactic", z=1., A=1.01*g/mole, density,
                   kStateGas,temperature,pressure);

  density = 2.03*g/cm3;
  G4Material* Concrete = new G4Material("Concrete", density, 10);
  Concrete->AddElement(elH , fractionmass= 0.01);
  Concrete->AddElement(elO , fractionmass= 0.529);
  Concrete->AddElement(elNa , fractionmass= 0.016);
  Concrete->AddElement(elHg , fractionmass= 0.002);
  Concrete->AddElement(elAl , fractionmass= 0.034);
  Concrete->AddElement(elSi , fractionmass= 0.337);
  Concrete->AddElement(elK , fractionmass= 0.013);
  Concrete->AddElement(elCa , fractionmass= 0.044);
  Concrete->AddElement(elFe , fractionmass= 0.014);
  Concrete->AddElement(elC , fractionmass= 0.001);

  density = 0.0203*g/cm3;
  G4Material* LightConcrete = new G4Material("LightConcrete", density, 10);
  LightConcrete->AddElement(elH , fractionmass= 0.01);
  LightConcrete->AddElement(elO , fractionmass= 0.529);
  LightConcrete->AddElement(elNa , fractionmass= 0.016);
  LightConcrete->AddElement(elHg , fractionmass= 0.002);
  LightConcrete->AddElement(elAl , fractionmass= 0.034);
  LightConcrete->AddElement(elSi , fractionmass= 0.337);
  LightConcrete->AddElement(elK , fractionmass= 0.013);
  LightConcrete->AddElement(elCa , fractionmass= 0.044);
  LightConcrete->AddElement(elFe , fractionmass= 0.014);
  LightConcrete->AddElement(elC , fractionmass= 0.001);
   
  G4Material *WorldMaterial = Galactic; 

  /////////////////////////////
  // world cylinder volume
  ////////////////////////////

  // world solid

  G4double innerRadiusCylinder = 0*cm;
  G4double outerRadiusCylinder = 100*cm;
  G4double hightCylinder       = 15.001*cm;
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

  G4std::vector< G4VPhysicalVolume * > physvolumes;
  physvolumes.push_back(worldCylinder_phys);

  ///////////////////////////////////////////////
  // shield cylinder for (cells 2-4)
  ////////////////////////////////////////////////

  G4double innerRadiusShield = 0*cm;
  G4double outerRadiusShield = 100*cm;
  G4double hightShield       = 5*cm;
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

  //        physical shields for cell 2 to 4
  for(G4int i=0; i<3; i++) {
    if (i+2<10)    sprintf(line,"cell: 0%d, shield",i+2);
    else sprintf(line,"cell: %d, shield",i+2);
    G4String name(line);
    pos_x = 0*cm;
    pos_y = 0*cm;
    pos_z = -10*cm+(10*cm)*i;
    physvolumes.
      push_back(new G4PVPlacement(0, G4ThreeVector(pos_x, pos_y, pos_z),
				  aShield_log, name, worldCylinder_log, 
				  false, 0));
  }

  fIStore = new G4IStore(*worldCylinder_phys);
  // for the world volume repnum is -1 !!!!!!!!!!!!!!!!!!!!!!!
  G4int n = 0;
  for (G4std::vector<G4VPhysicalVolume *>::iterator it = physvolumes.begin();
       it != physvolumes.end(); it++) {
    G4double i = pow(2., n++);
    G4cout << "Going to assign importance: " << i << ", to volume" 
	   << (*it)->GetName() << G4endl;
    if (*it == worldCylinder_phys) {
      // repnum -1 
      fIStore->AddImportanceGeometryCell(i, **it, -1); 
    }
    fIStore->AddImportanceGeometryCell(i, **it);
  }
  
  return worldCylinder_phys;
}
