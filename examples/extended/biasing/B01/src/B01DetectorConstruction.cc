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
// $Id: B01DetectorConstruction.cc,v 1.8 2003-02-19 08:02:36 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "g4std/strstream"
#include "globals.hh"

#include "B01DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

// for importance biasing
#include "G4IStore.hh"

B01DetectorConstruction::B01DetectorConstruction()
 : fIStore(0)
{;}

B01DetectorConstruction::~B01DetectorConstruction()
{;}

G4IStore* B01DetectorConstruction::GetIStore()
{
  if (!fIStore) G4Exception("B01DetectorConstruction::fIStore empty!");
  return fIStore;
}

G4VPhysicalVolume* B01DetectorConstruction::Construct()
{
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
   














  /////////////////////////////
  // world cylinder volume
  ////////////////////////////

  // world solid

  G4double innerRadiusCylinder = 0*cm;
  G4double outerRadiusCylinder = 101*cm; // dont't have scoring
                   // cells coinside eith world volume boundary
  G4double hightCylinder       = 105*cm;
  G4double startAngleCylinder  = 0*deg;
  G4double spanningAngleCylinder    = 360*deg;

  G4Tubs *worldCylinder = new G4Tubs("worldCylinder",
                                     innerRadiusCylinder,
                                     outerRadiusCylinder,
                                     hightCylinder,
                                     startAngleCylinder,
                                     spanningAngleCylinder);

  // logical world

  G4LogicalVolume *worldCylinder_log = 
    new G4LogicalVolume(worldCylinder, Galactic, "worldCylinder_log");

  name = "shieldWorld";
  G4VPhysicalVolume *pWorldVolume = new 
    G4PVPlacement(0, G4ThreeVector(0,0,0), worldCylinder_log,
		  name, 0, false, 0);


  G4std::vector< G4VPhysicalVolume * > physvolumes;
  physvolumes.push_back(pWorldVolume);



  // creating 18 slobs of 10 cm thick concrete

  G4double innerRadiusShield = 0*cm;
  G4double outerRadiusShield = 100*cm;
  G4double hightShield       = 5*cm;
  G4double startAngleShield  = 0*deg;
  G4double spanningAngleShield    = 360*deg;

  G4Tubs *aShield = new G4Tubs("aShield",
                               innerRadiusShield,
                               outerRadiusShield,
                               hightShield,
                               startAngleShield,
                               spanningAngleShield);
  
  // logical shield

  G4LogicalVolume *aShield_log = 
    new G4LogicalVolume(aShield, Concrete, "aShield_log");

  G4VisAttributes* pShieldVis = new 
    G4VisAttributes(G4Colour(0.0,0.0,1.0));
  pShieldVis->SetForceSolid(true);
  aShield_log->SetVisAttributes(pShieldVis);

  // physical shields

  G4int i;
  G4double startz = -85*cm; 
  for (i=1; i<=18; i++) {
    name = GetCellName(i);

    G4double pos_x = 0*cm;
    G4double pos_y = 0*cm;
    G4double pos_z = startz + (i-1) * (2*hightShield);
    G4VPhysicalVolume *pvol = 
      new G4PVPlacement(0, 
			G4ThreeVector(pos_x, pos_y, pos_z),
			aShield_log, 
			name, 
			worldCylinder_log, 
			false, 
			0);
    physvolumes.push_back(pvol);
  }

  // filling the rest of the world volumr behind the concrete with
  // another slob which should get the same importance value as the 
  // last slob
  innerRadiusShield = 0*cm;
  outerRadiusShield = 100*cm;
  hightShield       = 7.5*cm;
  startAngleShield  = 0*deg;
  spanningAngleShield    = 360*deg;

  G4Tubs *aRest = new G4Tubs("Rest",
			     innerRadiusShield,
			     outerRadiusShield,
			     hightShield,
			     startAngleShield,
			     spanningAngleShield);
  
  G4LogicalVolume *aRest_log = 
    new G4LogicalVolume(aRest, Galactic, "aRest_log");
  name = "rest";
    
  pos_x = 0*cm;
  pos_y = 0*cm;
  pos_z = 97.5*cm;
  G4VPhysicalVolume *pvol_rest = 
    new G4PVPlacement(0, 
		      G4ThreeVector(pos_x, pos_y, pos_z),
		      aRest_log, 
		      name, 
		      worldCylinder_log, 
		      false, 
		      0);


  // creating and filling the importance store

  fIStore = new G4IStore(*pWorldVolume);

  // for the world volume repnum is -1 !
  G4int n = 0;
  G4double imp =1;
  for (G4std::vector<G4VPhysicalVolume *>::iterator it =
	 physvolumes.begin();
       it != physvolumes.end(); it++)
  {
    imp = pow(2., n++);
    G4cout << "Going to assign importance: " << imp << ", to volume: " 
	   << (*it)->GetName() << G4endl;
    if (*it == pWorldVolume)
    {
      // repnum -1 
      fIStore->AddImportanceGeometryCell(imp, **it, -1); 
    }
    fIStore->AddImportanceGeometryCell(imp, **it);
  }

  // the remaining part pf the geometry (rest) gets the same
  // importance as the last conrete cell
  fIStore->AddImportanceGeometryCell(imp, *pvol_rest);


  return pWorldVolume;
}



G4String B01DetectorConstruction::GetCellName(G4int i) {
  char st[200];
  G4std::ostrstream os(st,200);
  os << "cell_";
  if (i<10) {
    os << "0";
  }
  os << i 
     << '\0';
  G4String name(st);
  return name;
}
