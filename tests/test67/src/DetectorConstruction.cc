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
// $Id: DetectorConstruction.cc,v 1.1 2009/03/21 18:37:27 vnivanch Exp $
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"


#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4Color.hh"

#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include <iomanip>

#include "G4SDManager.hh"
#include "SensitiveDetector.hh"

#include "DetectorMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
  //Initialization
  casingHeight=70.0*mm;
  geometryType = 2; //this is what is called geometry #2 in the paper

  waterDiameter = 90.0*mm;
  waterHeight = 40.0*mm;
  waterZPosition = casingHeight/2. + 5.0*mm + waterHeight/2.;

  theDeteMessenger = new DetectorMessenger(this);
}
  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete theDeteMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  DefineMaterials();

  //geometry #2 of the MC intercomparison test

  //************************************************************************
  //Mother volume of vacuum
  //************************************************************************
  G4double world_x = 40*cm;
  G4double world_y = 40*cm;
  G4double world_z = 40*cm;
  G4Box* world_box
    = new G4Box("world_box",world_x/2.,world_y/2.,world_z/2.);
  G4LogicalVolume* world_log = new G4LogicalVolume(world_box,
						   G4Material::GetMaterial("Galactic"),
						   "world_log",0,0,0);
  G4VPhysicalVolume* world_phys = new G4PVPlacement(0,G4ThreeVector(),"world",
                                            world_log,0,false,0);
  world_log->SetVisAttributes(G4VisAttributes::Invisible);

  G4Colour  grey(0.5, 0.5, 0.5) ;
  G4Colour red (1.0,0.0,0.0);
  G4Colour green (0.0,1.0,0.0);
  G4Colour blue(0.0, 0.0, 1.0);
  G4Colour magenta(1.0, 0.0, 1.0);

  //************************************************************************
  //Al casing
  //************************************************************************
  G4double casingDiameter = 70.0*mm;
  G4double casingThickness = 1.0*mm;

  G4Tubs* alCasingTube = new G4Tubs("alCasingTube",0,casingDiameter/2.,
				    casingHeight/2.,0,twopi);
  G4LogicalVolume* alCasingLog = new G4LogicalVolume(alCasingTube,
						     G4Material::GetMaterial("Aluminum"),
						     "alCasingLog",0,0,0);
  alCasingLog->SetVisAttributes(new G4VisAttributes(grey));
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    alCasingLog,"AlCasing",world_log,false,0);

  

  //************************************************************************
  //Vacuum inside the Al casing
  //************************************************************************
  G4double vacDiameter = casingDiameter - 2.0*casingThickness;
  G4double vacHeight = casingHeight - 2.0*casingThickness;
  
  G4Tubs* vacCasingTube = new G4Tubs("vacCasingTube",0,vacDiameter/2.,
				     vacHeight/2.,0,twopi);
  G4LogicalVolume* vacCasingLog = new G4LogicalVolume(vacCasingTube,
						      G4Material::GetMaterial("Galactic"),
						      "vacCasingLog",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    vacCasingLog,"CasingVacuum",alCasingLog,false,0);
  vacCasingLog->SetVisAttributes(new G4VisAttributes(grey));
  
  //************************************************************************
  // Ge block (dead layer)
  //************************************************************************
  G4double geDiameter = 60.0*mm;
  G4double geHeight = 60.0*mm;
  G4Tubs* geTube = new G4Tubs("geTube",0,geDiameter/2.,geHeight/2.,0,twopi);
  G4LogicalVolume* geLog = new G4LogicalVolume(geTube,
					       G4Material::GetMaterial("Germanium"),
					       "geLog",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    geLog,"DeadLayer",vacCasingLog,false,0);
  geLog->SetVisAttributes(new G4VisAttributes(magenta));
  
  
  //************************************************************************
  // Sensitive Ge block 
  //************************************************************************
  G4double DLThickness = 1.0*mm;
  G4double sensDiameter = geDiameter - 2.0*DLThickness;
  G4double sensHeight = geHeight - 1.0*DLThickness;
  G4Tubs* geSensTube = new G4Tubs("geSensTube",0,sensDiameter/2.,sensHeight/2.,
					   0,twopi);
  G4LogicalVolume* geSensLog = new G4LogicalVolume(geSensTube,
						   G4Material::GetMaterial("Germanium"),
						   "geSensLog",0,0,0);
  G4double tolerance = 0.1*micrometer;
  G4double zshift = -0.5*DLThickness+tolerance;
  new G4PVPlacement(0,G4ThreeVector(0,0,zshift),
		    geSensLog,"SensitiveGe",
		    geLog,false,0);
  geSensLog->SetVisAttributes(new G4VisAttributes(red));
  
  //************************************************************************
  // Inner hole
  //************************************************************************
  G4double holeHeight = 40.0*mm;
  G4double holeDiameter = 10.0*mm;
  
  G4Tubs* holeTubs = new G4Tubs("holeTube",0,holeDiameter/2.,holeHeight/2.,
				0,twopi);
  G4LogicalVolume* holeLog = new G4LogicalVolume(holeTubs,
						 G4Material::GetMaterial("Galactic"),
						 "holeLog",0,0,0);
  zshift = -sensHeight/2. + holeHeight/2. + tolerance;
  new G4PVPlacement(0,G4ThreeVector(0,0,zshift),
		    holeLog,"InnerBoreHole",geSensLog,false,0);
  holeLog->SetVisAttributes(new G4VisAttributes(magenta));

  //************************************************************************
  // block of water, if it is the case
  //************************************************************************
  if (geometryType == 3)
    {
      G4cout << "Using geometry#3" << G4endl;
      G4Tubs* waterTube = new G4Tubs("waterTube",0,waterDiameter/2.,waterHeight/2.,
				     0,twopi);
      G4LogicalVolume* waterLog = new G4LogicalVolume(waterTube,
						      G4Material::GetMaterial("HeavyWater"),
						      "waterLog",0,0,0);
      waterLog->SetVisAttributes(new G4VisAttributes(blue));
      new G4PVPlacement(0,G4ThreeVector(0,0,waterZPosition),
			waterLog,"WaterSource",world_log,false,0);

    }
  else
    G4cout << "Using geometry#2" << G4endl;


  //definisco lo scintillatore come un materiale sensibile
  G4SDManager* SDmanager=G4SDManager::GetSDMpointer();
  G4String nomescint = "scintillator";
  SensitiveDetector *pointerscint = new SensitiveDetector(nomescint);
  SDmanager->AddNewDetector(pointerscint);
  geSensLog->SetSensitiveDetector(pointerscint);
  

  return world_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // This function illustrates the possible ways to define materials using 
  // G4 database on G4Elements
  G4NistManager* manager = G4NistManager::Instance();
  manager->SetVerbose(0);
  //
  // define Elements
  //
  G4double z,a;

  G4Element* H  = manager->FindOrBuildElement(1);
  //G4Element* C  = manager->FindOrBuildElement(6);
  //G4Element* N  = manager->FindOrBuildElement(7);
  G4Element* O  = manager->FindOrBuildElement(8);
  //G4Element* Na = manager->FindOrBuildElement(11);
  G4Element* Al = manager->FindOrBuildElement(13);
  //G4Element* Si = manager->FindOrBuildElement(14);
  G4Element* Ge = manager->FindOrBuildElement(32);
  //G4Element* Sb = manager->FindOrBuildElement(51);
  //G4Element* I  = manager->FindOrBuildElement(53);
  //G4Element* Cs = manager->FindOrBuildElement(55);
  //G4Element* Pb = manager->FindOrBuildElement(82);
  //G4Element* Bi = manager->FindOrBuildElement(83);

  G4int   ncomponents;				     
  G4double abundance;				     
  //
  // define simple materials
  //
  G4double density;

  G4Material* metGe = new G4Material("Germanium",density=5.323*g/cm3,ncomponents=1);
  metGe->AddElement(Ge,abundance=100.0*perCent);


  G4Material* metAl = new G4Material("Aluminum",density=2.7*g/cm3,ncomponents=1);
  metAl->AddElement(Al,abundance=100.0*perCent);

  //This is water with artificial density of 3 g/cm3
  G4Material* heaWater = new G4Material("HeavyWater",density=3.0*g/cm3,ncomponents=2);
  heaWater->AddElement(O,1);
  heaWater->AddElement(H,2);

  //
  // examples of vacuum
  //

  density     = universe_mean_density;    //from PhysicalConstants.h
  G4double pressure    = 3.e-18*pascal;
  G4double temperature = 2.73*kelvin;
  new G4Material("Galactic", z=1., a=1.008*g/mole, density,
                             kStateGas,temperature,pressure);


  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetGeometryType(G4int gt)
{
  if (gt < 2 || gt > 3)
    {
      G4cout << "Valid geometry types are only 2 or 3 (geoetry #2 and geometry#3, respectively)" << G4endl;
      G4cout << "Command ignored" << G4endl;
      return;
    }
  geometryType = gt;
  return;
}
