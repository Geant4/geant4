// Em6DetectorConstruction.cc

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Em6DetectorConstruction.hh"
#include "Em6DetectorMessenger.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include "G4RunManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UnitsTable.hh"
#include "G4ios.hh"
#include "G4UserLimits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em6DetectorConstruction::Em6DetectorConstruction()
:nLtot(1),dLradl(1.),
 defaultMaterial(0),magField(0)  ,
 AbsorberLength(0.),AbsorberSizeXY(1.*meter),
 solidAbsorber(0),logicAbsorber(0),physiAbsorber(0),
 solidLayer(0)  ,logicLayer(0)  ,physiLayer(0)
{
  detectorMessenger = new Em6DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em6DetectorConstruction::~Em6DetectorConstruction()
{ delete detectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* Em6DetectorConstruction::Construct()
{
  DefineMaterials();
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em6DetectorConstruction::DefineMaterials()
{
  G4String name, symbol;
  G4double a, z, density;
  G4int ncomponents, natoms;
  G4double fractionmass;

  //
  // define few Elements
  //

//    a = 9.012182*g/mole;
//    G4Element* Be = new G4Element(name="Beryllium", symbol="Be", z=4., a);

  //
  // define materials
  //

  //Be
    a = 9.012182*g/mole;
    density = 1.848*g/cm3;
    G4Material* Be = new G4Material(name="Be", z=4., a, density);

  //Fe
    a = 55.85*g/mole;
    density = 7.87*g/cm3;
    G4Material* Fe = new G4Material(name="Fe", z=26., a, density);

  //choose material
  defaultMaterial = Fe;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* Em6DetectorConstruction::ConstructVolumes()
{
  G4double Radl = defaultMaterial->GetRadlen();

  G4double dL = dLradl*Radl;
  AbsorberLength = nLtot*dL;

  //
  // Absorber ( Absorber) block, cube of nLtot rad length
  //
  solidAbsorber = new G4Box("Absorber", // its name
	    AbsorberSizeXY/2,AbsorberSizeXY/2,AbsorberLength/2); // its x, y, z size

  logicAbsorber = new G4LogicalVolume(solidAbsorber,  //its solid
                                   defaultMaterial,	//its material
                                   "Absorber");		//its name
  physiAbsorber = new G4PVPlacement(0,               //no rotation
                                 G4ThreeVector(),   //at (0,0,0)
                                 "Absorber",         //its name
                                 logicAbsorber,      //its logical volume
                                 0,			        //its mother  volume
                                 false,			    //no boolean operation
                                 0);			    //copy number

  // subdivide in Layers in z-direction
  solidLayer = new G4Box("Layer",
	           AbsorberSizeXY/2,AbsorberSizeXY/2,dL/2); // x, y, z size
  logicLayer = new G4LogicalVolume(solidLayer,   //its solid
                              defaultMaterial,   //its material
                                       "Layer"); //its name
  if (nLtot >1)
    physiLayer = new G4PVReplica("Layer",		//its name
      				     logicLayer,	//its logical volume
      				     physiAbsorber,	//its mother
                         kZAxis,		//axis of replication
                         nLtot,	//number of replica
                         dL);	//width of replica
  else
    physiLayer = new G4PVPlacement(0, //no rotation
                                   G4ThreeVector(), //at (0,0,0)
                                   "Layer",     //its name
                                   logicLayer,  //its logical volume
                                   physiAbsorber,  //its mother  volume
                                   false,       //no boolean operation
                                   0);          //copy number
  G4cout
	<< "Absorber length in z-directions is " << nLtot*dLradl << " rad. length or "
    << G4BestUnit(AbsorberLength,"Length")
	<< " and SizeXY=" << G4BestUnit(AbsorberSizeXY,"Length") << '\n'
    << "Absorber material is " << defaultMaterial->GetName()
	<< " which has a rad. length of " << G4BestUnit(Radl,"Length") << '\n'
	<< "The Absorber is subdivided in z in nLtot=" << nLtot << " Layers" << '\n'
	<< " with Layer thickness dLradl=" << dLradl << " rad.leng."
	<< " or dL=" << G4BestUnit(dL,"Length") << '\n'
    << G4endl;

  G4VisAttributes* VisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  VisAtt->SetVisibility(true);
  logicAbsorber->SetVisAttributes(VisAtt);

  //
  //always return the physical World
  //
  return physiAbsorber;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em6DetectorConstruction::SetMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);
  if (pttoMaterial)
     {defaultMaterial = pttoMaterial;
      logicAbsorber->SetMaterial(defaultMaterial);
     }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em6DetectorConstruction::SetLBining(G4ThreeVector Value)
{
  nLtot = (G4int)Value(0); dLradl = Value(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em6DetectorConstruction::SetMagField(G4double fieldValue)
{
  //apply a global uniform magnetic field along Z axis
  G4FieldManager* fieldMgr
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  if(magField) delete magField;		//delete the existing magn field

  if(fieldValue!=0.)			// create a new one if non nul
  { magField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);
  } else {
    magField = 0;
    fieldMgr->SetDetectorField(magField);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em6DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
