
// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em0DetectorConstruction.cc,v 1.1 1999-01-08 16:32:35 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em0DetectorConstruction.hh"
#include "Em0DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4RunManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UnitsTable.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em0DetectorConstruction::Em0DetectorConstruction()
:pBox     (NULL), lBox (NULL),
 aMaterial(NULL)
{
  // create commands for interactive definition of the detector  
  detectorMessenger = new Em0DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em0DetectorConstruction::~Em0DetectorConstruction()
{ delete detectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* Em0DetectorConstruction::Construct()
{
  DefineMaterials();
  BoxSize = 20*m;
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em0DetectorConstruction::DefineMaterials()
{ 
 //This function illustrates the possible ways to define materials
 
G4String name, symbol;             //a=mass of a mole;
G4double a, z, density;            //z=mean number of protons;  

G4int ncomponents, natoms;
G4double fractionmass;

//
// define Elements
//

//a = 1.01*g/mole;
//G4Element* H = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);

//a = 14.01*g/mole;
//G4Element* N = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

//a = 16.00*g/mole;
//G4Element* O = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

//
// define materials
//

//density = 1.290*mg/cm3;
//G4Material* Air = new G4Material(name="Air"  , density, ncomponents=2);
//Air->AddElement(N, fractionmass=70.*perCent);
//Air->AddElement(O, fractionmass=30.*perCent);

//density = 70.8*mg/cm3;
//G4Material* H2l = new G4Material(name="H2liquid", density, ncomponents=1);
//H2l->AddElement(H, fractionmass=1.);
 
//density = 1.000*g/cm3;
//G4Material* H2O = new G4Material(name="Water", density, ncomponents=2);
//H2O->AddElement(H, natoms=2);
//H2O->AddElement(O, natoms=1);

//density = 1.390*g/cm3;
//a = 39.95*g/mole;
//G4Material* lAr = new G4Material(name="liquidArgon", z=18., a, density);

density = 2.700*g/cm3;
a = 26.98*g/mole;
G4Material* Al  = new G4Material(name="Aluminium"  , z=13., a, density);

//density = 7.870*g/cm3;
//a = 55.85*g/mole;
//G4Material* Fe  = new G4Material(name="Iron"       , z=26., a, density);

//density = 11.35*g/cm3;
//a = 207.19*g/mole;
//G4Material* Pb = new G4Material(name="Lead"        , z=82., a, density);

//density = 18.95*g/cm3;
//a = 238.03*g/mole;
//G4Material*  U = new G4Material(name="Uranium"     , z=82., a, density);


G4cout << *(G4Material::GetMaterialTable()) << endl;

  //default materials of the calorimeter
  aMaterial = Al;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
G4VPhysicalVolume* Em0DetectorConstruction::ConstructVolumes()
{    
  G4Box*
  sBox = new G4Box("Container",				//its name
                   BoxSize/2,BoxSize/2,BoxSize/2);	//its dimensions
                 
  lBox = new G4LogicalVolume(sBox,			//its shape
                             aMaterial,			//its material
                             "Container");		//its name
                                   
  pBox = new G4PVPlacement(0,				//no rotation
  			   G4ThreeVector(),		//at (0,0,0)
                           "Container",			//its name
                           lBox,			//its logical volume
                           NULL,			//its mother  volume
                           false,			//no boolean operation
                           0);				//copy number
                           
  //                                        
  // Visualization attributes
  //

  G4VisAttributes* visAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  visAtt->SetVisibility(true);
  lBox->SetVisAttributes(visAtt);
  
  //
  //always return the root volume
  //
  
  PrintParameters();
  return pBox;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em0DetectorConstruction::PrintParameters()
{
  G4cout << "\n The Box is " << G4BestUnit(BoxSize,"Length")
         << " of " << aMaterial->GetName() << endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em0DetectorConstruction::SetMaterial(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial)
     {aMaterial = pttoMaterial;
      lBox->SetMaterial(aMaterial); 
     }             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em0DetectorConstruction::SetSize(G4double value)
{
  BoxSize = value;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em0DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructVolumes());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
