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
// $Id: Tst50DetectorConstruction.cc,v 1.5 2003-01-07 15:29:39 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Tst50DetectorConstruction.hh"
#include "Tst50DetectorMessenger.hh"

#include "Tst50TrackerSD.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst50DetectorConstruction::Tst50DetectorConstruction()
:TargetMaterial(0),defaultMaterial(0),
 solidWorld(0),logicWorld(0),physiWorld(0),
 solidTarget(0),logicTarget(0),physiTarget(0)
 
{
  // default parameter values of the calorimeter
  TargetThickness = 5.*mm;
  targetX=20. *cm;
  targetY=20. *cm;
  // create commands for interactive definition of the calorimeter  
  detectorMessenger = new Tst50DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst50DetectorConstruction::~Tst50DetectorConstruction()
{ delete detectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* Tst50DetectorConstruction::Construct()
{
  DefineMaterials();
  return ConstructWorld();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst50DetectorConstruction::DefineMaterials()
{ 
 //This function illustrates the possible ways to define materials
 
G4String name, symbol;             //a=mass of a mole;
G4double a, z, density;            //z=mean number of protons;  
G4int iz, n;                       //iz=number of protons  in an isotope; 
                                   // n=number of nucleons in an isotope;

G4int ncomponents, natoms;
G4double abundance, fractionmass;
G4double temperature, pressure;

//
// define Elements
//

a = 1.01*g/mole;
G4Element* H  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);

a = 12.01*g/mole;
G4Element* C  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);

a = 14.01*g/mole;
G4Element* N  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

a = 16.00*g/mole;
G4Element* O  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

a = 28.09*g/mole;
G4Element* Si = new G4Element(name="Silicon",symbol="Si" , z= 14., a);

//
// define an Element from isotopes, by relative abundance 
//

G4Isotope* U5 = new G4Isotope(name="U235", iz=92, n=235, a=235.01*g/mole);
G4Isotope* U8 = new G4Isotope(name="U238", iz=92, n=238, a=238.03*g/mole);

G4Element* U  = new G4Element(name="enriched Uranium",symbol="U",ncomponents=2);
U->AddIsotope(U5, abundance= 90.*perCent);
U->AddIsotope(U8, abundance= 10.*perCent);

//
// define simple materials
//



density = 1.390*g/cm3;
a = 39.95*g/mole;
G4Material* lAr = new G4Material(name="liquidArgon", z=18., a, density);

density = 11.35*g/cm3;
a = 207.19*g/mole;
G4Material* Pb = new G4Material(name="Lead"     , z=82., a, density);

//
// define a material from elements.   case 1: chemical molecule
//
 
density = 1.000*g/cm3;
G4Material* H2O = new G4Material(name="Water", density, ncomponents=2);
H2O->AddElement(H, natoms=2);
H2O->AddElement(O, natoms=1);


density = 2.200*g/cm3;
G4Material* SiO2 = new G4Material(name="quartz", density, ncomponents=2);
SiO2->AddElement(Si, natoms=1);
SiO2->AddElement(O , natoms=2);

//
// define a material from elements.   case 2: mixture by fractional mass
//

density = 1.290*mg/cm3;
G4Material* Air = new G4Material(name="Air"  , density, ncomponents=2);
Air->AddElement(N, fractionmass=0.7);
Air->AddElement(O, fractionmass=0.3);

//
// define a material from elements and/or others materials (mixture of mixtures)
//

//
// examples of gas in non STP conditions
//

density     = 27.*mg/cm3;
pressure    = 50.*atmosphere;
temperature = 325.*kelvin;
G4Material* CO2 = new G4Material(name="CarbonicGas", density, ncomponents=2,
                                     kStateGas,temperature,pressure);
CO2->AddElement(C, natoms=1);
CO2->AddElement(O, natoms=2);
 
density     = 0.3*mg/cm3;
pressure    = 2.*atmosphere;
temperature = 500.*kelvin;
G4Material* steam = new G4Material(name="WaterSteam", density, ncomponents=1,
                                      kStateGas,temperature,pressure);
steam->AddMaterial(H2O, fractionmass=1.);

//
// examples of vacuum
//

density     = universe_mean_density;
pressure    = 3.e-18*pascal;
temperature = 2.73*kelvin;
G4Material* Vacuum = new G4Material(name="Galactic", z=1., a=1.01*g/mole,
                                    density,kStateGas,temperature,pressure);

density     = 1.e-5*g/cm3;
pressure    = 2.e-2*bar;
temperature = STP_Temperature;         //from PhysicalConstants.h


  TargetMaterial = Pb;


  defaultMaterial  = Vacuum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
G4VPhysicalVolume* Tst50DetectorConstruction::ConstructWorld()
{
  // complete the Calor parameters definition 
  
  WorldSizeX= 20.*m;
    WorldSizeYZ=20.*m;
  //     
  // World
  //
  solidWorld = new G4Box("World",				//its name
                   WorldSizeX/2,WorldSizeYZ/2,WorldSizeYZ/2);	//its size
                         
  logicWorld = new G4LogicalVolume(solidWorld,		//its solid
                                   defaultMaterial,	//its material
                                   "World");		//its name
                                   
  physiWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 "World",		//its name
                                 logicWorld,		//its logical volume
                                 0,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number
  

  solidTarget=0; logicTarget=0; physiTarget=0;  
  
  if (TargetThickness > 0.) 
    { solidTarget = new G4Box("Target",		//its name
                          targetX/2,targetY/2, TargetThickness/2 ); 
                          
      logicTarget = new G4LogicalVolume(solidTarget,    //its solid
      			                  TargetMaterial, //its material
      			                  "Target");      //its name
      			                  
      physiTarget = new G4PVPlacement(0,		   //no rotation
      		    G4ThreeVector(0.,0.,0.),  //its position
                                        "Target",        //its name
                                        logicTarget,     //its logical volume
                                        physiWorld,        //its mother
                                        false,             //no boulean operat
                                        0);                //copy number
                                        
    }
  

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  if (pTargetSD ==0)
    {
  G4String targetSD_name = "target";
  pTargetSD = new Tst50TrackerSD( targetSD_name  );
  if(pTargetSD)
   { SDman->AddNewDetector( pTargetSD );
  logicTarget->SetSensitiveDetector( pTargetSD );}
    }

 
  // Visualization attributes
  //
  
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  logicWorld->SetVisAttributes(simpleBoxVisAtt);
  logicTarget ->SetVisAttributes(simpleBoxVisAtt);
 
  //
  //always return the physical World
  //
 PrintParameters();
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst50DetectorConstruction::PrintParameters()
{
  G4cout << " Target zdimension is: "
        
         << TargetThickness/mm <<G4endl; 

  G4cout<<" Target xdimension is: "<<targetX/mm<<" mm"<<G4endl;
 G4cout<<" Target ydimension is: "<<targetY/mm<<" mm"<<G4endl;
 G4cout<<"Target Material is: "<<TargetMaterial->GetName()<<G4endl; 
    }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst50DetectorConstruction::SetTargetMaterial(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial)
     {TargetMaterial = pttoMaterial;
      logicTarget->SetMaterial(pttoMaterial); 
      PrintParameters();
     }             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst50DetectorConstruction::SetTargetThickness(G4double val)
{
  // change Target thickness
  TargetThickness = val;
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Tst50DetectorConstruction::SetTargetX(G4double valX)
{
  targetX=valX;

}
void Tst50DetectorConstruction::SetTargetY(G4double valY) 
{
  targetY=valY;

}



void Tst50DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructWorld());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
