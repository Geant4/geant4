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
// $Id: Tst50DetectorConstruction.cc,v 1.32 2010-06-25 09:46:57 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// author: Susanna Guatelli (guatelli@ge.infn.it)
// 
// History:
// -----------
// 17 May  2003   S. Guatelli   1st implementation
//
// -------------------------------------------------------------------

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
#include "G4UserLimits.hh"
#include "G4ios.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

Tst50DetectorConstruction::Tst50DetectorConstruction()
  :isRegisteredUserLimits(false), hydrogen(0),beryllium(0),graphite(0), 
   magnesium(0), aluminium(0),silicon(0),liquidArgon(0),titanium(0),iron(0),
   cobalt(0),nickel(0),indium(0), tin(0),
    copper(0),zinc(0),
   gallium(0),germanium(0), zirconium(0),
   molybdenium(0), silver(0),cadmium(0),
   cesium(0),samarium(0), ytterbium(0),tantalum(0),tungsten(0),
   gold(0),
   lead(0),uranium(0), water(0), quartz(0), air(0),vacuum(0),nytrogen(0),
   targetMaterial(0),defaultMaterial(0),
   solidWorld(0),logicWorld(0),physiWorld(0),
   solidTarget(0),logicTarget(0),physiTarget(0), 
   targetSD(0)
{
  // default parameter values of the target
  targetThickness = 0.1*mm;
  targetX=20. *cm;
  targetY=20. *cm;
 
  theUserLimitsForTarget = 0; 
  fUseUserLimits = false;
  theMaxStepInTarget = 1.*micrometer;

  messenger = new Tst50DetectorMessenger(this);
}

Tst50DetectorConstruction::~Tst50DetectorConstruction()
{
  delete messenger;
  delete nytrogen;
  delete vacuum;
  delete air;
  delete quartz;
  delete water;
  delete uranium; 
  delete lead;
  delete gold;
  delete tungsten;
  delete tantalum;
  delete ytterbium;
  delete samarium;
  delete cesium;
  delete tin;
  delete indium;
  delete cadmium;
  delete silver;
  delete molybdenium;
  delete zirconium;
  delete germanium;
  delete gallium;
  delete zinc;
  delete copper;
  delete nickel;
  delete cobalt;
  delete iron; 
  delete titanium;
  delete liquidArgon;
  delete silicon;
  delete aluminium;
  delete  magnesium;
  delete graphite;
  delete beryllium;
  delete hydrogen;
}

G4VPhysicalVolume* Tst50DetectorConstruction::Construct()
{
  DefineMaterials();
  return ConstructWorld();
}

void Tst50DetectorConstruction::DefineMaterials()
{ 
 
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
  a = 1.01*g/mole;
  density = 8.3748e-5 *g/cm3; 
  hydrogen = new G4Material(name="Hydrogen", z=1., a, density);

  a = 9.012*g/mole;
  density = 1.848*g/cm3; 
  beryllium = new G4Material(name="Beryllium", z=4., a, density);

  a = 12.01*g/mole;
  density = 1.7*g/cm3; 
  graphite = new G4Material(name="Graphite", z=6., a, density);
  graphite -> GetIonisation()->SetMeanExcitationEnergy(78.*eV);
 
  a = 16.00*g/mole;
  density = 0.00133*g/cm3;
  ossigeno = new G4Material(name="Oxygen", z=8.,a,density);
  ossigeno->GetIonisation()->SetMeanExcitationEnergy(95.*eV);

  a = 24.312*g/mole;
  density = 1.738*g/cm3; 
  magnesium = new G4Material(name="Magnesium", z=12., a, density);
  

  a = 26.981*g/mole;
  density = 2.6989*g/cm3; 
  aluminium = new G4Material(name="Aluminium", z=13., a, density);
  aluminium->GetIonisation()->SetMeanExcitationEnergy(166.0*eV);

  a = 28.085*g/mole;
  density = 2.336*g/cm3; 
  silicon = new G4Material(name="Silicon", z=14., a, density);
  silicon->GetIonisation()->SetMeanExcitationEnergy(173.0*eV);


  density = 1.390*g/cm3;
  a = 39.95*g/mole;
  liquidArgon = new G4Material(name="liquidArgon", z=18., a, density);


  density = 4.54 *g/cm3;
  a = 47.867 *g/mole;
  titanium = new G4Material(name="Titanium", z=22., a, density);


  a = 55.845*g/mole;
  density = 7.874*g/cm3;
  iron = new G4Material(name="Iron", z=26., a, density);
  iron -> GetIonisation()->SetMeanExcitationEnergy(286.*eV);
 
  a = 58.933*g/mole;
  density = 8.9 *g/cm3;
  cobalt = new G4Material(name="Cobalt", z=27., a, density);

  a = 58.69*g/mole;
  density = 8.902 *g/cm3;
  nickel = new G4Material(name="Nickel", z=28., a, density);

  a = 63.546*g/mole;
  density = 8.96 *g/cm3;
  copper = new G4Material(name="Copper", z=29., a, density);

  a = 65.409*g/mole;
  density = 7.133 *g/cm3;
  zinc = new G4Material(name="Zinc", z=30., a, density);
 
  density = 5.904*g/cm3;
  a = 69.723*g/mole;
  gallium = new G4Material(name="Gallium", z=31., a, density);

  density = 5.323*g/cm3;
  a = 72.64*g/mole;
  germanium = new G4Material(name="Germanium", z=32., a, density);

  density = 6.506*g/cm3;
  a = 91.224 *g/mole;
  zirconium = new G4Material(name="Zirconium", z=40., a, density);

  density = 10.22 *g/cm3;
  a = 95.94 *g/mole;
  molybdenium = new G4Material(name="Molybdenium", z=42., a, density);

  density = 10.5*g/cm3;
  a = 107.8682*g/mole;
  silver = new G4Material(name="Silver", z=47., a, density);

  density = 8.65*g/cm3;
  a = 112.411*g/mole;
  cadmium = new G4Material(name="Cadmium", z=48., a, density);
 
  density = 7.310 *g/cm3;
  a = 114.818*g/mole;
  indium = new G4Material(name="Indium", z=49., a, density);

  density = 7.310 *g/cm3;
  a = 118.71*g/mole;
  tin = new G4Material(name="Tin", z=50., a, density);

  density = 1.873*g/cm3;
  a = 132.90545*g/mole;
  cesium = new G4Material(name="Cesium", z=55., a, density);

  density = 7.46*g/cm3;
  a = 150.36*g/mole;
  samarium = new G4Material(name="Samarium", z=62., a, density);

  density = 6.73*g/cm3;
  a = 173.04*g/mole;
  ytterbium = new G4Material(name="Ytterbium", z=70., a, density);
  
  density = 16.65 *g/cm3;
  a = 180.9947*g/mole;
  tantalum = new G4Material(name="Tantalum", z=73., a, density);

  density = 19.3 *g/cm3;
  a = 183.85*g/mole;
  tungsten = new G4Material(name="Tungsten", z=74., a, density);

  density = 19.32*g/cm3;
  a = 196.966*g/mole;
  gold = new G4Material(name="Gold", z=79., a, density);

  density = 11.35*g/cm3;
  a = 207.19*g/mole;
  lead = new G4Material(name="Lead", z=82., a, density);

  density = 18.95*g/cm3;
  a = 238.*g/mole;
  uranium = new G4Material(name="Uranium", z=92., a, density);


  //
  // define a material from elements.   case 1: chemical molecule
  //
 
  // density = 1.000*g/cm3;
//   water = new G4Material(name="Water", density, ncomponents=2);
//   water->SetChemicalFormula("H_2O");
//   water->AddElement(H, natoms=2);
//   water->AddElement(O, natoms=1);
//   water->GetIonisation()->SetMeanExcitationEnergy(75.0*eV);

  density = 2.200*g/cm3;
  quartz = new G4Material(name="Quartz", density, ncomponents=2);
  quartz->AddElement(Si, natoms=1);
  quartz->AddElement(O , natoms=2);

  density = 1.290*mg/cm3;
  air = new G4Material(name="Air"  , density, ncomponents=2);
  air->AddElement(N, fractionmass=0.7);
  air->AddElement(O, fractionmass=0.3);
  air->GetIonisation()->SetMeanExcitationEnergy(85.7*eV);

  //
  // examples of vacuum
  //

  density = universe_mean_density;
  pressure = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  vacuum = new G4Material(name="Galactic", z=1., a=1.01*g/mole,
				      density,kStateGas,temperature,pressure);

  density = 1.16*mg/cm3;
  a = 14.*g/mole;
  nytrogen = new G4Material(name="Nytrogen", z=7., a, density);
  nytrogen->GetIonisation()->SetMeanExcitationEnergy(82.0*eV);

 
  nistMan = G4NistManager::Instance();
  nistMan->SetVerbose(2);
  water_nist = nistMan -> FindOrBuildMaterial("G4_WATER");
  
  targetMaterial = water_nist;
  defaultMaterial  = vacuum;

}
  
G4VPhysicalVolume* Tst50DetectorConstruction::ConstructWorld()
{
  //     
  // World
  //
  solidWorld = new G4Box("World",				//its name
			 200000.*m,200000.*m,200000.*m);	//its size
                         
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
 
  // Slab filled with absorber material ...
  solidTarget = 0; 
  logicTarget = 0; 
  physiTarget = 0;  
  
  if (targetThickness > 0.) 
    { solidTarget = new G4Box("Target",		//its name
			      targetX/2,targetY/2, targetThickness/2 ); 
                          
    logicTarget = new G4LogicalVolume(solidTarget,    //its solid
				      targetMaterial, //its material
				      "Target");      //its name
      			                  
    // create UserLimits
    if (theUserLimitsForTarget != 0) delete theUserLimitsForTarget;
    theUserLimitsForTarget = new G4UserLimits(theMaxStepInTarget);

    // attach UserLimits   
    if (fUseUserLimits) {
      logicTarget->SetUserLimits(theUserLimitsForTarget);
    }



    physiTarget = new G4PVPlacement(0,		   //no rotation
				    G4ThreeVector(0.,0.,0.),  //its position
				    "Target",        //its name
				    logicTarget,     //its logical volume
				    physiWorld,        //its mother
				    false,             //no boulean operat
				    0);                //copy number
                                        
    }
  

  //  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  /*  
  if (targetSD == 0)
    {
      G4String targetSD_name = "target";
      targetSD = new Tst50TrackerSD( targetSD_name,this );
      SDman -> AddNewDetector( targetSD );
      G4cout<<"SD initialised in detector"<<G4endl;
    }
     logicTarget -> SetSensitiveDetector( targetSD );
  */
   
  // Visualization attributes
  //
  G4VisAttributes* simpleBoxVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt -> SetVisibility(true);
  logicWorld -> SetVisAttributes(simpleBoxVisAtt);
  logicTarget -> SetVisAttributes(simpleBoxVisAtt);
 
  //
  //always return the physical World
  //
  PrintParameters();
  return physiWorld;
}


void Tst50DetectorConstruction::PrintParameters()
{
  G4cout << " Target zdimension is: "
        
         << targetThickness/mm<<" mm" <<G4endl; 

  G4cout<<" Target xdimension is: "<<targetX/mm<<" mm"<<G4endl;
  G4cout<<" Target ydimension is: "<<targetY/mm<<" mm"<<G4endl;
  G4cout<<"Target Material is: "<<targetMaterial->GetName()<<G4endl; 


}

void Tst50DetectorConstruction::SetTargetMaterial(G4String materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if (pttoMaterial)
    {targetMaterial = pttoMaterial;
    logicTarget->SetMaterial(pttoMaterial); 
    PrintParameters();
    }             
}
G4double Tst50DetectorConstruction::GetDensity()
{
  G4double Density=(targetMaterial->GetDensity());
  return Density;
}
G4String Tst50DetectorConstruction::GetMaterialName()
{
  G4String name;
  name=targetMaterial->GetName();
  return name; 

}


void Tst50DetectorConstruction::SetTargetThickness(G4double val)
{ 
 
  // change Target thickness
  targetThickness = val;
 
}  
G4double Tst50DetectorConstruction::GetTargetThickness()
{ 
 
  // Get target thickness
  return  targetThickness;
 
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

void  Tst50DetectorConstruction::SetMaxStepInTarget(G4double value)
{ 
  theMaxStepInTarget = value; 
  if (theUserLimitsForTarget != 0) 
    {
      theUserLimitsForTarget->SetMaxAllowedStep(value);
    
    }
}

void  Tst50DetectorConstruction::UseUserLimits(G4bool isUse) 
{
  fUseUserLimits = isUse;
  if ( fUseUserLimits && (theUserLimitsForTarget!= 0)) 
    {logicTarget->SetUserLimits(theUserLimitsForTarget);
    }    
} 
void  Tst50DetectorConstruction::SetUserLimits(G4bool isUse) 
{
  if (isRegisteredUserLimits == false) {

    fUseUserLimits = isUse;
 
    if( fUseUserLimits && (theUserLimitsForTarget!= 0))   
      {logicTarget->SetUserLimits(theUserLimitsForTarget);
      } 
        
    isRegisteredUserLimits = true; }
  else 
    {
      G4cout<< "UseLimits is registered!  "<<G4endl;}  
}
