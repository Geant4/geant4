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
// $Id: Tst50DetectorConstruction.cc,v 1.15 2003-02-10 11:21:04 guatelli Exp $
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
#include "G4UserLimits.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst50DetectorConstruction::Tst50DetectorConstruction(G4bool  max_Step)
:TargetMaterial(0),defaultMaterial(0),
 solidWorld(0),logicWorld(0),physiWorld(0),
 solidTarget(0),logicTarget(0),physiTarget(0)
 
{
  // default parameter values of the calorimeter
  TargetThickness = 0.1*mm;
  targetX=20. *cm;
  targetY=20. *cm;
  // create commands for interactive definition of the calorimeter  


 
theUserLimitsForTarget = NULL; 
 fUseUserLimits = max_Step;//non fa nulla se c'e false ,devo mettere true 
 theMaxStepInTarget = 0.000000001*micrometer;

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
a = 1.01*g/mole;
 density=8.3748e-5 *g/cm3; 
G4Material* HMat = new G4Material(name="Hydrogen", z=1., a, density);

density = 0.178e-3*g/cm3;
a = 4.0026*g/mole;
G4Material* HeMat = new G4Material(name="Helium", z=2., a, density);

a = 6.941*g/mole;
density=5.3400e-1*g/cm3; 
G4Material* LiMat = new G4Material(name="Lithium", z=3., a, density);

a = 9.012*g/mole;
density= 1.848*g/cm3; 
G4Material* BeMat = new G4Material(name="Beryllium", z=4., a, density);

a = 10.811*g/mole;
density= 2.37*g/cm3; 
G4Material* BMat = new G4Material(name="Boron", z=5., a, density);

a = 12.01*g/mole;
density= 2.00*g/cm3; 
G4Material* ACMat = new G4Material(name="AmorphousCarbon", z=6., a, density);

a = 12.01*g/mole;
density= 1.7*g/cm3; 
G4Material* GraphMat = new G4Material(name="Graphite", z=6., a, density);

a = 14.008*g/mole;
density= 1.165e-3*g/cm3; 
G4Material* NMat = new G4Material(name="Nitrogen", z=7., a, density);

a = 15.998*g/mole;
density= 1.331e-3*g/cm3; 
G4Material* OMat = new G4Material(name="Oxygen", z=8., a, density);

a = 18.998*g/mole;
density= 1.580e-3*g/cm3; 
G4Material* FMat = new G4Material(name="Fluorine", z=9., a, density);

a = 20.179*g/mole;
density= 8.385e-4*g/cm3; 
G4Material* NeMat = new G4Material(name="Neon", z=10., a, density);

a = 22.989*g/mole;
density= 9.71e-1*g/cm3; 
G4Material* NaMat = new G4Material(name="Sodium", z=11., a, density);

a = 24.312*g/mole;
density= 1.738*g/cm3; 
G4Material* MgMat = new G4Material(name="Magnesium", z=12., a, density);

a = 26.981*g/mole;
density= 2.6989*g/cm3; 
G4Material* AlMat = new G4Material(name="Aluminium", z=13., a, density);

a = 28.085*g/mole;
density= 2.336*g/cm3; 
G4Material* SiMat = new G4Material(name="Silicon", z=14., a, density);

a = 30.974*g/mole;
density= 2.20*g/cm3; 
G4Material* PMat = new G4Material(name="Phosphorus", z=15., a, density);

a = 32.065*g/mole;
density= 2.070*g/cm3; 
G4Material* SMat = new G4Material(name="Sulphur", z=16., a, density);

a =35.453 *g/mole;
density= 2.994e-3*g/cm3; 
G4Material* ClMat = new G4Material(name="Chlorine", z=17., a, density);

density = 1.390*g/cm3;
a = 39.95*g/mole;
G4Material* lArMat = new G4Material(name="liquidArgon", z=18., a, density);

a = 39.098*g/mole;
density= 8.62e-1*g/cm3;
G4Material* KMat = new G4Material(name="Potassium", z=19., a, density);

a = 40.078*g/mole;
density= 1.55*g/cm3;
G4Material* CaMat = new G4Material(name="Calcium", z=20., a, density);

a = 44.956*g/mole;
density= 2.989*g/cm3;
G4Material* ScMat = new G4Material(name="Scandium", z=21., a, density);

a = 47.88*g/mole;
density= 4.54*g/cm3;
G4Material* TiMat = new G4Material(name="Titanium", z=22., a, density);

a = 50.9415*g/mole;
density= 6.11*g/cm3;
G4Material* VMat = new G4Material(name="Vanadium", z=23., a, density);

a = 51.996*g/mole;
density= 7.18*g/cm3;
G4Material* CrMat = new G4Material(name="Chromium", z=24., a, density);

a = 54.9380*g/mole;
density= 7.44*g/cm3;
G4Material* MnMat = new G4Material(name="Manganese", z=25., a, density);

a = 55.845*g/mole;
density= 7.874*g/cm3;
G4Material* FeMat = new G4Material(name="Iron", z=26., a, density);

a = 58.9332*g/mole;
density= 8.9*g/cm3;
G4Material* CoMat = new G4Material(name="Cobalt", z=27., a, density);

a = 58.6934*g/mole;
density= 8.902*g/cm3;
G4Material* NiMat = new G4Material(name="Nickel", z=28., a, density);

density = 8.96*g/cm3;
a = 63.546*g/mole;
G4Material* CuMat = new G4Material(name="Copper"  , z=29., a, density);

density = 7.133*g/cm3;
a = 65.409*g/mole;
G4Material* ZnMat = new G4Material(name="Zinc"  , z=30., a, density);

density = 5.904*g/cm3;
a = 69.723*g/mole;
G4Material* GaMat = new G4Material(name="Gallium"  , z=31., a, density);

density = 5.323*g/cm3;
a = 72.64*g/mole;
G4Material* GeMat = new G4Material(name="Germanium"  , z=32., a, density);

density = 5.727*g/cm3;
a = 74.92160*g/mole;
G4Material* AsMat = new G4Material(name="Arsenic"  , z=33., a, density);

density = 4.819*g/cm3;
a = 78.96*g/mole;
G4Material* SeMat = new G4Material(name="Selenium"  , z=34., a, density);

density = 7.07218e-3*g/cm3;
a = 79.904*g/mole;
G4Material* BrMat = new G4Material(name="Bromine"  , z=35., a, density);

density = 3.47832e-3*g/cm3;
a = 83.798*g/mole;
G4Material* KrMat = new G4Material(name="Krypton"  , z=36., a, density);

density = 1.532*g/cm3;
a = 85.4678*g/mole;
G4Material* RbMat = new G4Material(name="Rubidium"  , z=37., a, density);

density = 2.54*g/cm3;
a = 87.62*g/mole;
G4Material* SrMat = new G4Material(name="Strontium"  , z=38., a, density);

density = 4.469*g/cm3;
a = 88.90585*g/mole;
G4Material* YMat = new G4Material(name="Yttrium"  , z=39., a, density);

density = 6.506*g/cm3;
a =91.224*g/mole;
G4Material* ZrMat = new G4Material(name="Zirconium"  , z=40., a, density);

density = 8.57*g/cm3;
a =92.90638*g/mole;
G4Material* NbMat = new G4Material(name="Niobium"  , z=41., a, density);

density = 10.22*g/cm3;
a =95.94*g/mole;
G4Material* MoMat = new G4Material(name="Molybdenum"  , z=42., a, density);

density = 11.5*g/cm3;
a =98.*g/mole;
G4Material* TcMat = new G4Material(name="Technetium"  , z=43., a, density);

density = 12.41*g/cm3;
a =101.07*g/mole;
G4Material* RuMat = new G4Material(name="Ruthenium"  , z=44., a, density);

density = 12.41*g/cm3;
a =102.905*g/mole;
G4Material* RhMat = new G4Material(name="Rhodium"  , z=45., a, density);

density =12.02*g/cm3;
a =106.42*g/mole;
 G4Material* PdMat = new G4Material(name="Palladium"  , z=46., a, density);

density =10.5*g/cm3;
a =107.8682*g/mole;
 G4Material* AgMat = new G4Material(name="Silver"  , z=47., a, density);

density =8.65*g/cm3;
a =112.411*g/mole;
 G4Material* CdMat = new G4Material(name="Cadmium"  , z=48., a, density);

density =7.31*g/cm3;
a =114.818*g/mole;
 G4Material* InMat = new G4Material(name="Indium"  , z=49., a, density);

density =7.31*g/cm3;
a =118.710*g/mole;
 G4Material* SnMat = new G4Material(name="Tin"  , z=50., a, density);

density =6.691*g/cm3;
a =121.760*g/mole;
 G4Material* SbMat = new G4Material(name="Antimony"  , z=51., a, density);

density =6.24*g/cm3;
a =127.60*g/mole;
 G4Material* TeMat = new G4Material(name="Tellurium"  , z=52., a, density);

density =4.93*g/cm3;
a =126.90447*g/mole;
 G4Material* IMat = new G4Material(name="Iodine"  , z=53., a, density);


density =5.48536e-3*g/cm3;
a =131.293*g/mole;
 G4Material* XeMat = new G4Material(name="Xenon"  , z=54., a, density);

density =1.873*g/cm3;
a =132.90545*g/mole;
 G4Material* CsMat = new G4Material(name="Cesium"  , z=55., a, density);

density =3.5*g/cm3;
a =137.327*g/mole;
 G4Material* BaMat = new G4Material(name="Barium"  , z=56., a, density);


density =6.154*g/cm3;
a =138.9055*g/mole;
 G4Material* LaMat = new G4Material(name="Lanthanum" , z=57., a, density);

density =6.657*g/cm3;
a =140.116*g/mole;
 G4Material* CeMat = new G4Material(name="Cerium" , z=58., a, density);

density =6.71*g/cm3;
a =140.90765*g/mole;
 G4Material* PrMat = new G4Material(name="Praseodymium" , z=59., a, density);

density =6.9*g/cm3;
a =144.24*g/mole;
 G4Material* NdMat = new G4Material(name="Neodymium" , z=60., a, density);

density =7.22*g/cm3;
a =145.*g/mole;
 G4Material* PmMat = new G4Material(name="Promethium" , z=61., a, density);

density =7.46*g/cm3;
a =150.36*g/mole;
 G4Material* SmMat = new G4Material(name="Samarium"  , z=62., a, density);

density =5.243*g/cm3;
a =151.964*g/mole;
 G4Material* EuMat = new G4Material(name="Europium"  , z=63., a, density);


density =7.9004*g/cm3;
a =157.25*g/mole;
 G4Material* GdMat = new G4Material(name="Gadolinium"  , z=64., a, density);


density =8.229*g/cm3;
a =158.92534*g/mole;
 G4Material* TbMat = new G4Material(name="Terbium"  , z=65., a, density);

density =8.55*g/cm3;
a =162.500*g/mole;
 G4Material* DyMat = new G4Material(name="Dysprosium"  , z=66., a, density);

density =8.795*g/cm3;
a =164.93032*g/mole;
 G4Material* HoMat = new G4Material(name="Holmium"  , z=67., a, density);

density =9.066*g/cm3;
a =167.259*g/mole;
 G4Material* ErMat = new G4Material(name="Erbium"  , z=68., a, density);


density =9.321*g/cm3;
a =168.93421*g/mole;
 G4Material* TmMat = new G4Material(name="Thulium"  , z=69., a, density);

density =6.73*g/cm3;
a =173.04*g/mole;
 G4Material* YbMat = new G4Material(name="Ytterbium"  , z=70., a, density);

density =9.84*g/cm3;
a =174.967*g/mole;
 G4Material* LuMat = new G4Material(name="Lutetium"  , z=71., a, density);

density =13.31*g/cm3;
a =178.49*g/mole;
 G4Material* HfMat = new G4Material(name="Hafnium"  , z=72., a, density);

density =16.6540*g/cm3;
a =180.9479*g/mole;
 G4Material* TaMat = new G4Material(name="Tantalum"  , z=73., a, density);

density =19.3*g/cm3;
a =183.84*g/mole;
 G4Material* WMat = new G4Material(name="Tungsten"  , z=74., a, density);

density =21.020*g/cm3;
a =186.207*g/mole;
 G4Material* ReMat = new G4Material(name="Rhenium"  , z=75., a, density);

density =22.57*g/cm3;
a =190.23*g/mole;
 G4Material* OsMat = new G4Material(name="Osmium"  , z=76., a, density);

density =22.42*g/cm3;
a =192.217*g/mole;
 G4Material* IrMat = new G4Material(name="Iridium"  , z=77., a, density);

density =21.45*g/cm3;
a =195.078*g/mole;
 G4Material* PtMat = new G4Material(name="Platinum"  , z=78., a, density);

density = 19.32*g/cm3;
a = 196.966*g/mole;
G4Material* AuMat = new G4Material(name="Gold"     , z=79., a, density);

density =13.54*g/cm3;
a =200.59*g/mole;
G4Material* HgMat = new G4Material(name="Mercury"     , z=80., a, density);

density =11.72*g/cm3;
a =204.3833*g/mole;
G4Material* TlMat = new G4Material(name="Thallium"     , z=81., a, density);

density = 11.35*g/cm3;
a = 207.19*g/mole;
G4Material* PbMat = new G4Material(name="Lead"     , z=82., a, density);



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


  TargetMaterial =lArMat;


  defaultMaterial  = Vacuum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
G4VPhysicalVolume* Tst50DetectorConstruction::ConstructWorld()
{
  // complete the Calor parameters definition 
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
  

  solidTarget=0; logicTarget=0; physiTarget=0;  
  
  if (TargetThickness > 0.) 
    { solidTarget = new G4Box("Target",		//its name
                          targetX/2,targetY/2, TargetThickness/2 ); 
                          
      logicTarget = new G4LogicalVolume(solidTarget,    //its solid
      			                  TargetMaterial, //its material
      			                  "Target");      //its name
      			                  
  // create UserLimits
  if (theUserLimitsForTarget != NULL) delete theUserLimitsForTarget;
  theUserLimitsForTarget = new G4UserLimits(//DBL_MAX,  //step max
					      //DBL_MAX,  // track max
					      theMaxStepInTarget);
					     
					
 

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
        
         << TargetThickness/mm<<" mm" <<G4endl; 

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
G4double Tst50DetectorConstruction::GetDensity()
{
  G4double Density=(TargetMaterial->GetDensity());
  return Density;
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

void  Tst50DetectorConstruction::SetMaxStepInTarget(G4double value)
{ 
  theMaxStepInTarget = value; 
  if (theUserLimitsForTarget != NULL) 
  {
    theUserLimitsForTarget->SetMaxAllowedStep(value);
  }
}

void  Tst50DetectorConstruction::UseUserLimits(G4bool isUse) 
{
  fUseUserLimits = isUse;
  if ( fUseUserLimits && (theUserLimitsForTarget!= NULL)) 
  {logicTarget->SetUserLimits(theUserLimitsForTarget);
  }    
} 
