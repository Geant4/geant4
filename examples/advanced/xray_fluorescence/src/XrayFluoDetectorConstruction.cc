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
// $Id: XrayFluoDetectorConstruction.cc
// GEANT4 tag $Name: xray_fluo-V03-02-00
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
//
// -------------------------------------------------------------------

#include "XrayFluoDetectorConstruction.hh"
#include "XrayFluoDetectorMessenger.hh"
#include "XrayFluoHPGeSD.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "G4PVReplica.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoDetectorConstruction::XrayFluoDetectorConstruction()
  : DeviceSizeX(0),DeviceSizeY(0),DeviceThickness(0),
    solidWorld(0),logicWorld(0),physiWorld(0),
    solidHPGe(0),logicHPGe(0),physiHPGe(0),
    solidSample (0),logicSample(0),physiSample (0),
    solidDia1(0),logicDia1(0),physiDia1(0),
    solidDia3(0),logicDia3(0),physiDia3(0),
    solidOhmicPos(0),logicOhmicPos(0), physiOhmicPos(0),
    solidOhmicNeg(0),logicOhmicNeg(0), physiOhmicNeg(0),
    solidPixel(0),logicPixel(0), physiPixel(0),
    OhmicPosMaterial(0), OhmicNegMaterial(0),
    pixelMaterial(0),sampleMaterial(0),
    Dia1Material(0),Dia3Material(0),
    defaultMaterial(0)
  ,HPGeSD(0)
  
{ 
  NbOfPixelRows     =  1;
  NbOfPixelColumns  =  1;
  NbOfPixels        =  NbOfPixelRows*NbOfPixelColumns;
    PixelSizeXY       = 5 * cm;
  PixelThickness =  3.5 * mm;
  ContactSizeXY     = 0.005*mm;
  SampleThickness = 2.5 * mm;
  SampleSizeXY = 3. * cm;
  Dia1Thickness = 1. *mm;
  Dia3Thickness = 1. *mm;
  Dia1SizeXY = 3. *cm;
  Dia3SizeXY = 3. *cm;


  DiaInnerSize = 2.99 * cm;


  OhmicNegThickness = 0.005*mm;
  OhmicPosThickness = 0.005*mm;
  ThetaHPGe = 135. * deg;
  Dia3Dist =  66.5 * mm;
  Dia3InnerSize = 1. * mm;
  PhiHPGe = 225. * deg;
  ThetaDia1 = 135. * deg;
  ThetaDia3 = 180. * deg;
  PhiDia3 = 90. * deg;
  DistDia = 66.5 * mm;
  DistDe =DistDia+ (Dia1Thickness
		    +PixelThickness)/2+OhmicPosThickness ;
  PhiDia1 = 90. * deg;
  AlphaDia1 = 225. * deg;
  AlphaDia3 = 180. * deg;
  PixelCopyNb=0;
  ComputeApparateParameters();
  
  // create commands for interactive definition of the apparate
  
  detectorMessenger = new XrayFluoDetectorMessenger(this);
  G4cout << "XrayFluoDetectorConstruction created" << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoDetectorConstruction::~XrayFluoDetectorConstruction()
{ 
  delete detectorMessenger;
  G4cout << "XrayFluoDetectorConstruction deleted" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VPhysicalVolume* XrayFluoDetectorConstruction::Construct()
{
  DefineMaterials();
  return ConstructApparate();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::DefineMaterials()
{
  //define elements
  G4String name, symbol;             //a=mass of a mole;
  G4double a, z, density;            //z=mean number of protons;  
  G4int  natoms,ncomponents;
  G4double temperature, pressure;  
  G4double fractionmass;

  // Elements Definitions

  //define Niobium
  
  a = 92.906*g/mole;
  G4Element* Nb  = new G4Element(name="Niobium"  ,symbol="Nb" , z= 41., a);

  //define Zirconium
  
  a = 91.22*g/mole;
  G4Element* Zr  = new G4Element(name="Zirconium"  ,symbol="Zr" , z= 40., a);

  //define Yttrium
  
  a = 88.905*g/mole;
  G4Element* Y  = new G4Element(name="Yttrium"  ,symbol="Y" , z= 39., a);

  //define Stronzium
  
  a = 87.62*g/mole;
  G4Element* Sr  = new G4Element(name="Stronzium"  ,symbol="Sr" , z= 38., a);

  //define Rubidium
  
  a = 85.47*g/mole;
  G4Element* Rb  = new G4Element(name="Rubidium"  ,symbol="Rb" , z= 37., a);

  //define Zinc
  
  a = 65.37*g/mole;
  G4Element* Zn  = new G4Element(name="Zinc"  ,symbol="Zn" , z= 30., a);

  //define Nichel
  
  a = 58.71*g/mole;
  G4Element* Ni  = new G4Element(name="Nichel"  ,symbol="Ni" , z= 28., a);

  //define Scandio
  
  a = 44.956*g/mole;
  G4Element* Sc  = new G4Element(name="Scandium"  ,symbol="Sc" , z= 21., a);

  //define Vanadium
  
  a = 50.942*g/mole;
  G4Element* V  = new G4Element(name="Vanadium"  ,symbol="V" , z= 39., a);

  //define Cromium
 
  a = 51.996*g/mole;
  G4Element* Cr  = new G4Element(name="Cromium"  ,symbol="Cr" , z= 24., a);

  //define Cobalt
  
  a = 58.933*g/mole;
  G4Element* Co  = new G4Element(name="Cobalt"  ,symbol="Co" , z= 27., a);

  //define Copper
  
  a = 63.54*g/mole;
  G4Element* elCu  = new G4Element(name="Copper"  ,symbol="Cu" , z= 29., a);

  //define Barium
  
  a = 137.34*g/mole;
  G4Element* Ba  = new G4Element(name="Barium"  ,symbol="Ba" , z= 56., a);

  //define Cerium
  
  a = 140.12*g/mole;
  G4Element* Ce  = new G4Element(name="Cerium"  ,symbol="Ce" , z= 58., a);

  //define Neodimuim
  
  a = 144.24*g/mole;
  G4Element* Nd  = new G4Element(name="Neodimuim"  ,symbol="Nd" , z= 60., a);


  //Define Zolfo

  a = 32.064*g/mole;
  G4Element* elS  = new G4Element(name="Zolfo"  ,symbol="S" , z= 16., a);



  //define carbon
  
  a = 12.0107*g/mole;
  G4Element* C  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);

  //define Nitrogen

  a = 14.01*g/mole;
  G4Element* N  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

  // define Oxigen
  a = 15.9994*g/mole;
  G4Element* O  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

  //define Arsenic
  a = 74.9216 * g/mole;
  G4Element * As = new G4Element( name="arsenic",symbol="As",z= 33.,a);

  //Define Gallium  
  a = 69.72* g/mole;
  G4Element * Ga = new G4Element(name="gallium",symbol="Ga",z= 31.,a);

  //define Iron  
  a = 55.847*g/mole;
  G4Element* Fe = new G4Element(name="Iron"  ,symbol="Fe", z=26., a);
  //define hydrogen
  
  a = 1.01*g/mole;
  G4Element* H  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);
  
  //define germanium
  a = 72.61*g/mole;
  G4Element* Ge = new G4Element(name="Germanium",symbol="Ge", z= 32.,a);
  //define phosporus

  a = 30.97*g/mole;
  G4Element* P =  new G4Element(name="Phosporus",symbol="P", z= 15., a);

  // define Titanium
  a = 47.88*g/mole;
  G4Element* elTi = new G4Element(name="Titanium",symbol="Ti" , z= 22., a);

  // define Calcium
  a = 40.078*g/mole;
  G4Element* Ca = new G4Element(name="Calcium",symbol="Ca" , z= 20., a);

  // define silicon
  a = 28.0855*g/mole;
  G4Element* elSi = new G4Element(name="Silicon",symbol="Si" , z= 14., a);

  // define Aluminium
  a = 26.98154*g/mole;
  G4Element* elAl = new G4Element(name="Aluminium",symbol="Al" , z= 13., a);

  // Define Magnesium
  a = 24.305*g/mole;
  G4Element* Mg = new G4Element(name="Magnesium",symbol="Mg" , z= 12., a);

  // Define Manganese
  a = 54.938*g/mole;
  G4Element* Mn = new G4Element(name="Manganese",symbol="Mn" , z= 25., a);

  // Define Sodium
  a = 22.989*g/mole;
  G4Element* Na = new G4Element(name="Sodium",symbol="Na" , z= 11., a);

  // Define Potassium
  a = 39.10*g/mole;
  G4Element* K = new G4Element(name="Potassium",symbol="K" , z= 19., a); 
  // Define lead

  a=207.19*g/mole;
  G4Element* elPb = new G4Element(name="Lead",symbol="pb", z=82.,a);

  //define Uranium 

  a =  238.02891*g/mole;
  G4Element* elU  = new G4Element(name="Uranium",symbol="U", z=92.,a);


  G4cout << "Elementi creati" << G4endl;



  // Materials Definitions

  // Define dolorite main components

  density = 3*g/cm3;
  G4Material* diorite = new G4Material(name="Diorite", density, ncomponents=11);
  diorite->AddElement(Fe,   fractionmass=0.1750);
  diorite->AddElement(elTi, fractionmass=0.0082);
  diorite->AddElement(Ca,   fractionmass=0.0753);
  diorite->AddElement(elSi, fractionmass=0.2188);
  diorite->AddElement(elAl, fractionmass=0.0676);
  diorite->AddElement(Mg,   fractionmass=0.0008);
  diorite->AddElement(O ,   fractionmass=0.4377);
  diorite->AddElement(Mn ,  fractionmass=0.0015);
  diorite->AddElement(Na ,  fractionmass=0.0134);
  diorite->AddElement(K ,   fractionmass=0.0011);
  diorite->AddElement(P ,   fractionmass=0.0006);

  // define traces in dolorite

  density = 3*g/cm3;
  G4Material* tracesOfDolorite = new G4Material(name="TracesOfDolorite", density, ncomponents=16);
  tracesOfDolorite->AddElement(Nb,   natoms=5);
  tracesOfDolorite->AddElement(Zr,   natoms=91);
  tracesOfDolorite->AddElement(Y,    natoms=29);
  tracesOfDolorite->AddElement(Sr,   natoms=140);
  tracesOfDolorite->AddElement(Rb,   natoms=3);
  tracesOfDolorite->AddElement(Ga,   natoms=20);
  tracesOfDolorite->AddElement(Zn,   natoms=99);
  tracesOfDolorite->AddElement(Ni,   natoms=77);
  tracesOfDolorite->AddElement(Sc,   natoms=32);
  tracesOfDolorite->AddElement(V,    natoms=314);
  tracesOfDolorite->AddElement(Cr,   natoms=130);
  tracesOfDolorite->AddElement(Co,   natoms=56);
  tracesOfDolorite->AddElement(elCu,   natoms=119);
  tracesOfDolorite->AddElement(Ba,   natoms=38);
  tracesOfDolorite->AddElement(Ce,   natoms=15);
  tracesOfDolorite->AddElement(Nd,   natoms=9);

  // define dolorite (full)

  density = 3*g/cm3;
  dolorite = new G4Material(name="Dolorite", density, ncomponents=2);
  dolorite->AddMaterial(tracesOfDolorite, fractionmass=0.0027842352);
  dolorite->AddMaterial(diorite, fractionmass=0.9972157648);

  // define mars1

  density = 3*g/cm3;
  G4Material* mars1Main = new G4Material(name="Mars1 Main components", density, ncomponents=11);
  mars1Main->AddElement(Fe,   fractionmass=0.100916);
  mars1Main->AddElement(elTi, fractionmass=0.0186804);
  mars1Main->AddElement(Ca,   fractionmass=0.0404091);
  mars1Main->AddElement(elSi, fractionmass=0.196378);
  mars1Main->AddElement(elAl, fractionmass=0.103282);
  mars1Main->AddElement(Mg,   fractionmass=0.0241622);
  mars1Main->AddElement(Mn ,  fractionmass=0.00184331);
  mars1Main->AddElement(Na ,  fractionmass=0.0177908);
  mars1Main->AddElement(K ,   fractionmass=0.00574498);
  mars1Main->AddElement(P ,   fractionmass=0.00280169);
  mars1Main->AddElement(O ,   fractionmass=0.48799152);

  density = 3*g/cm3;
  G4Material* tracesOfMars1 = new G4Material(name="TracesOfMars1", density, ncomponents=17);
  tracesOfMars1->AddElement(Nb,   natoms=55);
  tracesOfMars1->AddElement(Zr,   natoms=433);
  tracesOfMars1->AddElement(Y,    natoms=58);
  tracesOfMars1->AddElement(Sr,   natoms=968);
  tracesOfMars1->AddElement(Rb,   natoms=16);
  tracesOfMars1->AddElement(Ga,   natoms=24);
  tracesOfMars1->AddElement(Zn,   natoms=109);
  tracesOfMars1->AddElement(Ni,   natoms=70);
  tracesOfMars1->AddElement(Sc,   natoms=21);
  tracesOfMars1->AddElement(V,    natoms=134);
  tracesOfMars1->AddElement(Cr,   natoms=141);
  tracesOfMars1->AddElement(Co,   natoms=30);
  tracesOfMars1->AddElement(elCu, natoms=19);
  tracesOfMars1->AddElement(Ba,   natoms=580);
  tracesOfMars1->AddElement(elPb,   natoms=4);
  tracesOfMars1->AddElement(elS,  natoms=444);
  tracesOfMars1->AddElement(elU,    natoms=2);


  density = 3*g/cm3;
  G4Material* mars1 = new G4Material(name="Mars1", density, ncomponents=2);
  mars1->AddMaterial(tracesOfMars1, fractionmass=0.0044963163);
  mars1->AddMaterial(mars1Main, fractionmass=0.9955036837);

  //define Neodimuim
  
  density =  6800*kg/m3;
              materialNd  = new G4Material(name="Neodimuim"  ,density , ncomponents=1);
  materialNd ->AddElement(Nd,natoms=1);

  // Define Magnesium
  density = 1738*kg/m3;
              materialMg = new G4Material(name="Magnesium",density , ncomponents=1);
  materialMg->AddElement(Mg,natoms=1);

  //define iron 

  density = 7.86 * g/cm3;
              FeMaterial = new G4Material(name="Iron",density,ncomponents=1);
  FeMaterial->AddElement(Fe,natoms=1);



  //define gallium arsenide
  
  density = 5.32 * g/cm3;
  G4Material * GaAs = new G4Material(name ="gallium arsenide",density,ncomponents=2);
  GaAs->AddElement(Ga,natoms=1);
  GaAs->AddElement(As,natoms=1);


  // define germanium
  
  density = 5.32 * g/cm3;
              HPGe = new G4Material(name="HPGe",density,ncomponents=1);
  HPGe ->AddElement(Ge,natoms=1);
  
  //define silicon
  
  density = 2.333*g/cm3;
  a = 28.0855*g/mole;
	      Si = new G4Material(name="Silicon",z=14., a,density);
  
  //define copper
  
  density = 8.960*g/cm3;
  a = 63.55*g/mole;
              Cu = new G4Material(name="Copper"   , z=29., a, density);


  
  //define Oxigen
  density = 1*g/cm3;
  a=16*g/mole;
  G4Material* matOx = new G4Material(name="Oxigen", z=8., a, density);
  
  //define aluminium
  
  density = 2.700*g/cm3;
  a = 26.98*g/mole;
	      Al = new G4Material(name="Aluminium", z=13., a, density);

  //define titanium 
  density = 4.54 *g/cm3;
  a = 47.867*g/mole;
              Ti  = new G4Material(name="Titanium",z=22.,a,density);


  //define Uranium 
  density = 19050*kg/m3;
  a =  238.02891*g/mole;
  G4Material* U  = new G4Material(name="Uranium",z=92.,a,density);

  //define Tin
  density = 7310*kg/m3;
  a =  118.710*g/mole;
              Sn  = new G4Material(name="Tin",z=50.,a,density);

  //define lead
  
  density = 11.35*g/cm3;
  a=207.19*g/mole;
  G4Material* Pb = new G4Material(name="Lead",z=82.,a,density);

  //define scintillator
  
  density = 1.032*g/cm3;
  G4Material* Sci = new G4Material(name="Scintillator", density, ncomponents=2);
  Sci->AddElement(C, natoms=9);
  Sci->AddElement(H, natoms=10);  

  
  //define air


  density = 1.290*mg/cm3;
  G4Material* Air = new G4Material(name="Air"  , density, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);
  
  //define vacuum
  
  density     = universe_mean_density;    //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  G4Material * Vacuum = new G4Material(name="Galactic", z=1., a=1.01*g/mole, density,
				       kStateGas,temperature,pressure);




  
  //define basalt
  density = 3.*g/cm3; 
  G4Material* basalt = new G4Material(name="Basalt", density, ncomponents=7);
  basalt->AddElement(Fe, fractionmass=0.1200);
  basalt->AddElement(elTi, fractionmass=0.0160);
  basalt->AddElement(Ca, fractionmass=0.0750);
  basalt->AddElement(elSi, fractionmass=0.2160);
  basalt->AddElement(elAl, fractionmass=0.0710);
  basalt->AddElement(Mg, fractionmass=0.0590);
  basalt->AddElement(O , fractionmass=0.4430);

  // end basalt
  
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  
  //default materials of the apparate
  
  sampleMaterial = mars1;
  Dia1Material = Pb;
  Dia3Material = Pb;
  pixelMaterial = Si;
  OhmicPosMaterial = Cu;
  OhmicNegMaterial = Pb;
  defaultMaterial = Vacuum;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

G4VPhysicalVolume* XrayFluoDetectorConstruction::ConstructApparate()
{
  // complete the apparate parameters definition 
  
  ComputeApparateParameters();
  
  //world
  
  solidWorld = new G4Box("World",	      		        //its name
			 WorldSizeXY/2,WorldSizeXY/2,WorldSizeZ/2);	//its size
  
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
 
  //HPGeDetector
  
  solidHPGe = 0;  physiHPGe = 0;  logicHPGe=0;
  solidPixel=0; logicPixel=0; physiPixel=0;
  
  if (DeviceThickness > 0.)  
    {
      solidHPGe = new G4Box("HPGeDetector",		//its name
			    DeviceSizeX/2,DeviceSizeY/2,DeviceThickness/2);//size
      
      
      logicHPGe = new G4LogicalVolume(solidHPGe,	//its solid
				      defaultMaterial,	//its material
				      "HPGeDetector");	//its name
      
      zRotPhiHPGe.rotateX(PhiHPGe);
      G4double x,y,z;
      z = DistDe * cos(ThetaHPGe);
      y =DistDe * sin(ThetaHPGe);
      x = 0.*cm;
      physiHPGe = new G4PVPlacement(G4Transform3D(zRotPhiHPGe,G4ThreeVector(x,y,z)),                                           "HPGeDetector",	//its name
				    logicHPGe,	//its logical volume
				    physiWorld,	//its mother  volume
				    false,		//no boolean operation
				    0);		//copy number
    }
  // Pixel   
  
  for ( G4int j=0; j < NbOfPixelColumns ; j++ )
    { for ( G4int i=0; i < NbOfPixelRows ; i++ )
      { 
	solidPixel=0; logicPixel=0;   physiPixel=0;
	if (PixelThickness > 0.)
	  solidPixel = new G4Box("Pixel",			
				 PixelSizeXY/2,PixelSizeXY/2, PixelThickness/2);
	
	logicPixel = new G4LogicalVolume(solidPixel,	
					 pixelMaterial,	//its material
					 "Pixel");	        //its name 
	zRotPhiHPGe.rotateX(PhiHPGe);
	G4double x,y,z;
	z = DistDe * cos(ThetaHPGe);
	y =DistDe * sin(ThetaHPGe);
	x = 0.*cm; 
	
	physiPixel = new G4PVPlacement(0,	       
				       G4ThreeVector(0,
						     i*PixelSizeXY,
						     j*PixelSizeXY ),
				       "Pixel",  
				       logicPixel,	 //its logical volume
				       physiHPGe, //its mother  volume
				       false,	 //no boolean operation
				       PixelCopyNb);//copy number
	
	// OhmicNeg
	
	solidOhmicNeg=0; logicOhmicNeg=0; physiOhmicNeg=0;  
	
	if (OhmicNegThickness > 0.) 
	  { solidOhmicNeg = new G4Box("OhmicNeg",		//its name
				      PixelSizeXY/2,PixelSizeXY/2,OhmicNegThickness/2); 
		
	  logicOhmicNeg = new G4LogicalVolume(solidOhmicNeg,    //its solid
						    OhmicNegMaterial, //its material
					      "OhmicNeg");      //its name
		
	  physiOhmicNeg = new G4PVPlacement(0,
					    G4ThreeVector
					    (0.,
					     0.,
					      (PixelThickness+OhmicNegThickness)/2),
					    "OhmicNeg",        //its name
					    logicOhmicNeg,     //its logical volume
					    physiHPGe,        //its mother
					    false,             //no boulean operat
					    PixelCopyNb);                //copy number
	  
	  }
	// OhmicPos
	solidOhmicPos=0; logicOhmicPos=0; physiOhmicPos=0;  
	
	if (OhmicPosThickness > 0.) 
	  { solidOhmicPos = new G4Box("OhmicPos",		//its name
				      ContactSizeXY/2,ContactSizeXY/2,OhmicPosThickness/2); 
	  
	  logicOhmicPos = new G4LogicalVolume(solidOhmicPos,    //its solid
					      OhmicPosMaterial, //its material
					      "OhmicPos");      //its name
	  
	  physiOhmicPos = new G4PVPlacement(0,	
					    G4ThreeVector(0.,
							  0.,
							  (-PixelThickness-OhmicPosThickness)/2),  
					    "OhmicPos",  
					    logicOhmicPos,
					    physiHPGe,  
					    false,     
					    PixelCopyNb); 
	  
	  }
	PixelCopyNb += PixelCopyNb;
      }
    } 
  
    //Sample
    
    solidSample=0;  logicSample=0;  physiSample=0;
    
    if (SampleThickness > 0.)  
      {
	solidSample = new G4Box("Sample",		//its name
				SampleSizeXY/2,SampleSizeXY/2,SampleThickness/2);//size
	
	logicSample= new G4LogicalVolume(solidSample,	//its solid
					 sampleMaterial,	//its material
      				       "Sample");	//its name
	
	physiSample = new G4PVPlacement(0,			//no rotation
					G4ThreeVector(),	//at (0,0,0)
					"Sample",	//its name
					logicSample,	//its logical volume
				      physiWorld,	//its mother  volume
					false,		//no boolean operation
					0);		//copy number
	
      }
    
    //Diaphragm1
    
  solidDia1 = 0;  physiDia1 = 0;  logicDia1=0;
  
  if (Dia1Thickness > 0.)  
    {
      solidDia1 = new G4Tubs("Diaphragm1",		//its name
			     DiaInnerSize/2,
			     Dia1SizeXY/2,
			     Dia1Thickness/2,
			     0,
			     360);//size
      
   
      logicDia1 = new G4LogicalVolume(solidDia1,	//its solid
				      Dia1Material,	//its material
				      "Diaphragm1");	//its name
      
      zRotPhiDia1.rotateX(AlphaDia1);
      G4double x,y,z;
      z = DistDia * cos(ThetaDia1);
      y =DistDia * sin(ThetaDia1);
      x = 0.*cm;
      physiDia1 = new G4PVPlacement(G4Transform3D(zRotPhiDia1,G4ThreeVector(x,y,z)),                                           "Diaphragm1",	//its name
				    logicDia1,	//its logical volume
				    physiWorld,	//its mother  volume
				    false,		//no boolean operation
                                   0);		//copy number
    }  
 
  //Diaphragm3
  
  solidDia3 = 0;  physiDia3 = 0;  logicDia3 =0;
  
  if (Dia3Thickness > 0.)  
    {
      solidDia3 = new G4Tubs("Diaphragm3",
			     Dia3InnerSize/2,
			     Dia3SizeXY/2,
			     Dia3Thickness/2,
			     0,
			     360);
      
      
      logicDia3 = new G4LogicalVolume(solidDia3,	//its solid
				      Dia3Material,	//its material
				      "Diaphragm3");	//its name
      
      zRotPhiDia3.rotateX(AlphaDia3);
      G4double x,y,z;
      z = Dia3Dist * cos(ThetaDia3);
      y =Dia3Dist * sin(ThetaDia3);
      x = 0.*cm;
      physiDia3 = new G4PVPlacement(G4Transform3D(zRotPhiDia3,G4ThreeVector(x,y,z)),                                           "Diaphragm3",	//its name
				    logicDia3,	//its logical volume
				    physiWorld,	//its mother  volume
				    false,		//no boolean operation
				    0);		//copy number
    }    
    
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  
  if(!HPGeSD)
    {
      HPGeSD = new XrayFluoHPGeSD ("HPGeSD",this);
      SDman->AddNewDetector(HPGeSD);
    }
  
  
  if (logicPixel)
    {
      logicPixel->SetSensitiveDetector(HPGeSD);
    }
  
  // Visualization attributes
  
  logicWorld->SetVisAttributes (G4VisAttributes::Invisible);
   G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
   G4VisAttributes * yellow= new G4VisAttributes( G4Colour(255/255. ,255/255. ,51/255. ));
  yellow->SetVisibility(true);
  yellow->SetForceSolid(true);
  simpleBoxVisAtt->SetVisibility(true);
 
  logicPixel->SetVisAttributes(simpleBoxVisAtt);
  logicHPGe->SetVisAttributes(G4VisAttributes::Invisible );
  logicSample->SetVisAttributes(simpleBoxVisAtt);
  
  logicDia1->SetVisAttributes(simpleBoxVisAtt);
  logicDia3->SetVisAttributes(simpleBoxVisAtt);
  
  logicOhmicNeg->SetVisAttributes(yellow);
  logicOhmicPos->SetVisAttributes(yellow);
  //always return the physical World
  
  
  PrintApparateParameters();
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::PrintApparateParameters()
{
  G4cout << "-----------------------------------------------------------------------"
	 << G4endl
	 << "The sample is a box whose size is: "
	 << G4endl      
	 << SampleThickness/cm
	 << " cm * "
	 << SampleSizeXY/cm
	 << " cm * "
	 << SampleSizeXY/cm
	 << " cm"
	 << G4endl
	 <<" Material: " << sampleMaterial->GetName() 
	 <<G4endl
	  <<"The HPGeDetector is a slice  " << DeviceThickness/(1.e-6*m) <<  " micron thick"
	 <<G4endl
	 

<<"-------------------------------------------------------------------------"
	 << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoDetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructApparate());
}

void XrayFluoDetectorConstruction::SetSampleMaterial(G4String newMaterial)
{

  if( newMaterial == "dolorite") {
    sampleMaterial = dolorite;
  }

  else if ( newMaterial == "iron") {
    sampleMaterial = FeMaterial;
}

  else if ( newMaterial == "silicon") {
    sampleMaterial = Si;
}

  else if ( newMaterial == "aluminium") {
    sampleMaterial = Al;
}

  else if ( newMaterial == "copper") {
    sampleMaterial = Cu;
}
  else if ( newMaterial == "germanium") {
    sampleMaterial = HPGe;
}
  else if ( newMaterial == "magnesium") {
    sampleMaterial = materialMg;
}

  else if ( newMaterial == "neodimium") {
    sampleMaterial = materialNd;
}
  else if ( newMaterial == "tin") {
    sampleMaterial = Sn;
}
  else if ( newMaterial == "titanium") {
    sampleMaterial = Ti;
}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....











