//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

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
    solidSi(0),logicSi(0),physiSi(0),
    solidHPGe(0),logicHPGe(0),physiHPGe(0),
    solidSample (0),logicSample(0),physiSample (0),
    solidDia1(0),logicDia1(0),physiDia1(0),
    solidDia2(0),logicDia2(0),physiDia2(0),
    solidDia3(0),logicDia3(0),physiDia3(0),
    solidOhmicPos(0),logicOhmicPos(0), physiOhmicPos(0),
    solidOhmicNeg(0),logicOhmicNeg(0), physiOhmicNeg(0),
    solidPixel(0),logicPixel(0), physiPixel(0),
    OhmicPosMaterial(0), OhmicNegMaterial(0),
    SiMaterial(0),pixelMaterial(0),sampleMaterial(0),
    Dia1Material(0),Dia2Material(0),Dia3Material(0),
    defaultMaterial(0)
  ,HPGeSD(0)
  
{ 
  NbOfPixelRows     =  1;
  NbOfPixelColumns  =  1;
  NbOfPixels        =  NbOfPixelRows*NbOfPixelColumns;
  SiSizeXY = 1. * cm;
  SiThickness = 1.01* mm;
    PixelSizeXY       = 0.7 * cm; 
  //PixelSizeXY       = 2. * cm; 
  PixelThickness =  1. * mm;
  ContactSizeXY     = 0.005*mm;
  //ContactSizeXY = 0.5 * cm;
  SampleThickness = 0.25 * mm;
  //SampleThickness = 0.25 * cm;
  SampleSizeXY = 3. * cm;
  Dia1Thickness = 1. *mm;
  Dia2Thickness = 1. *mm;
  Dia3Thickness = 1. *mm;
  Dia1SizeXY = 3. *cm;
  Dia2SizeXY = 3. *cm;
  Dia3SizeXY = 3. *cm;
  //DiaInnerSize = 1.8 * cm;
  DiaInnerSize = 1.4 * mm;
  OhmicNegThickness = 0.005*mm;
  OhmicPosThickness = 0.005*mm;
  ThetaHPGe = 135. * deg;
  ThetaSi = 210. * deg;
  Dia3Dist =  66.5 * mm;
  Dia3InnerSize = 1. * mm;
  PhiHPGe = 225. * deg;
  PhiSi = 150. * deg;
  ThetaDia1 = 135. * deg;
  ThetaDia2 = 210. * deg;
  ThetaDia3 = 180. * deg;
  PhiDia3 = 90. * deg;
  DistDia = 66.5 * mm;
  //DistDia = 26.5 * mm;
  DistDe =DistDia+ (Dia1Thickness
		    +PixelThickness)/2+OhmicPosThickness ;
 DistSi = DistDia + SiThickness/2;
  PhiDia1 = 90. * deg;
  PhiDia2 = 90. * deg;
  AlphaDia1 = 225. * deg;
  AlphaDia2 = 150. * deg;
  AlphaDia3 = 180. * deg;
  PixelCopyNb=0;
  ComputeApparateParameters();
  
  // create commands for interactive definition of the apparate
  
  detectorMessenger = new XrayFluoDetectorMessenger(this);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoDetectorConstruction::~XrayFluoDetectorConstruction()
{ delete detectorMessenger;}

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
  
 
  a = 74.9216 * g/mole;
  G4Element * As = new G4Element( name="arsenic",symbol="As",z= 33.,a);
  
  a = 69.72* g/mole;
  G4Element * Ga = new G4Element(name="gallium",symbol="Ga",z= 31.,a);
  
  a = 55.85*g/mole;
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
 
  density = 7.86 * g/cm3;
  
  G4int  natoms,ncomponents;
  G4double temperature, pressure;
  
  G4Material * FeMaterial = new G4Material(name="Iron",density,ncomponents=1);
  FeMaterial->AddElement(Fe,natoms=1);
  
  density = 5.32 * g/cm3;
  G4Material * HPGe = new G4Material(name="HPGe",density,ncomponents=1);
  HPGe ->AddElement(Ge,natoms=1);
  //define gallium arsenide
  
  density = 5.32 * g/cm3;
  G4Material * GaAs = new G4Material(name ="gallium arsenide",density,ncomponents=2);
  GaAs->AddElement(Ga,natoms=1);
  GaAs->AddElement(As,natoms=1);
  
  //define silicon
  
  density = 2.333*g/cm3;
  a = 28.09*g/mole;
  G4Material* Si = new G4Material(name="Silicon",z=14., a,density);
  
  //define copper
  
  density = 8.960*g/cm3;
  a = 63.55*g/mole;
  G4Material* Cu = new G4Material(name="Copper"   , z=29., a, density);
  
  
  //define carbon
  
  a = 12.01*g/mole;
  G4Element* C  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);
  
  
  //define scintillator
  
  density = 1.032*g/cm3;
  G4Material* Sci = new G4Material(name="Scintillator", density, ncomponents=2);
  Sci->AddElement(C, natoms=9);
  Sci->AddElement(H, natoms=10);
  
  //define aluminium
  
  density = 2.700*g/cm3;
  a = 26.98*g/mole;
  G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);
  
  //define titanium 
  density = 4.54 *g/cm3;
 a = 47.867*g/mole;
 G4Material* Ti  = new G4Material(name="Titanium",z=22.,a,density);



//define lead
  
  density = 11.35*g/cm3;
  a=207.19*g/mole;
  G4Material* Pb = new G4Material(name="Lead",z=82.,a,density);
  
  //define air

  a = 14.01*g/mole;
  G4Element* N  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);
  
  a = 16.00*g/mole;
  G4Element* O  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);
  
  G4double fractionmass;
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
  
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  
  //default materials of the apparate
  
  SiMaterial = Si;
  sampleMaterial = Ti;
  Dia1Material = Pb;
  Dia2Material = Pb;
  Dia3Material = Pb;
  pixelMaterial = HPGe;
  //pixelMaterial =Al;
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
  /*
  //SiDetector
  
  solidSi = 0;  physiSi = 0;  logicSi=0;
  
  if (SiThickness > 0.)  
    {
      solidSi = new G4Box("SiDetector",		//its name
			  SiSizeXY/2,SiSizeXY/2,SiThickness/2);//size
      
      
      logicSi = new G4LogicalVolume(solidSi,	//its solid
				    SiMaterial,	//its material
				    "SiDetector");	//its name
      
      zRotPhiSi.rotateX(PhiSi);
      G4double x,y,z;
      z = DistSi * cos(ThetaSi);
      y =DistSi * sin(ThetaSi);
      x = 0.*cm;
      physiSi = new G4PVPlacement(G4Transform3D(zRotPhiSi,G4ThreeVector(x,y,z)),                                           "SiDetector",	//its name
				  logicSi,	//its logical volume
				  physiWorld,	//its mother  volume
				  false,		//no boolean operation
				  0);		//copy number
    }
  */
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
  //solidPixel=0; logicPixel=0;   physiPixel=0;
  
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
	
	physiPixel = new G4PVPlacement(
				       // G4Transform3D(zRotPhiHPGe,
				       //    G4ThreeVector(x,y+(i*PixelSizeXY),z+(j*PixelSizeXY))),
				       0,	       
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
  /*   
  //Diaphragm2
  
  solidDia2 = 0;  physiDia2 = 0;  logicDia2=0;
  
  if (Dia2Thickness > 0.)  
    {
      solidDia2 = new G4Tubs("Diaphragm2",
			     DiaInnerSize/2,
			     Dia2SizeXY/2,
			     Dia2Thickness/2,
			     0,
			     360);
      
      
      logicDia2 = new G4LogicalVolume(solidDia2,	//its solid
				      Dia2Material,	//its material
				      "Diaphragm2");	//its name
      
      zRotPhiDia2.rotateX(AlphaDia2);
      G4double x,y,z;
      z = DistDia * cos(ThetaDia2);
      y =DistDia * sin(ThetaDia2);
      x = 0.*cm;
      physiDia2 = new G4PVPlacement(G4Transform3D(zRotPhiDia2,G4ThreeVector(x,y,z)),                                           "Diaphragm2",	//its name
				logicDia2,	//its logical volume
				    physiWorld,	//its mother  volume
				    false,		//no boolean operation
				    0);		//copy number
    } 
  */
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
  //logicSi->SetVisAttributes(simpleBoxVisAtt);
 
  logicPixel->SetVisAttributes(simpleBoxVisAtt);
  logicHPGe->SetVisAttributes(G4VisAttributes::Invisible );
  logicSample->SetVisAttributes(simpleBoxVisAtt);
  
  logicDia1->SetVisAttributes(simpleBoxVisAtt);
  
  // logicDia2->SetVisAttributes(simpleBoxVisAtt);
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....







