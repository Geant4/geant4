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
// $Id: Em8DetectorConstruction.cc,v 1.17 2004/12/03 10:51:12 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 

#include "Em8DetectorConstruction.hh"
#include "Em8DetectorMessenger.hh"
#include "Em8CalorimeterSD.hh"


#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4GeometryManager.hh"
#include "G4RunManager.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4ProductionCuts.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UnitsTable.hh"
#include "G4ios.hh"



const G4double Em8DetectorConstruction::fDelta = 0.0001*mm;




/////////////////////////////////////////////////////////////////////////////
//
//

Em8DetectorConstruction::Em8DetectorConstruction()
:
fWorldChanged(false),
fWorldMaterial(NULL),fSolidWorld(NULL),fLogicWorld(NULL),fPhysicsWorld(NULL),
fAbsorberMaterial(NULL),fSolidAbsorber(NULL),fLogicAbsorber(NULL),
fPhysicsAbsorber(NULL),fDetectorMessenger(NULL),
fCalorimeterSD(NULL),fRegGasDet(NULL)
{
  // default parameter values of the calorimeter

  //  G4double inch = 2.54*cm ;
  // G4double  mil = inch/1000.0 ;

  // G4double delta = 0.0001*mm;

  // fAbsorberThickness = 85.*mm;
  fAbsorberThickness = 23.0*mm;

  fAbsorberRadius    = 10.*cm;
  fAbsorberZ         = 0.*cm ;

  fWindowThick       = 51.0*micrometer ;

  fGammaCut    = 23*mm; 
  fElectronCut = 23*mm; 
  fPositronCut = 23*mm; 

  fDetectorMessenger = new Em8DetectorMessenger(this);
}

//////////////////////////////////////////////////////////////////////////
//
//

Em8DetectorConstruction::~Em8DetectorConstruction()
{ 
  delete fDetectorMessenger;
}

//////////////////////////////////////////////////////////////////////////
//
//

G4VPhysicalVolume* Em8DetectorConstruction::Construct()
{
  DefineMaterials();
  return ConstructCalorimeter();
}

//////////////////////////////////////////////////////////////////////////////
//
//

void Em8DetectorConstruction::DefineMaterials()
{ 
 //This function illustrates the possible ways to define materials
 
G4String name, symbol ;             //a=mass of a mole;
G4double a, z, density ;            //z=mean number of protons;  
// G4int iz, n; 
G4int nel ;                       //iz=number of protons  in an isotope; 
                                   // n=number of nucleons in an isotope;

 G4int ncomponents; 
  // G4int natoms;
// G4double abundance; 
G4double fractionmass;
//G4double temperature, pressure;

//
// define Elements
//

  a = 1.01*g/mole;
  G4Element* elH  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);

  a = 12.01*g/mole;
  G4Element* elC = new G4Element(name="Carbon", symbol="C", z=6., a);

  a = 14.01*g/mole;
  G4Element* elN  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

  a = 16.00*g/mole;
  G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

  a = 39.948*g/mole;
  G4Element* elAr = new G4Element(name="Argon", symbol="Ar", z=18., a);

  //  a = 131.29*g/mole;
  //  G4Element* elXe = new G4Element(name="Xenon", symbol="Xe", z=54., a);
  
  //  a = 19.00*g/mole;
  //  G4Element* elF  = new G4Element(name="Fluorine", symbol="F", z=9., a);


//
// define simple materials
//

     /* ******************************************************************



density = 1.390*g/cm3;
a = 39.95*g/mole;
G4Material* lAr = new G4Material(name="liquidArgon", z=18., a, density);

density = 7.870*g/cm3;
a = 55.85*g/mole;
G4Material* Fe = new G4Material(name="Iron"   , z=26., a, density);

density = 8.960*g/cm3;
a = 63.55*g/mole;
G4Material* Cu = new G4Material(name="Copper"   , z=29., a, density);

density = 19.32*g/cm3;
a =196.97*g/mole;
G4Material* Au = new G4Material(name="Gold"   , z=79., a, density);

density = 11.35*g/cm3;
a = 207.19*g/mole;
G4Material* Pb = new G4Material(name="Lead"     , z=82., a, density);

//
// define a material from elements.   case 1: chemical molecule
//

density = 1.000*g/cm3;
G4Material* H2O = new G4Material(name="Water", density, ncomponents=2);
H2O->AddElement(elH, natoms=2);
H2O->AddElement(elO, natoms=1);

  // Kapton (polyimide) ??? since = Mylar C5H4O2

  density = 1.39*g/cm3;
  G4Material* Kapton = new G4Material(name="Kapton", density, nel=3);
  Kapton->AddElement(elO,2);
  Kapton->AddElement(elC,5);
  Kapton->AddElement(elH,4);


  // Carbon dioxide

  density = 1.977*mg/cm3;
  G4Material* CO2 = new G4Material(name="CO2", density, nel=2,
				       kStateGas,273.15*kelvin,1.*atmosphere);
  CO2->AddElement(elC,1);
  CO2->AddElement(elO,2);


  // TRT_CH2
      
  density = 0.935*g/cm3;
  G4Material* TRT_CH2 = new G4Material(name="TRT_CH2",density, nel=2);
  TRT_CH2->AddElement(elC,1);
  TRT_CH2->AddElement(elH,2);

  // Radiator

  density = 0.059*g/cm3;
  G4Material* Radiator = new G4Material(name="Radiator",density, nel=2);
  Radiator->AddElement(elC,1);
  Radiator->AddElement(elH,2);

  // Carbon Fiber

  density = 0.145*g/cm3;
  G4Material* CarbonFiber = new G4Material(name="CarbonFiber",density, nel=1);
  CarbonFiber->AddElement(elC,1);

  density = 1.290*mg/cm3;  // old air from elements
  G4Material* air = new G4Material(name="air"  , density, ncomponents=2);
  Air->AddElement(elN, fractionmass=0.7);
  Air->AddElement(elO, fractionmass=0.3);


  density = 1.25053*mg/cm3 ;       // STP
  a = 14.01*g/mole ;       // get atomic weight !!!
  //  a = 28.016*g/mole;
  G4Material* N2  = new G4Material(name="Nitrogen", z= 7.,a,density) ;

  density = 1.25053*mg/cm3 ;       // STP
  G4Material* anotherN2 = new G4Material(name="anotherN2", density,ncomponents=2);
  anotherN2->AddElement(elN, 1);
  anotherN2->AddElement(elN, 1);


  // liquid hydrogen for muon cooling target

  density = 0.071*g/cm3;
  a = 1.01*g/mole;
  G4Material* lH2 = new G4Material(name="liquidHydrigen", z=1., a, density);

  // Beryllium

  density = 1.848*g/cm3;
  a = 9.01*g/mole;
  G4Material* Be = new G4Material(name="Beryllium", z=4., a, density);
  // Al for electrodes

  density = 2.700*g/cm3;
  a = 26.98*g/mole;
  G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);


  // Polypropelene

  G4Material* CH2 = new G4Material ("Polypropelene" , 0.91*g/cm3, 2);
  CH2->AddElement(elH,2);
  CH2->AddElement(elC,1);

************************ */
  // Aluminium

  a = 26.98*g/mole;
  density = 2.7*g/cm3;
  G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);
  if(Al);
  // Mylar

  density = 1.39*g/cm3;
  G4Material* Mylar = new G4Material(name="Mylar", density, nel=3);
  Mylar->AddElement(elO,2);
  Mylar->AddElement(elC,5);
  Mylar->AddElement(elH,4);

  // Silicon as detector material

  density = 2.330*g/cm3;
  a = 28.09*g/mole;
  G4Material* Si = new G4Material(name="Silicon", z=14., a, density);
  if(Si);
  // Krypton as detector gas, STP

  density = 3.700*mg/cm3 ;
  a = 83.80*g/mole ;
  G4Material* Kr  = new G4Material(name="Kr",z=36., a, density );

  // Metane, STP

  //  density = 0.7174*mg/cm3 ;
  //  G4Material* metane = new G4Material(name="CH4",density,nel=2) ;
  //  metane->AddElement(elC,1) ;
  //  metane->AddElement(elH,4) ;


  // Dry air (average composition)

  density = 1.7836*mg/cm3 ;       // STP
  G4Material* Argon = new G4Material(name="Argon"  , density, ncomponents=1);
  Argon->AddElement(elAr, 1);

  density = 1.25053*mg/cm3 ;       // STP
  G4Material* Nitrogen = new G4Material(name="N2"  , density, ncomponents=1);
  Nitrogen->AddElement(elN, 2);

  density = 1.4289*mg/cm3 ;       // STP
  G4Material* Oxygen = new G4Material(name="O2"  , density, ncomponents=1);
  Oxygen->AddElement(elO, 2);


  density = 1.2928*mg/cm3 ;       // STP
  G4Material* Air = new G4Material(name="Air"  , density, ncomponents=3);
  Air->AddMaterial( Nitrogen, fractionmass = 0.7557 ) ;
  Air->AddMaterial( Oxygen,   fractionmass = 0.2315 ) ;
  Air->AddMaterial( Argon,    fractionmass = 0.0128 ) ;


  // 93% Kr + 7% CH4, STP

  //  density = 3.491*mg/cm3 ;      
  //  G4Material* Kr7CH4 = new G4Material(name="Kr7CH4"  , density, 
  //                                                   ncomponents=2);
  //  Kr7CH4->AddMaterial( Kr,       fractionmass = 0.986 ) ;
  //  Kr7CH4->AddMaterial( metane,   fractionmass = 0.014 ) ;

  /* **************

  G4double TRT_Xe_density = 5.485*mg/cm3;
  G4Material* TRT_Xe = new G4Material(name="TRT_Xe", TRT_Xe_density, nel=1,
				      kStateGas,293.15*kelvin,1.*atmosphere);
  TRT_Xe->AddElement(elXe,1);

  G4double TRT_CO2_density = 1.842*mg/cm3;
  G4Material* TRT_CO2 = new G4Material(name="TRT_CO2", TRT_CO2_density, nel=2,
				       kStateGas,293.15*kelvin,1.*atmosphere);
  TRT_CO2->AddElement(elC,1);
  TRT_CO2->AddElement(elO,2);

  G4double TRT_CF4_density = 3.9*mg/cm3;
  G4Material* TRT_CF4 = new G4Material(name="TRT_CF4", TRT_CF4_density, nel=2,
                                       kStateGas,293.15*kelvin,1.*atmosphere);
  TRT_CF4->AddElement(elC,1);
  TRT_CF4->AddElement(elF,4);

  // ATLAS TRT straw tube gas mixture (20 C, 1 atm)

  G4double XeCO2CF4_density = 4.76*mg/cm3;
  G4Material* XeCO2CF4 = new G4Material(name="XeCO2CF4", XeCO2CF4_density,
					ncomponents=3,
					kStateGas,293.15*kelvin,1.*atmosphere);
  XeCO2CF4->AddMaterial(TRT_Xe,0.807);
  XeCO2CF4->AddMaterial(TRT_CO2,0.039);
  XeCO2CF4->AddMaterial(TRT_CF4,0.154);

  *********** */

  // Xenon as detector gas, STP

  density = 5.858*mg/cm3 ;
  a = 131.29*g/mole ;
  G4Material* Xe  = new G4Material(name="Xenon",z=54., a, density );

  // Metane, STP

  density = 0.7174*mg/cm3 ;
  G4Material* metane = new G4Material(name="CH4",density,nel=2) ;
  metane->AddElement(elC,1) ;
  metane->AddElement(elH,4) ;

  // C3H8,20 C, 2 atm

  density = 3.758*mg/cm3 ;
  G4Material* C3H8 = new G4Material(name="C3H8",density,nel=2) ;
  C3H8->AddElement(elC,3) ;
  C3H8->AddElement(elH,8) ;

  // Propane, STP

  density = 2.005*mg/cm3 ;
  G4Material* propane = new G4Material(name="propane",density,nel=2) ;
  propane->AddElement(elC,3) ;
  propane->AddElement(elH,8) ;

  // 87.5% Xe + 7.5% CH4 + 5% C3H8, 20 C, 1 atm 

  density = 4.9196*mg/cm3 ;

  G4Material* XeCH4C3H8 = new G4Material(name="XeCH4C3H8"  , 
                                  density,  ncomponents=3);
  XeCH4C3H8->AddMaterial( Xe,       fractionmass = 0.971 ) ;
  XeCH4C3H8->AddMaterial( metane,   fractionmass = 0.010 ) ;
  XeCH4C3H8->AddMaterial( propane,  fractionmass = 0.019 ) ;

  // 93% Ar + 7% CH4, STP

  density = 1.709*mg/cm3 ;      
  G4Material* Ar7CH4 = new G4Material(name="Ar7CH4", density, ncomponents=2);
  Ar7CH4->AddMaterial( Argon,    fractionmass = 0.971 ) ;
  Ar7CH4->AddMaterial( metane,   fractionmass = 0.029 ) ;

  // Carbon dioxide, STP

  density = 1.977*mg/cm3;
  G4Material* CarbonDioxide = new G4Material(name="CO2", density, nel=2);
  CarbonDioxide->AddElement(elC,1);
  CarbonDioxide->AddElement(elO,2);

  // 80% Ar + 20% CO2, STP

//  density = 1.8223*mg/cm3 ;      
//  G4Material* Ar_80CO2_20 = new G4Material(name="ArCO2"  , density, 
//                                    ncomponents=2);
//  Ar_80CO2_20->AddMaterial( Argon,           fractionmass = 0.783 ) ;
//  Ar_80CO2_20->AddMaterial( CarbonDioxide,   fractionmass = 0.217 ) ;

  // 80% Xe + 20% CO2, STP

  //  density = 5.0818*mg/cm3 ;      
  //  G4Material* Xe20CO2 = new G4Material(name="Xe20CO2"  , 
  // density, ncomponents=2);
  //  Xe20CO2->AddMaterial( Xe,              fractionmass = 0.922 ) ;
  //  Xe20CO2->AddMaterial( CarbonDioxide,   fractionmass = 0.078 ) ;

  // 80% Kr + 20% CO2, STP

  density = 3.601*mg/cm3 ;      
  G4Material* Kr20CO2 = new G4Material(name="Kr20CO2"  , density, 
                                        ncomponents=2);
  Kr20CO2->AddMaterial( Kr,              fractionmass = 0.89 ) ;
  Kr20CO2->AddMaterial( CarbonDioxide,   fractionmass = 0.11 ) ;

  // G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  
  // fWindowMat = Mylar ;
 
  fAbsorberMaterial = XeCH4C3H8;
  // Al; // Si; // Xe; // Ar7CH4; // C3H8; // XeCH4C3H8; 

  fWorldMaterial    = Mylar; // Air ;
}

/////////////////////////////////////////////////////////////////////////
//
// 
  
G4VPhysicalVolume* Em8DetectorConstruction::ConstructCalorimeter()
{
  // Cleanup old geometry

  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  //  G4RegionStore::GetInstance()->Clean();

  // complete the Calor parameters definition and print 

  ComputeCalorParameters();
  PrintCalorParameters();
      
  // World
  
  // if(fSolidWorld)   delete fSolidWorld ;
  // if(fLogicWorld)   delete fLogicWorld ;
  // if(fPhysicsWorld) delete fPhysicsWorld ;

  fSolidWorld = new G4Tubs("World",				//its name
                   0.,fWorldSizeR,fWorldSizeZ/2.,0.,twopi)       ;//its size
                         
  fLogicWorld = new G4LogicalVolume(fSolidWorld,		//its solid
                                   fWorldMaterial,	//its material
                                   "World");		//its name
                                   
  fPhysicsWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 "World",		//its name
                                 fLogicWorld,		//its logical volume
                                 NULL,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number

  // Absorber

  if (fAbsorberThickness > 0.) 
  { 
    // if(fSolidAbsorber)   delete fSolidAbsorber ;
    // if(fLogicAbsorber)   delete fLogicAbsorber ;
    // if(fPhysicsAbsorber) delete fPhysicsAbsorber ;

      fSolidAbsorber = new G4Tubs("Absorber",		
                          0.,fAbsorberRadius,fAbsorberThickness/2.,0.,twopi); 
                          
      fLogicAbsorber = new G4LogicalVolume(fSolidAbsorber,    
      			                  fAbsorberMaterial, 
      			                  "Absorber");     
      			                  
      fPhysicsAbsorber = new G4PVPlacement(0,		   
      		    G4ThreeVector(0.,0.,fAbsorberZ),        
                                        "Absorber",        
                                        fLogicAbsorber,     
                                        fPhysicsWorld,       
                                        false,             
                                        0);                
                                        
  }
  if( fRegGasDet != 0 )  // remove obsolete root logical volume
  {
    fRegGasDet->RemoveRootLogicalVolume(fLogicAbsorber);
  }
  G4ProductionCuts* cuts = 0;

  if( fRegGasDet == 0 ) // First time - instantiate a region and a cut objects
  {    
    fRegGasDet = new G4Region("VertexDetector");
    cuts = new G4ProductionCuts();
    fRegGasDet->SetProductionCuts(cuts);
  }
  else  // Second time - get a cut object from region
  {   
    cuts = fRegGasDet->GetProductionCuts();
  }
  fRegGasDet->AddRootLogicalVolume(fLogicAbsorber);                               

  cuts->SetProductionCut(fGammaCut,"gamma");
  cuts->SetProductionCut(fElectronCut,"e-");
  cuts->SetProductionCut(fPositronCut,"e+");

  // Sensitive Detectors: Absorber 
  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  if(!fCalorimeterSD)
  {
    fCalorimeterSD = new Em8CalorimeterSD("CalorSD",this);
    SDman->AddNewDetector( fCalorimeterSD );
  }
  if (fLogicAbsorber)  fLogicAbsorber->SetSensitiveDetector(fCalorimeterSD);

  // Parameterisation

  //   G4VXrayTRmodel* pTRModel = new G4IrregularXrayTRmodel(logicRadiator,
  //        fRadThickness,fGasGap);

  //  G4VXrayTRmodel* pTRModel = new G4FoamXrayTRmodel(logicRadiator,
  //                                    fRadThickness,fGasGap);

  //   G4VXrayTRmodel* pTRModel = new G4RegularXrayTRmodel(logicRadiator,
  //                                                       fRadThickness,fGasGap);

  //  G4double alphaPlate = 160.0 ;
  //  G4double alphaGas   = 160.0 ;

  //  G4VXrayTRmodel* pTRModel = new G4GamDistrXrayTRmodel(logicRadiator,
  //					       fRadThickness,alphaPlate,
  //                                              fGasGap,alphaGas);

  //  G4VXrayTRmodel* pTRModel = new G4PlateIrrGasXrayTRmodel(logicRadiator,
  //       fRadThickness,fGasGap);

  //  pTRModel->GetPlateZmuProduct() ;
  //  pTRModel->GetGasZmuProduct() ;

  //  pTRModel->GetNumberOfPhotons() ;  
  
  // always return physics world

  return fPhysicsWorld;
}

////////////////////////////////////////////////////////////////////////////
//
//

void Em8DetectorConstruction::PrintCalorParameters()
{
  G4cout << "\n The  WORLD   is made of " 
       << fWorldSizeZ/mm << "mm of " << fWorldMaterial->GetName() ;
  G4cout << ", the transverse size (R) of the world is " << fWorldSizeR/mm << " mm. " << G4endl;
  G4cout << " The ABSORBER is made of " 
       << fAbsorberThickness/mm << "mm of " << fAbsorberMaterial->GetName() ;
  G4cout << ", the transverse size (R) is " << fAbsorberRadius/mm << " mm. " << G4endl;
  G4cout << " Z position of the (middle of the) absorber " << fAbsorberZ/mm << "  mm." << G4endl;
  G4cout << G4endl;
}

///////////////////////////////////////////////////////////////////////////
//
//

void Em8DetectorConstruction::SetAbsorberMaterial(G4String materialChoice)
{
  // get the pointer to the material table
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  // search the material by its name   
  G4Material* pttoMaterial;

  for (size_t J = 0 ; J < theMaterialTable->size() ; J++)
  { 
    pttoMaterial = (*theMaterialTable)[J];
     
    if(pttoMaterial->GetName() == materialChoice)
    {
      fAbsorberMaterial = pttoMaterial;
      fLogicAbsorber->SetMaterial(pttoMaterial); 

        // PrintCalorParameters();
    }             
  }
}

////////////////////////////////////////////////////////////////////////////
//
//

void Em8DetectorConstruction::SetWorldMaterial(G4String materialChoice)
{
  // get the pointer to the material table
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  // search the material by its name   
  G4Material* pttoMaterial;
  for (size_t J=0 ; J<theMaterialTable->size() ; J++)
  { 
    pttoMaterial = (*theMaterialTable)[J]; 
    
    if(pttoMaterial->GetName() == materialChoice)
    {
      fWorldMaterial = pttoMaterial;
      fLogicWorld->SetMaterial(pttoMaterial);
 
       //  PrintCalorParameters();
    }             
  }
}

///////////////////////////////////////////////////////////////////////////
//
//

void Em8DetectorConstruction::SetAbsorberThickness(G4double val)
{
  // change Absorber thickness and recompute the calorimeter parameters
  fAbsorberThickness = val;
  ComputeCalorParameters();
}  

/////////////////////////////////////////////////////////////////////////////
//
//

void Em8DetectorConstruction::SetAbsorberRadius(G4double val)
{
  // change the transverse size and recompute the calorimeter parameters
  fAbsorberRadius = val;
  ComputeCalorParameters();
}  

////////////////////////////////////////////////////////////////////////////
//
//

void Em8DetectorConstruction::SetWorldSizeZ(G4double val)
{
  fWorldChanged=true;
  fWorldSizeZ = val;
  ComputeCalorParameters();
}  

///////////////////////////////////////////////////////////////////////////
//
//

void Em8DetectorConstruction::SetWorldSizeR(G4double val)
{
  fWorldChanged=true;
  fWorldSizeR = val;
  ComputeCalorParameters();
}  

//////////////////////////////////////////////////////////////////////////////
//
//

void Em8DetectorConstruction::SetAbsorberZpos(G4double val)
{
  fAbsorberZ  = val;
  ComputeCalorParameters();
}  


///////////////////////////////////////////////////////////////////////////////
//
//
  
void Em8DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructCalorimeter());
}

//
//
////////////////////////////////////////////////////////////////////////////

















