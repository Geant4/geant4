// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em10DetectorConstruction.cc,v 1.2 2001-03-19 17:59:03 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "Em10DetectorConstruction.hh"
#include "Em10DetectorMessenger.hh"
#include "Em10CalorimeterSD.hh"

#include "G4VXrayTRmodel.hh"
#include "G4VXrayTRadModel.hh"
#include "G4IrregularXrayTRmodel.hh"
#include "G4FoamXrayTRmodel.hh"
#include "G4RegularXrayTRmodel.hh"
#include "G4GamDistrXrayTRmodel.hh"
#include "G4PlateIrrGasXrayTRmodel.hh"

#include "G4VXTRdEdx.hh"
#include "G4IrregularXTRdEdx.hh"
#include "G4FoamXTRdEdx.hh"
#include "G4RegularXTRdEdx.hh"
#include "G4GamDistrXTRdEdx.hh"
#include "G4PlateIrrGasXTRdEdx.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"

#include "G4ios.hh"

/////////////////////////////////////////////////////////////////////////////
//
//

Em10DetectorConstruction::Em10DetectorConstruction()
:solidWorld(NULL),logicWorld(NULL),physiWorld(NULL),
 solidAbsorber(NULL),logicAbsorber(NULL),physiAbsorber(NULL),
 AbsorberMaterial(NULL),WorldMaterial(NULL),fRadiatorMat(NULL),
 magField(NULL),calorimeterSD(NULL),worldchanged(false),fXTRModel(NULL)
{
  // default parameter values of the calorimeter

  G4double inch = 2.54*cm ;
  G4double  mil = inch/1000.0 ;

  WorldSizeZ = 80.*cm;
  WorldSizeR = 20.*cm;

  // Radiator and detector parameters

  fRadThickness = 40.0*micrometer ; // 25*micrometer ;     
  fGasGap       = 0.126*mm ;          //  1500*micrometer  ;   
  fFoilNumber   = 300 ;             //  188 ;

  AbsorberThickness = 3.0*cm ;   // 40.0*mm ;

  AbsorberRadius   = 10.*cm;
  zAbsorber = 36.*cm ;

  fWindowThick = 51.0*micrometer ;
  fElectrodeThick = 10.0*micrometer ;
  fGapThick = 1.0*mm ;


  fDetThickness = 40.0*mm ;
  fDetLength    = 200.0*cm  ;
  fDetGap       = 1.0*mm ;

  fStartR       = 40*cm  ;
  fStartZ       = 10.0*mm  ;

  fModuleNumber = 1      ;  

  // create commands for interactive definition of the calorimeter  

  detectorMessenger = new Em10DetectorMessenger(this);
}

//////////////////////////////////////////////////////////////////////////
//
//

Em10DetectorConstruction::~Em10DetectorConstruction()
{ 
  delete detectorMessenger;
  if(fXTRModel) delete fXTRModel;
}

//////////////////////////////////////////////////////////////////////////
//
//

G4VPhysicalVolume* Em10DetectorConstruction::Construct()
{
  DefineMaterials();
  return ConstructCalorimeter();  
}

//////////////////////////////////////////////////////////////////////////////
//
//

void Em10DetectorConstruction::DefineMaterials()
{ 
 //This function illustrates the possible ways to define materials
 
G4String name, symbol ;             //a=mass of a mole;
G4double a, z, density ;            //z=mean number of protons;  
G4int iz, n, nel ;                       //iz=number of protons  in an isotope; 
                                   // n=number of nucleons in an isotope;

G4int ncomponents, natoms;
G4double abundance, fractionmass;
G4double temperature, pressure;

//
// define Elements
//

  a = 1.01*g/mole;
  G4Element* elH  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);

  a = 6.94*g/mole;
  G4Element* elLi  = new G4Element(name="Lithium",symbol="Li" , z= 3., a);

  a = 9.01*g/mole;
  G4Element* elBe  = new G4Element(name="Berillium",symbol="Be" , z= 4., a);

  a = 6.01*g/mole;
  G4Element* elC = new G4Element(name="Carbon", symbol="C", z=6., a);

  a = 14.01*g/mole;
  G4Element* elN  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

  a = 16.00*g/mole;
  G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

  a = 39.948*g/mole;
  G4Element* elAr = new G4Element(name="Argon", symbol="Ar", z=18., a);

  a = 131.29*g/mole;
  G4Element* elXe = new G4Element(name="Xenon", symbol="Xe", z=54., a);
  
  a = 19.00*g/mole;
  G4Element* elF  = new G4Element(name="Fluorine", symbol="F", z=9., a);

/////////////////////////////////////////////////////////////////
//
// Detector windows, electrodes 
  // Al for electrodes

  density = 2.700*g/cm3;
  a = 26.98*g/mole;
  G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);


////////////////////////////////////////////////////////////////////////////
//
// Materials for popular X-ray TR radiators
//

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

  // Lithium

  density = 0.534*g/cm3;
  G4Material* Li = new G4Material(name="Li",density, nel=1);
  Li->AddElement(elLi,1);

  // Beryllium

  density = 1.848*g/cm3;
  G4Material* Be = new G4Material(name="Be",density, nel=1);
  Be->AddElement(elBe,1);

  // Mylar

  density = 1.39*g/cm3;
  G4Material* Mylar = new G4Material(name="Mylar", density, nel=3);
  Mylar->AddElement(elO,2);
  Mylar->AddElement(elC,5);
  Mylar->AddElement(elH,4);

  // Kapton (polyimide) ??? since = Mylar C5H4O2

  density = 1.39*g/cm3;
  G4Material* Kapton = new G4Material(name="Kapton", density, nel=3);
  Kapton->AddElement(elO,2);
  Kapton->AddElement(elC,5);
  Kapton->AddElement(elH,4);

  // Polypropelene

  G4Material* CH2 = new G4Material ("Polypropelene" , 0.91*g/cm3, 2);
  CH2->AddElement(elH,2);
  CH2->AddElement(elC,1);

//////////////////////////////////////////////////////////////////////////
//
// Noble gases , STP conditions

  // Helium as detector gas, STP

  density = 0.178*mg/cm3 ;
  a = 4.0026*g/mole ;
  G4Material* He  = new G4Material(name="He",z=2., a, density );

  // Neon as detector gas, STP

  density = 0.900*mg/cm3 ;
  a = 20.179*g/mole ;
  G4Material* Ne  = new G4Material(name="Ne",z=10., a, density );

  // Argon as detector gas, STP

  density = 1.7836*mg/cm3 ;       // STP
  G4Material* Argon = new G4Material(name="Argon"  , density, ncomponents=1);
  Argon->AddElement(elAr, 1);

  // Krypton as detector gas, STP

  density = 3.700*mg/cm3 ;
  a = 83.80*g/mole ;
  G4Material* Kr  = new G4Material(name="Kr",z=36., a, density );

  // Xenon as detector gas, STP

  density = 5.858*mg/cm3 ;
  a = 131.29*g/mole ;
  G4Material* Xe  = new G4Material(name="Xenon",z=54., a, density );

/////////////////////////////////////////////////////////////////////////////
//
// Hydrocarbones, metane and others

  // Metane, STP

  density = 0.7174*mg/cm3 ;
  G4Material* metane = new G4Material(name="CH4",density,nel=2) ;
  metane->AddElement(elC,1) ;
  metane->AddElement(elH,4) ;

  // Propane, STP

  density = 2.005*mg/cm3 ;
  G4Material* propane = new G4Material(name="C3H8",density,nel=2) ;
  propane->AddElement(elC,3) ;
  propane->AddElement(elH,8) ;

  // iso-Butane (methylpropane), STP

  density = 2.67*mg/cm3 ;
  G4Material* isobutane = new G4Material(name="isoC4H10",density,nel=2) ;
  isobutane->AddElement(elC,4) ;
  isobutane->AddElement(elH,10) ;

///////////////////////////////////////////////////////////////////////////
//
// Molecular gases

  // Carbon dioxide

  density = 1.977*mg/cm3;
  G4Material* CO2 = new G4Material(name="CO2", density, nel=2,
				       kStateGas,273.15*kelvin,1.*atmosphere);
  CO2->AddElement(elC,1);
  CO2->AddElement(elO,2);

  // Carbon dioxide, STP

  density = 1.977*mg/cm3;
  G4Material* CarbonDioxide = new G4Material(name="CO2", density, nel=2);
  CarbonDioxide->AddElement(elC,1);
  CarbonDioxide->AddElement(elO,2);


  // Nitrogen, STP

  density = 1.25053*mg/cm3 ;       // STP
  G4Material* Nitrogen = new G4Material(name="N2"  , density, ncomponents=1);
  Nitrogen->AddElement(elN, 2);

 // Oxygen, STP

  density = 1.4289*mg/cm3 ;       // STP
  G4Material* Oxygen = new G4Material(name="O2"  , density, ncomponents=1);
  Oxygen->AddElement(elO, 2);

  /* *****************************
  density = 1.25053*mg/cm3 ;       // STP
  a = 14.01*g/mole ;       // get atomic weight !!!
  //  a = 28.016*g/mole;
  G4Material* N2  = new G4Material(name="Nitrogen", z= 7.,a,density) ;

  density = 1.25053*mg/cm3 ;       // STP
  G4Material* anotherN2 = new G4Material(name="anotherN2", density,ncomponents=2);
  anotherN2->AddElement(elN, 1);
  anotherN2->AddElement(elN, 1);

  // air made from oxigen and nitrogen only

  density = 1.290*mg/cm3;  // old air from elements
  G4Material* air = new G4Material(name="air"  , density, ncomponents=2);
  air->AddElement(elN, fractionmass=0.7);
  air->AddElement(elO, fractionmass=0.3);
  ******************************************** */

  // Dry Air (average composition), STP

  density = 1.2928*mg/cm3 ;       // STP
  G4Material* Air = new G4Material(name="Air"  , density, ncomponents=3);
  Air->AddMaterial( Nitrogen, fractionmass = 0.7557 ) ;
  Air->AddMaterial( Oxygen,   fractionmass = 0.2315 ) ;
  Air->AddMaterial( Argon,    fractionmass = 0.0128 ) ;

////////////////////////////////////////////////////////////////////////////
//
// MWPC mixtures

  // 90% Xe + 10% CH4, STP ; NIM A248 (1986) 379-388

  density = 5.344*mg/cm3 ;      
  G4Material* Xe10CH4 = new G4Material(name="Xe10CH4"  , density, 
                                                  ncomponents=2);
  Xe10CH4->AddMaterial( Xe,       fractionmass = 0.987 ) ;
  Xe10CH4->AddMaterial( metane,   fractionmass = 0.013 ) ;


  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  // default materials of the calorimeter and TR radiator

  fRadiatorMat =  Li ; // CH2 Mylar ; 
  
  fWindowMat = Mylar ;
  fElectrodeMat = Al ;

  AbsorberMaterial = Xe10CH4 ; 
  fGapMat          = Xe10CH4 ;

  WorldMaterial    = He ;
}

/////////////////////////////////////////////////////////////////////////
//
//
  
G4VPhysicalVolume* Em10DetectorConstruction::ConstructCalorimeter()
{
  G4int i, j ; 
  G4double zModule, zRadiator, rModule, rRadiator ; 

  // complete the Calor parameters definition and Print 

  ComputeCalorParameters();
  PrintCalorParameters();
      
  // World
  
  if(solidWorld) delete solidWorld ;
  if(logicWorld) delete logicWorld ;
  if(physiWorld) delete physiWorld ;

  solidWorld = new G4Box("World",				//its name
                   WorldSizeR,WorldSizeR,WorldSizeZ/2.)       ;//its size
                         
  logicWorld = new G4LogicalVolume(solidWorld,		//its solid
                                   WorldMaterial,	//its material
                                   "World");		//its name
                                   
  physiWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 "World",		//its name
                                 logicWorld,		//its logical volume
                                 NULL,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number

  // TR radiator envelope

  G4double radThick = fFoilNumber*(fRadThickness + fGasGap) + fDetGap   ;

  G4double zRad = fStartZ + 0.5*radThick ;
  G4cout<<"zRad = "<<zRad/mm<<" mm"<<G4endl ;

  radThick *= 1.2 ;
  G4cout<<"radThick = "<<radThick/mm<<" mm"<<G4endl ;

  G4Box* solidRadiator = new G4Box("Radiator",1.1*AbsorberRadius , 
                                              1.1*AbsorberRadius, 
                                              0.5*radThick             ) ; 
                         
  logicRadiator = new G4LogicalVolume(solidRadiator,	
                                                       WorldMaterial,      
                                                       "Radiator");	       
                                   
  G4VPhysicalVolume* physiRadiator = new G4PVPlacement(0,
                                     G4ThreeVector(0,0,zRad),	        
                                     "Radiator", logicRadiator,		
                                     physiWorld, false,	0       );  	

  

    fSolidRadSlice = new G4Box("RadSlice",AbsorberRadius,
                                AbsorberRadius,0.5*fRadThickness ) ;

    fLogicRadSlice = new G4LogicalVolume(fSolidRadSlice,fRadiatorMat,
                                          "RadSlice",0,0,0);

    //   fPhysicRadSlice = new G4PVPlacement(0,
    //      G4ThreeVector(0.,0.,fStartZ+1.2*fDetThickness),
    //             "RadSlice",fLogicRadSlice,
    //               physiWorld,false,0);
    
  for(i=0;i<fModuleNumber;i++)
  {
    //   rModule = fStartR + fDetThickness + fDetGap + 
    //           (i-1)*(fFoilNumber*(fRadThickness + fGasGap) + 
    //           fDetThickness + fDetGap) ;

    zModule = fStartZ + fRadThickness + 
              i*( fFoilNumber*(fRadThickness + fGasGap) + 
              fDetThickness + fDetGap )  ;
    G4cout<<"zModule = "<<zModule/mm<<" mm"<<G4endl ;
    G4cout<<"i = "<<i<<"\t"<<G4endl ; 

    for(j=0;j<fFoilNumber;j++)
    {  
      //   rRadiator = rModule + j*(fRadThickness + fGasGap) ;

      zRadiator = zModule + j*(fRadThickness + fGasGap) ;
      G4cout<<zRadiator/mm<<" mm"<<"\t" ;
      //   G4cout<<"j = "<<j<<"\t" ;         
      // RadRing

            
      // fSolidRadRing = new G4Box("RadRing",rRadiator,
      //        rRadiator + fRadThickness,
      //     fDetLength,0.0,360*deg     ) ;

      //  fLogicRadRing = new G4LogicalVolume(fSolidRadRing,fRadiatorMat,
      //                         "radRing",0,0,0);
      
      //  fPhysicRadRing = new G4PVPlacement(0,G4ThreeVector(),
      //       "RadRing",fLogicRadRing,
      //               physiWorld,false,j)  ; 
                                                            
      // We put slice relatively of Radiator, so zRadiator-zRad
      
      fPhysicRadSlice = new G4PVPlacement(0,G4ThreeVector(0.,0.,zRadiator-zRad),
                                         "RadSlice",fLogicRadSlice,
                                          physiRadiator,false,j);
     }                                 
    //   fPhysicDetSlice = new G4PVPlacement(0,
    //          G4ThreeVector(0.,0.,zRadiator+
    //                        fDetGap +0.5*fDetThickness),"DetSlice",
    //                        fLogicDetSlice,physiWorld,false,i); 
  }                                            
  G4cout<<G4endl ;

  G4Box* solidElectrode = new G4Box("Electrode",AbsorberRadius,
                                                AbsorberRadius,
                                                fElectrodeThick/2. ); 
                          
  G4LogicalVolume* logicElectrode = new G4LogicalVolume(solidElectrode,
                                        fElectrodeMat, "Electrode"); 

  G4double zElectrode = zAbsorber - AbsorberThickness/2. - 
                        fElectrodeThick/2. - 0.01*mm;    
      			                  
  //  G4VPhysicalVolume*    physiElectrode = new G4PVPlacement(0,		   
  //  		                       G4ThreeVector(0.,0.,zElectrode),        
  //                                    "Electrode",logicElectrode,
  //                                     physiWorld,false,0);    
  


  G4Box* solidGap = new G4Box("Gap",AbsorberRadius,
                                    AbsorberRadius,
                                    fGapThick/2.     ) ; 
                          
  G4LogicalVolume* logicGap = new G4LogicalVolume(solidGap,fGapMat, "Gap"); 

  G4double zGap = zElectrode - fElectrodeThick/2. - fGapThick/2. - 0.01*mm ;    
      			                  
  // G4VPhysicalVolume*    physiGap = new G4PVPlacement(0,		   
  // 		                       G4ThreeVector(0.,0.,zGap),        
  //                                    "Gap",logicGap,physiWorld,false,0); 

  G4Box* solidWindow = new G4Box("Window",AbsorberRadius,
                                          AbsorberRadius,
                                          fWindowThick/2.  ); 
                          
  G4LogicalVolume* logicWindow = new G4LogicalVolume(solidWindow,
                                     fWindowMat, "Window"); 

  G4double zWindow = zGap - fGapThick/2. - fWindowThick/2. - 0.01*mm ;    
      			                  
  // G4VPhysicalVolume*    physiWindow = new G4PVPlacement(0,		   
  // 		                       G4ThreeVector(0.,0.,zWindow),        
  //                                    "Window",logicWindow,physiWorld,false,0); 
                             
  // Absorber

  if (AbsorberThickness > 0.) 
  { 
      if(solidAbsorber) delete solidAbsorber ;
      if(logicAbsorber) delete logicAbsorber ;
      if(physiAbsorber) delete physiAbsorber ;

      solidAbsorber = new G4Box("Absorber",AbsorberRadius,		
                                           AbsorberRadius,
                                           AbsorberThickness/2.   ); 
                          
      logicAbsorber = new G4LogicalVolume(solidAbsorber,    
      			                  AbsorberMaterial, 
      			                  "Absorber");     
      			                  
      physiAbsorber = new G4PVPlacement(0,		   
      		    G4ThreeVector(0.,0.,zAbsorber),        
                                        "Absorber",        
                                        logicAbsorber,     
                                        physiWorld,       
                                        false,             
                                        0);                
                                        
  }
                                 
  // Sensitive Detectors: Absorber 
  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  if(!calorimeterSD)
  {
    calorimeterSD = new Em10CalorimeterSD("CalorSD",this);
    SDman->AddNewDetector( calorimeterSD );
  }
  if (logicAbsorber)  logicAbsorber->SetSensitiveDetector(calorimeterSD);

  // Parameterisation

  ParametrisationModel();

  // always return physics world

  return physiWorld;
}

////////////////////////////////////////////////////////////////////////////
//
// Construct parametrisation model depending on the current value of fModelNumber

void Em10DetectorConstruction::ParametrisationModel()
{
  G4cout<<"Em10DetectorConstruction::ParametrisationModel() is called"<<G4endl;
  G4cout<<"fModelNumber = "<<fModelNumber<<G4endl;

  G4double alphaPlate = 160.0 ;
  G4double alphaGas   = 160.0 ;

  if(fXTRModel) delete fXTRModel;  

  switch(fModelNumber)
  {
  case 1:
    fXTRModel = new G4FoamXrayTRmodel(logicRadiator,fRadThickness,fGasGap);
    break;
  case 2:
    fXTRModel = new G4FoamXTRdEdx(logicRadiator,fRadThickness,fGasGap);
    break;
  case 3:
    fXTRModel = new G4GamDistrXrayTRmodel(logicRadiator,
  				       fRadThickness,alphaPlate,
                                               fGasGap,alphaGas);
    break;
  case 4:
    fXTRModel = new G4GamDistrXTRdEdx(logicRadiator,
  				       fRadThickness,alphaPlate,
                                               fGasGap,alphaGas);
    break;
  case 5:
    fXTRModel = new G4IrregularXrayTRmodel(logicRadiator,fRadThickness,fGasGap);
    break;
  case 6:
    fXTRModel = new G4IrregularXTRdEdx(logicRadiator,fRadThickness,fGasGap);
    break;
  case 7:
    fXTRModel = new G4PlateIrrGasXrayTRmodel(logicRadiator,
        fRadThickness,fGasGap);
    break;
  case 8:
    fXTRModel = new G4PlateIrrGasXTRdEdx(logicRadiator,fRadThickness,fGasGap);
    break;
  case 9:
    fXTRModel = new G4RegularXrayTRmodel(logicRadiator,fRadThickness,fGasGap);
    break;
  case 10:
    fXTRModel = new G4RegularXTRdEdx(logicRadiator,fRadThickness,fGasGap);
    break;
  default:
    fXTRModel = NULL;
    G4cout<<"Warning: No parametrisation model was defined?"<<G4endl ;
    break ;
  }
  //  fXTRModel->GetPlateZmuProduct() ;
  //  fXTRModel->GetGasZmuProduct() ;

  //  fXTRModel->GetNumberOfPhotons() ;  
  return;
}
////////////////////////////////////////////////////////////////////////////
//
//

void Em10DetectorConstruction::PrintCalorParameters()
{
  G4cout << "\n The  WORLD   is made of " 
       << WorldSizeZ/mm << "mm of " << WorldMaterial->GetName() ;
  G4cout << ", the transverse size (R) of the world is " << WorldSizeR/mm << " mm. " << G4endl;
  G4cout << " The ABSORBER is made of " 
       << AbsorberThickness/mm << "mm of " << AbsorberMaterial->GetName() ;
  G4cout << ", the transverse size (R) is " << AbsorberRadius/mm << " mm. " << G4endl;
  G4cout << " Z position of the (middle of the) absorber " << zAbsorber/mm << "  mm." << G4endl;
  G4cout << G4endl;
}

///////////////////////////////////////////////////////////////////////////
//
//

void Em10DetectorConstruction::SetAbsorberMaterial(G4String materialChoice)
{
  // get the pointer to the material table
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  // search the material by its name   
  G4Material* pttoMaterial;
  for (G4int J=0 ; J<theMaterialTable->length() ; J++)
   { pttoMaterial = (*theMaterialTable)(J);     
     if(pttoMaterial->GetName() == materialChoice)
        {AbsorberMaterial = pttoMaterial;
         logicAbsorber->SetMaterial(pttoMaterial); 
        // PrintCalorParameters();
        }             
   }
}
///////////////////////////////////////////////////////////////////////////
//
//

void Em10DetectorConstruction::SetRadiatorMaterial(G4String materialChoice)
{
  // get the pointer to the material table

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  // search the material by its name
   
  G4Material* pttoMaterial;
  for (G4int J=0 ; J<theMaterialTable->length() ; J++)
  { 
    pttoMaterial = (*theMaterialTable)(J);
     
    if(pttoMaterial->GetName() == materialChoice)
    {
      AbsorberMaterial = pttoMaterial;
      logicRadiator->SetMaterial(pttoMaterial); 
      // PrintCalorParameters();
    }             
   }
}

////////////////////////////////////////////////////////////////////////////
//
//

void Em10DetectorConstruction::SetWorldMaterial(G4String materialChoice)
{
  // get the pointer to the material table
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  // search the material by its name   
  G4Material* pttoMaterial;
  for (G4int J=0 ; J<theMaterialTable->length() ; J++)
   { pttoMaterial = (*theMaterialTable)(J);     
     if(pttoMaterial->GetName() == materialChoice)
        {WorldMaterial = pttoMaterial;
         logicWorld->SetMaterial(pttoMaterial); 
       //  PrintCalorParameters();
        }             
   }
}

///////////////////////////////////////////////////////////////////////////
//
//

void Em10DetectorConstruction::SetAbsorberThickness(G4double val)
{
  // change Absorber thickness and recompute the calorimeter parameters
  AbsorberThickness = val;
  ComputeCalorParameters();
}  

///////////////////////////////////////////////////////////////////////////
//
//

void Em10DetectorConstruction::SetRadiatorThickness(G4double val)
{
  // change XTR radiator thickness and recompute the calorimeter parameters
  fRadThickness = val;
  // ComputeCalorParameters();
}
  
///////////////////////////////////////////////////////////////////////////
//
//

void Em10DetectorConstruction::SetGasGapThickness(G4double val)
{
  // change XTR gas gap thickness and recompute the calorimeter parameters
  fGasGap = val;
  // ComputeCalorParameters();
}  

/////////////////////////////////////////////////////////////////////////////
//
//

void Em10DetectorConstruction::SetAbsorberRadius(G4double val)
{
  // change the transverse size and recompute the calorimeter parameters
  AbsorberRadius = val;
  ComputeCalorParameters();
}  

////////////////////////////////////////////////////////////////////////////
//
//

void Em10DetectorConstruction::SetWorldSizeZ(G4double val)
{
  worldchanged=true;
  WorldSizeZ = val;
  ComputeCalorParameters();
}  

///////////////////////////////////////////////////////////////////////////
//
//

void Em10DetectorConstruction::SetWorldSizeR(G4double val)
{
  worldchanged=true;
  WorldSizeR = val;
  ComputeCalorParameters();
}  

//////////////////////////////////////////////////////////////////////////////
//
//

void Em10DetectorConstruction::SetAbsorberZpos(G4double val)
{
  zAbsorber  = val;
  ComputeCalorParameters();
}  

//////////////////////////////////////////////////////////////////////////////
//
//

void Em10DetectorConstruction::SetMagField(G4double fieldValue)
{
  //apply a global uniform magnetic field along X axis

  /* *********************************************************

  G4FieldManager* fieldMgr 
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    
  if(magField) delete magField;		//delete the existing magn field
  
  if(fieldValue!=0.)			// create a new one if non nul
  { 
    magField = new G4UniformMagField(G4ThreeVector(fieldValue,0.,0.));        
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);
  } 
  else 
  {
    magField = NULL;
    fieldMgr->SetDetectorField(magField);
  }

  *************************************************************** */

}

///////////////////////////////////////////////////////////////////////////////
//
//
  
void Em10DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructCalorimeter());
}

//
//
////////////////////////////////////////////////////////////////////////////

















