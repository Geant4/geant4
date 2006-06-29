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
// 
//
//  
//
//  Test routine for G4VXTRenergyLoss class and inherited from it code
//
// History:
//
// 20.01.05, V. Grichine 

#include "G4ios.hh"
#include <fstream>
#include <cmath>
#include "globals.hh"
#include "Randomize.hh"
#include "G4UnitsTable.hh"


#include <iomanip>

#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4RegionStore.hh"
#include "G4MaterialTable.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"

#include "G4VXTRenergyLoss.hh"
#include "G4RegularXTRadiator.hh"
#include "G4TransparentRegXTRadiator.hh"
#include "G4GammaXTRadiator.hh"
#include "G4StrawTubeXTRadiator.hh"

#include "G4XTRGammaRadModel.hh"
#include "G4XTRRegularRadModel.hh"
#include "G4XTRTransparentRegRadModel.hh"

#include "G4SynchrotronRadiation.hh"

#include "G4ParticleDefinition.hh"
#include "G4Proton.hh"



int main()
{

  /*
  std::ofstream outdEdx("XTRdEdx.out", std::ios::out ) ;
  outdEdx.setf( std::ios::scientific, std::ios::floatfield );

  std::ofstream outdNdx("XTRdNdx.out", std::ios::out ) ;
  outdNdx.setf( std::ios::scientific, std::ios::floatfield );

  std::ofstream outXsc("InitXsc.out", std::ios::out ) ;
  outXsc.setf( std::ios::scientific, std::ios::floatfield );

  //  std::ifstream fileRead("exp.dat", std::ios::out ) ;
  //  fileRead.setf( std::ios::scientific, std::ios::floatfield );

  std::ofstream fileWrite1("mpXTR.dat", std::ios::out ) ;
  fileWrite1.setf( std::ios::scientific, std::ios::floatfield );
  */
 

  /////////////////////////////////////////////////////////////////
  //
  // Create materials  

 
  G4String name, symbol ;                            //a =mass of a mole;
  G4double a, z ;   //z =mean number of protons;  
  G4double density, foilDensity, gasDensity, totDensity ;
  G4double fractionFoil, fractionGas ;
  G4int nel ; 
  
  //G4int ncomponents, natoms;
  G4int  ncomponents;
  //G4double abundance, fractionmass;
  G4double fractionmass;
  //G4double temperature, pressure;
  
  /////////////////////////////////////
  //
  // define Elements
  
   
  a = 1.01*g/mole;
  G4Element* elH  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);
  
  a = 6.94*g/mole;
  G4Element* elLi = new G4Element(name="Lithium",symbol="Li" , z= 3., a);
  
  a = 9.01*g/mole;
  G4Element* elBe = new G4Element(name="Berillium",symbol="Be" , z= 4., a);
  
  a = 12.01*g/mole;
  G4Element* elC  = new G4Element(name="Carbon", symbol="C", z=6., a);
  
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

  G4Material* CH2 = new G4Material ("CH2" , 0.91*g/cm3, 2);
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

  // Carbon dioxide, STP 

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
  
  // 80% Xe + 20% CO2, STP
  
  density = 5.0818*mg/cm3 ;      
  G4Material* Xe20CO2 = new G4Material(name="Xe20CO2"  , density, ncomponents=2);
  Xe20CO2->AddMaterial( Xe,              fractionmass = 0.922 ) ;
  Xe20CO2->AddMaterial( CarbonDioxide,   fractionmass = 0.078 ) ;
  
  // 80% Kr + 20% CO2, STP
  
  density = 3.601*mg/cm3 ;      
  G4Material* Kr20CO2 = new G4Material(name="Kr20CO2", density, 
                                       ncomponents=2);
  Kr20CO2->AddMaterial( Kr,              fractionmass = 0.89 ) ;
  Kr20CO2->AddMaterial( CarbonDioxide,   fractionmass = 0.11 ) ;
  
  // Xe + 55% He + 15% CH4 ; NIM A294 (1990) 465-472; STP
  
  density = 1.963*mg/cm3;
  G4Material* Xe55He15CH4 = new G4Material(name="Xe55He15CH4",density,
                                           ncomponents=3);
  Xe55He15CH4->AddMaterial(Xe, 0.895);
  Xe55He15CH4->AddMaterial(He, 0.050);
  Xe55He15CH4->AddMaterial(metane,0.055);

  // 90% Xe + 10% CH4, STP ; NIM A248 (1986) 379-388
  
  density = 5.344*mg/cm3 ;      
  G4Material* Xe10CH4 = new G4Material(name="Xe10CH4"  , density, 
                                       ncomponents=2);
  Xe10CH4->AddMaterial( Xe,       fractionmass = 0.987 ) ;
  Xe10CH4->AddMaterial( metane,   fractionmass = 0.013 ) ;
  
  // 95% Xe + 5% CH4, STP ; NIM A214 (1983) 261-268
  
  density = 5.601*mg/cm3 ;      
  G4Material* Xe5CH4 = new G4Material(name="Xe5CH4"  , density, 
                                      ncomponents=2);
  Xe5CH4->AddMaterial( Xe,       fractionmass = 0.994 ) ;
  Xe5CH4->AddMaterial( metane,   fractionmass = 0.006 ) ;
  
  // 80% Xe + 20% CH4, STP ; NIM A253 (1987) 235-244
  
  density = 4.83*mg/cm3 ;      
  G4Material* Xe20CH4 = new G4Material(name="Xe20CH4"  , density, 
                                       ncomponents=2);
  Xe20CH4->AddMaterial( Xe,       fractionmass = 0.97 ) ;
  Xe20CH4->AddMaterial( metane,   fractionmass = 0.03 ) ;

  // 93% Ar + 7% CH4, STP ; NIM 107 (1973) 413-422

  density = 1.709*mg/cm3 ;      
  G4Material* Ar7CH4 = new G4Material(name="Ar7CH4"  , density, 
                                                  ncomponents=2);
  Ar7CH4->AddMaterial( Argon,       fractionmass = 0.971 ) ;
  Ar7CH4->AddMaterial( metane,   fractionmass = 0.029 ) ;

  // 93% Kr + 7% CH4, STP ; NIM 107 (1973) 413-422

  density = 3.491*mg/cm3 ;      
  G4Material* Kr7CH4 = new G4Material(name="Kr7CH4"  , density, 
                                                  ncomponents=2);
  Kr7CH4->AddMaterial( Kr,       fractionmass = 0.986 ) ;
  Kr7CH4->AddMaterial( metane,   fractionmass = 0.014 ) ;

  // 0.5*(95% Xe + 5% CH4)+0.5*(93% Ar + 7% CH4), STP ; NIM A214 (1983) 261-268

  density = 3.655*mg/cm3 ;      
  G4Material* XeArCH4 = new G4Material(name="XeArCH4"  , density, 
                                                  ncomponents=2);
  XeArCH4->AddMaterial( Xe5CH4,       fractionmass = 0.766 ) ;
  XeArCH4->AddMaterial( Ar7CH4,   fractionmass = 0.234 ) ;


  ////////////////////////////////////////////////////////////
  //
  // Geometry 


  G4double foilThick = 0.02*mm ; // 25*micrometer ;     
  G4double gasGap       = 0.50*mm ;          //  1500*micrometer  ;   
  G4double foilGasRatio = foilThick/(foilThick+gasGap) ;
  G4int    foilNumber   = 120 ;             //  188 ;
  G4double detGap       = 0.01*mm ;


  G4double alphaPlate   = 2.0 ;
  G4double alphaGas     = 10.0 ;
  G4int    modelNumber  = 0 ;

// TR radiator envelope

  G4double radThick = foilNumber*(foilThick + gasGap) - gasGap + detGap;

  G4double absorberRadius    = 10.*cm;

  G4double absorberThickness = 15.0*mm ;   // 40.0*mm ;

  // Preparation of mixed radiator material

  foilDensity = 0.91*g/cm3 ;// CH2 1.39*g/cm3; // Mylar 0.534*g/cm3; //Li        
  gasDensity  = 1.2928*mg/cm3 ; // Air 0.178*mg/cm3 ; // He 1.977*mg/cm3; // CO2
 
  totDensity  = foilDensity*foilGasRatio + gasDensity*(1.0-foilGasRatio) ;

  fractionFoil =  foilDensity*foilGasRatio/totDensity ; 
  fractionGas  =  gasDensity*(1.0-foilGasRatio)/totDensity ;  
    
  G4Material* radiatorMat = new G4Material(name="radiatorMat"  , totDensity, 
                                                  ncomponents=2);
  radiatorMat->AddMaterial( CH2, fractionFoil ) ;
  radiatorMat->AddMaterial( Air, fractionGas  ) ;

  //  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  // default materials of the calorimeter and TR radiator

  G4Material* fRadiatorMat = radiatorMat ; // CH2 Mylar ; 
  G4Material* foilMat     = CH2 ; // Li ; // CH2 ;  
  G4Material* gasMat      = Air ; // He;  //  CO2 ; 
  G4Material* absMat      = Xe20CO2 ; // He ;// CO2 ; 
  
  // fWindowMat = Mylar ;
  // fElectrodeMat = Al ;

  // AbsorberMaterial = Ar7CH4; // Xe10CH4; // Xe55He15CH4;
 

  // fGapMat          = Xe10CH4 ;

  // WorldMaterial    = Air ; // CO2 ;  



  ///////////////////////

  G4int i, j, k, nBin, numOfMaterials, iSan, nbOfElements, sanIndex, row ;

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable() ;

  numOfMaterials = theMaterialTable->size();

  G4String testName;

  for( k = 0; k < numOfMaterials; k++ )
  {
    //    if((*theMaterialTable)[k]->GetName() != testName) continue ;

    // outFile << "Material : " <<(*theMaterialTable)[k]->GetName() << G4endl ;
    // G4cout <<k<<"\t"<< "Material : " <<(*theMaterialTable)[k]->GetName() << G4endl ;
  }

  //  G4cout<<"Enter material name for test : "<<std::flush ;
  //  G4cin>>testName ;

    
  // G4Region* regGasDet = new G4Region("VertexDetector");
  // regGasDet->AddRootLogicalVolume(logicAbsorber);
   
  G4ProductionCuts* cuts = new G4ProductionCuts();
  cuts->SetProductionCut(10.*mm,"gamma");
  cuts->SetProductionCut(1.*mm,"e-");
  cuts->SetProductionCut(1.*mm,"e+");

   // regGasDet->SetProductionCuts(cuts);
   
  G4cout.precision(4);

  //  G4MaterialCutsCouple* matCC = new G4MaterialCutsCouple( 
  //                                  (*theMaterialTable)[k], cuts);

  G4cout<<"radThick = "     <<radThick/mm<<" mm"<<G4endl ;
  G4cout<<"foilNumber = "  <<foilNumber<<G4endl ;
  G4cout<<"radiatorMat = " <<radiatorMat->GetName()<<G4endl ;
 

  G4Box* solidRadiator = new G4Box("Radiator",1.1*absorberRadius , 
                                              1.1*absorberRadius, 
                                              0.5*radThick             ) ; 
                         
  G4LogicalVolume* logicRadiator = new G4LogicalVolume(solidRadiator,    
                                                       radiatorMat,      
                                                       "Radiator");            
     

  // const G4RegionStore* theRegionStore = G4RegionStore::GetInstance();
  // G4Region* gas = theRegionStore->GetRegion("XTRdEdxDetector");

  G4VXTRenergyLoss* processXTR;

  const G4ParticleDefinition proton(
                 name,   0.9382723*GeV,       0.0*MeV,       eplus,
                    1,              +1,             0,
                    1,              +1,             0,
             "baryon",               0,            +1,        2212,
                 true,            -1.0,          NULL,
             false,           "neucleon"
              );

  G4ParticleDefinition* theProton = G4Proton::ProtonDefinition();
  // *proton = theProton;  

  // G4String fXTRModel = "transpR";
  G4String fXTRModel = "transpM";

  // G4String fXTRModel = "regR";
  // G4String fXTRModel = "regM";

  // G4String fXTRModel = "gammaR";
  // G4String fXTRModel = "gammaM";

  // G4String fXTRModel = "strawR";

  if(fXTRModel == "gammaR" )
  {
    // G4GammaXTRadiator*
    processXTR = new G4GammaXTRadiator(logicRadiator,
                                       1000.,
                                       100.,
                                       foilMat,
                                       gasMat,
                                       foilThick,
                                       gasGap,
                                       foilNumber,
                                       "GammaXTRadiator");
  }
  else if(fXTRModel == "gammaM" )
  {
    // G4XTRGammaRadModel*
    processXTR = new G4XTRGammaRadModel(logicRadiator,
                                       1000.,
                                       100.,
                                       foilMat,
                                       gasMat,
                                       foilThick,
                                       gasGap,
                                       foilNumber,
                                       "GammaXTRmodel");
  }
  else if(fXTRModel == "strawR" )
  {

    // G4StrawTubeXTRadiator*
    processXTR = new G4StrawTubeXTRadiator(logicRadiator,
                                           foilMat,
                                           gasMat,
                                0.53,      // foilThick,
                                3.14159,   // gasGap,
                                         absMat,
                                         true,
                                         "strawXTRadiator");
  }
  else if(fXTRModel == "regR" )
  {
    // G4RegularXTRadiator*
    processXTR = new G4RegularXTRadiator(logicRadiator,
                                         foilMat,
                                         gasMat,
                                         foilThick,
                                         gasGap,
                                         foilNumber,
                                         "RegularXTRadiator");
  }
  else if(fXTRModel == "transpR" )
  {
    // G4TransparentRegXTRadiator*
    processXTR = new G4TransparentRegXTRadiator(logicRadiator,
                                                foilMat,
                                                gasMat,
                                                foilThick,
                                                gasGap,
                                                foilNumber,
                                                "TranspRegXTRadiator");
  }
  else if(fXTRModel == "regM" )
  {
    // G4XTRRegularRadModel*
    processXTR = new G4XTRRegularRadModel(logicRadiator,
                                          foilMat,
                                          gasMat,
                                          foilThick,
                                          gasGap,
                                          foilNumber,
                                         "RegularXTRmodel");

  }
  else if(fXTRModel == "transpM" )
  {
    // G4XTRTransparentRegRadModel*
    processXTR = new G4XTRTransparentRegRadModel(logicRadiator,
                                                 foilMat,
                                                 gasMat,
                                                 foilThick,
                                                 gasGap,
                                                 foilNumber,
                                         "TranspRegXTRmodel");
  }
  else
  {
    G4Exception("Invalid XTR model name", "InvalidSetup",
                 FatalException, "XTR model name is out of the name list");
  }
  processXTR->SetVerboseLevel(1);
  // processXTR->SetAngleRadDistr(true);

  // processXTR->BuildPhysicsTable(proton);

  // processXTR->SetVerboseLevel(1);

  static G4int totBin = processXTR->GetTotBin();
  nBin = totBin;
  G4cout<<"totBin = "<<totBin<<G4endl;

  // test of XTR table step do-it


  G4double energyTR = 10*keV, cofAngle = 5.1,  angle2, dNdA, xCompton, lambdaC;
  G4double charge = 1.0;
  G4double chargeSq  = charge*charge ;
  G4double gamma     = 4.e4; 
  G4cout<<"gamma = "<<gamma<<G4endl;
  G4cout<<"energyTR = "<<energyTR/keV<<" keV"<<G4endl;
  
  processXTR->SetGamma(gamma);

  // processXTR->GetAngleVector(energyTR,nBin);

  // xCompton = processXTR->GetGasCompton(energyTR);

  // lambdaC = 1./xCompton;

  // G4cout<<"lambdaC = "<<lambdaC/m <<" m; for energy = "<<energyTR/keV<<" keV"<<G4endl;

  // G4double dNdA = processXTR->SpectralXTRdEdx(energyTR);

  // G4double angle2 = cofAngle*cofAngle/gamma/gamma;

  // G4double dNdAngle = processXTR-> AngleXTRdEdx(angle2);
  /*
  for(i = 0; i < 40; i++ )
  {
    angle2 = processXTR->GetRandomAngle(energyTR,40);
    G4cout<<"random theta*gamma = "<<std::sqrt(angle2)*gamma<<G4endl;
  }
  */

  /*
  for(i = 0; i < 40; i++ )
  {
     cofAngle = 0.5*i;
     G4double angle2 = cofAngle*cofAngle/gamma/gamma;
     G4double dNdAngle = processXTR-> AngleXTRdEdx(angle2);
     dNdAngle *=fine_structure_const/pi;
     G4cout<<"cofAngle = "<<cofAngle<<"; angle = "<<cofAngle/gamma
           <<"; dNdAngle = "<<dNdAngle<<G4endl;
  }
  */


  G4int iTkin;
  G4cout<<"gamma = "<<gamma<<G4endl;

  G4double TkinScaled = (gamma - 1.)*proton_mass_c2;

  /*  
  for( iTkin = 0; iTkin < totBin; iTkin++ )
  {
    if(TkinScaled < processXTR->GetProtonVector()->
                                        GetLowEdgeEnergy(iTkin)) break;    
  }

  G4double xtrEnergy[100];
  G4int spectrum[100];


  for( k = 0; k < 100; k++ )
  {
    xtrEnergy[k] = (1.0+ 1.0*k)*keV;
    spectrum[k]  = 0;
  }
  

  for( i = 0; i < 10000; i++ )
  {
    energyTR = processXTR->GetXTRrandomEnergy(TkinScaled,iTkin);

    for( k = 0; k < 100; k++ )
    {
      if( energyTR <= xtrEnergy[k] ) break;    
    }
    spectrum[k] += 1;
  }

  // output to file


  if(fXTRModel == "gammaR" )
  {
    std::ofstream fileWrite("gammaR.dat", std::ios::out ) ;
    fileWrite.setf( std::ios::scientific, std::ios::floatfield );

    for( k = 0; k < 41; k++ )
    {    
      G4cout<<k<<"\t"<<xtrEnergy[k]/keV<<"\t"<<spectrum[k]<<G4endl;
      fileWrite<<xtrEnergy[k]/keV<<"\t"<<spectrum[k]<<G4endl;  
    }
  }
  else if(fXTRModel == "gammaM" )
  {
    std::ofstream fileWrite("gammaM.dat", std::ios::out ) ;
    fileWrite.setf( std::ios::scientific, std::ios::floatfield );

    for( k = 0; k < 41; k++ )
    {    
      G4cout<<k<<"\t"<<xtrEnergy[k]/keV<<"\t"<<spectrum[k]<<G4endl;
      fileWrite<<xtrEnergy[k]/keV<<"\t"<<spectrum[k]<<G4endl;  
    }
  }
  else if(fXTRModel == "strawR" )
  {
    std::ofstream fileWrite("strawR.dat", std::ios::out ) ;
    fileWrite.setf( std::ios::scientific, std::ios::floatfield );

    for( k = 0; k < 41; k++ )
    {    
      G4cout<<k<<"\t"<<xtrEnergy[k]/keV<<"\t"<<spectrum[k]<<G4endl;
      fileWrite<<xtrEnergy[k]/keV<<"\t"<<spectrum[k]<<G4endl;  
    }
  }
  else if(fXTRModel == "regR" )
  {
    std::ofstream fileWrite("regR.dat", std::ios::out ) ;
    fileWrite.setf( std::ios::scientific, std::ios::floatfield );

    for( k = 0; k < 41; k++ )
    {    
      G4cout<<k<<"\t"<<xtrEnergy[k]/keV<<"\t"<<spectrum[k]<<G4endl;
      fileWrite<<xtrEnergy[k]/keV<<"\t"<<spectrum[k]<<G4endl;  
    }
  }
  else if(fXTRModel == "transpR" )
  {
    std::ofstream fileWrite("transpR.dat", std::ios::out ) ;
    fileWrite.setf( std::ios::scientific, std::ios::floatfield );

    for( k = 0; k < 41; k++ )
    {    
      G4cout<<k<<"\t"<<xtrEnergy[k]/keV<<"\t"<<spectrum[k]<<G4endl;
      fileWrite<<xtrEnergy[k]/keV<<"\t"<<spectrum[k]<<G4endl;  
    }
  }
  else if(fXTRModel == "regM" )
  {
    std::ofstream fileWrite("regM.dat", std::ios::out ) ;
    fileWrite.setf( std::ios::scientific, std::ios::floatfield );

    for( k = 0; k < 41; k++ )
    {    
      G4cout<<k<<"\t"<<xtrEnergy[k]/keV<<"\t"<<spectrum[k]<<G4endl;
      fileWrite<<xtrEnergy[k]/keV<<"\t"<<spectrum[k]<<G4endl;  
    }
  }
  else if(fXTRModel == "transpM" )
  {
    std::ofstream fileWrite("transpM.dat", std::ios::out ) ;
    fileWrite.setf( std::ios::scientific, std::ios::floatfield );

    for( k = 0; k < 41; k++ )
    {    
      G4cout<<k<<"\t"<<xtrEnergy[k]/keV<<"\t"<<spectrum[k]<<G4endl;
      fileWrite<<xtrEnergy[k]/keV<<"\t"<<spectrum[k]<<G4endl;  
    }
  }
  else
  {
    G4Exception("Invalid XTR model name, no output file", "InvalidSetup",
                 FatalException, "XTR model name is out of the name list");
  }
*/



  G4cout.precision(12);
  G4double ksi, prob;
  G4SynchrotronRadiation* sr = new G4SynchrotronRadiation();
  // sr->SetRootNumber(100);
  // ksi = 1.e-8;
  // ksi = 0.;
  prob = sr->GetIntProbSR( ksi);
  G4cout<<"ksi = "<<ksi<<"; SR probability = "<<prob<<G4endl;
  return 1 ;
}

