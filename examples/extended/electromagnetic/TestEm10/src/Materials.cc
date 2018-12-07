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
/// \file electromagnetic/TestEm10/src/Materials.cc
/// \brief Implementation of the Materials class
//
//
// 
//      GEANT 4 class 
//
//      History: based on object model of
//       Materials
//     Originally Created in Test30 by Vladimir Ivanchenko, 12 March 2002 
// 
//    Modified for Test by V. Grichine, 29 Jan 2006
//    is filled with XTR related materials, plastics, gas mixtures, etc


#include "Materials.hh"

#include "G4UnitsTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4SystemOfUnits.hh"

Materials* Materials::fgInstance = 0;

Materials* Materials::GetInstance()
{
  if ( ! fgInstance ) {
    fgInstance = new Materials();
  }
  return fgInstance;
}


Materials::Materials()
{
  fgInstance = this;
  Initialise();
}

Materials::~Materials()
{}

void Materials::Initialise()
{
  G4String name, symbol;                           
  G4double a, z;  
  G4double density, fractionmass; 
  G4int nel, ncomponents;

  // define Elements
 
  a = 1.01*g/mole;
  G4Element* elH  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);

  a = 6.94*g/mole;
  G4Element* elLi  = new G4Element(name="Lithium",symbol="Li" , z= 3., a);

  a = 9.01*g/mole;
  G4Element* elBe  = new G4Element(name="Berillium",symbol="Be" , z= 4., a);

  a = 12.01*g/mole;
  G4Element* elC = new G4Element(name="Carbon", symbol="C", z=6., a);

  a = 14.01*g/mole;
  G4Element* elN  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

  a = 16.00*g/mole;
  G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

  a = 39.948*g/mole;
  G4Element* elAr = new G4Element(name="Argon", symbol="Ar", z=18., a);

  /*
  a = 131.29*g/mole;
  G4Element* elXe = new G4Element(name="Xenon", symbol="Xe", z=54., a);

  a = 19.00*g/mole;
  G4Element* elF  = new G4Element(name="Fluorine", symbol="F", z=9., a);
  */

  //////////////
  //
  // Detector windows, electrodes
  // Al for electrodes

  density = 2.700*g/cm3;
  a = 26.98*g/mole;
  new G4Material(name="Al", z=13., a, density);


  /////////
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

  // Kapton Dupont de Nemur (density: 1.396-1.430, get middle )

  density = 1.413*g/cm3;
  G4Material* Kapton = new G4Material(name="Kapton", density, nel=4);
  Kapton->AddElement(elO,5);
  Kapton->AddElement(elC,22);
  Kapton->AddElement(elN,2);
  Kapton->AddElement(elH,10);

  // Kapton (polyimide) ??? since = Mylar C5H4O2

  // density = 1.39*g/cm3;
  // G4Material* kapton = new G4Material(name="kapton", density, nel=3);
  // Kapton->AddElement(elO,2);
  // Kapton->AddElement(elC,5);
  // Kapton->AddElement(elH,4);

  // Polypropelene

  G4Material* CH2 = new G4Material ("CH2" , 0.91*g/cm3, 2);
  CH2->AddElement(elH,2);
  CH2->AddElement(elC,1);

  ////////////////////////////
  //
  // Noble gases , STP conditions

  // Helium as detector gas, STP

  density = 0.178*mg/cm3;
  a = 4.0026*g/mole;
  G4Material* He  = new G4Material(name="He",z=2., a, density );

  // Neon as detector gas, STP

  density = 0.900*mg/cm3;
  a = 20.179*g/mole;
  new G4Material(name="Ne",z=10., a, density );

  // Argon as detector gas, STP

  density = 1.7836*mg/cm3;       // STP
  G4Material* Argon = new G4Material(name="Argon"  , density, ncomponents=1);
  Argon->AddElement(elAr, 1);

  // Krypton as detector gas, STP

  density = 3.700*mg/cm3;
  a = 83.80*g/mole;
  G4Material* Kr  = new G4Material(name="Kr",z=36., a, density );

  // Xenon as detector gas, STP

  density = 5.858*mg/cm3;
  a = 131.29*g/mole;
  G4Material* Xe  = new G4Material(name="Xenon",z=54., a, density );

/////////////////////////////////
//
// Hydrocarbones, metane and others

  // Metane, STP

  density = 0.7174*mg/cm3;
  G4Material* metane  = new G4Material(name="CH4",density,nel=2);
  metane->AddElement(elC,1);
  metane->AddElement(elH,4);

  // Propane, STP

  density = 2.005*mg/cm3 ;
  G4Material* propane = new G4Material(name="C3H8",density,nel=2);
  propane->AddElement(elC,3);
  propane->AddElement(elH,8);

  // iso-Butane (methylpropane), STP

  density = 2.67*mg/cm3;
  G4Material* isobutane = new G4Material(name="isoC4H10",density,nel=2);
  isobutane->AddElement(elC,4);
  isobutane->AddElement(elH,10);

  /////////////////////////
  //
  // Molecular gases

  // Carbon dioxide, STP

  density = 1.977*mg/cm3;
  G4Material* CO2 = new G4Material(name="CO2", density, nel=2,
                                       kStateGas,273.15*kelvin,1.*atmosphere);
  CO2->AddElement(elC,1);
  CO2->AddElement(elO,2);

  // Carbon dioxide, STP

  density = 1.977*273.*mg/cm3/293.;
  G4Material* CarbonDioxide = new G4Material(name="CO2_2", density, nel=2);
  CarbonDioxide->AddElement(elC,1);
  CarbonDioxide->AddElement(elO,2);

  // Nitrogen, STP

  density = 1.25053*mg/cm3;       // STP
  G4Material* Nitrogen = new G4Material(name="N2"  , density, ncomponents=1);
  Nitrogen->AddElement(elN, 2);

  // Oxygen, STP

  density = 1.4289*mg/cm3;       // STP
  G4Material* Oxygen = new G4Material(name="O2"  , density, ncomponents=1);
  Oxygen->AddElement(elO, 2);

  /* *****************************
  density = 1.25053*mg/cm3;       // STP
  a = 14.01*g/mole ;       // get atomic weight !!!
  //  a = 28.016*g/mole;
  G4Material* N2  = new G4Material(name="Nitrogen", z= 7.,a,density) ;

  density = 1.25053*mg/cm3;       // STP
  G4Material* anotherN2 = new G4Material(name="anotherN2", density,ncomponents=2);
  anotherN2->AddElement(elN, 1);
  anotherN2->AddElement(elN, 1);

  // air made from oxigen and nitrogen only

  density = 1.290*mg/cm3;  // old air from elements
  G4Material* air = new G4Material(name="air"  , density, ncomponents=2);
  air->AddElement(elN, fractionmass=0.7);
  air->AddElement(elO, fractionmass=0.3);
  ******************************************** */

  // Dry Air (average composition with Ar), STP

  density = 1.2928*mg/cm3 ;       // STP
  G4Material* Air = new G4Material(name="Air"  , density, ncomponents=3);
  Air->AddMaterial( Nitrogen, fractionmass = 0.7557 );
  Air->AddMaterial( Oxygen,   fractionmass = 0.2315 );
  Air->AddMaterial( Argon,    fractionmass = 0.0128 );

  ////////////////////////////////////////////////////////////////////////////
  //
  // MWPC mixtures

  // 85% Xe + 15% CO2, STP

  density = 4.9*mg/cm3;
  G4Material* Xe15CO2 = new G4Material(name="Xe15CO2"  , density, ncomponents=2);
  Xe15CO2->AddMaterial( Xe,              fractionmass = 0.979);
  Xe15CO2->AddMaterial( CarbonDioxide,   fractionmass = 0.021);

  // 80% Xe + 20% CO2, STP

  density = 5.0818*mg/cm3;
  G4Material* Xe20CO2 = new G4Material(name="Xe20CO2"  , density, ncomponents=2);
  Xe20CO2->AddMaterial( Xe,              fractionmass = 0.922 );
  Xe20CO2->AddMaterial( CarbonDioxide,   fractionmass = 0.078 );

  // 70% Xe + 27% CO2 + 3% O2, 20 1 atm ATLAS straw tube mixture

  density = 4.358*mg/cm3;
  G4Material* Xe27CO23O2 = new G4Material(name="Xe27CO23O2"  , density, ncomponents=3);
  Xe27CO23O2->AddMaterial( Xe,            fractionmass = 0.87671);
  Xe27CO23O2->AddMaterial( CarbonDioxide, fractionmass = 0.11412);
  Xe27CO23O2->AddMaterial( Oxygen,        fractionmass = 0.00917);

  // 80% Kr + 20% CO2, STP

  density = 3.601*mg/cm3;
  G4Material* Kr20CO2 = new G4Material(name="Kr20CO2", density,
                                       ncomponents=2);
  Kr20CO2->AddMaterial( Kr,              fractionmass = 0.89 );
  Kr20CO2->AddMaterial( CarbonDioxide,   fractionmass = 0.11 );

  // Xe + 55% He + 15% CH4 ; NIM A294 (1990) 465-472; STP

  density = 1.963*273.*mg/cm3/293.;
  G4Material* Xe55He15CH4 = new G4Material(name="Xe55He15CH4",density,
                                           ncomponents=3);
  Xe55He15CH4->AddMaterial(Xe, 0.895);
  Xe55He15CH4->AddMaterial(He, 0.050);
  Xe55He15CH4->AddMaterial(metane,0.055);

  // 90% Xe + 10% CH4, STP ; NIM A248 (1986) 379-388

  density = 5.344*mg/cm3;
  G4Material* Xe10CH4 = new G4Material(name="Xe10CH4"  , density,
                                       ncomponents=2);
  Xe10CH4->AddMaterial( Xe,       fractionmass = 0.987 ) ;
  Xe10CH4->AddMaterial( metane,   fractionmass = 0.013 ) ;

  // 95% Xe + 5% CH4, STP ; NIM A214 (1983) 261-268

  density = 5.601*mg/cm3;
  G4Material* Xe5CH4 = new G4Material(name="Xe5CH4"  , density,
                                      ncomponents=2);
  Xe5CH4->AddMaterial( Xe,       fractionmass = 0.994 );
  Xe5CH4->AddMaterial( metane,   fractionmass = 0.006 );

  // 80% Xe + 20% CH4, STP ; NIM A253 (1987) 235-244

  density = 4.83*mg/cm3;
  G4Material* Xe20CH4 = new G4Material(name="Xe20CH4"  , density,
                                       ncomponents=2);
  Xe20CH4->AddMaterial( Xe,       fractionmass = 0.97 );
  Xe20CH4->AddMaterial( metane,   fractionmass = 0.03 );

  // 93% Ar + 7% CH4, STP ; NIM 107 (1973) 413-422

  density = 1.709*mg/cm3;
  G4Material* Ar7CH4 = new G4Material(name="Ar7CH4"  , density,
                                                  ncomponents=2);
  Ar7CH4->AddMaterial( Argon,       fractionmass = 0.971 );
  Ar7CH4->AddMaterial( metane,   fractionmass = 0.029 );

  // 93% Kr + 7% CH4, STP ; NIM 107 (1973) 413-422

  density = 3.491*mg/cm3;
  G4Material* Kr7CH4 = new G4Material(name="Kr7CH4"  , density,
                                                  ncomponents=2);
  Kr7CH4->AddMaterial( Kr,       fractionmass = 0.986 );
  Kr7CH4->AddMaterial( metane,   fractionmass = 0.014 );

  // 0.5*(95% Xe + 5% CH4)+0.5*(93% Ar + 7% CH4), STP ; NIM A214 (1983) 261-268

  density = 3.655*mg/cm3;
  G4Material* XeArCH4 = new G4Material(name="XeArCH4"  , density,
                                                  ncomponents=2);
  XeArCH4->AddMaterial( Xe5CH4,       fractionmass = 0.766 );
  XeArCH4->AddMaterial( Ar7CH4,   fractionmass = 0.234 );

  // Silicon as detector material

  density = 2.330*g/cm3;
  a = 28.09*g/mole;
  new G4Material(name="Si", z=14., a, density);




  /*
  G4Material* ma;
  ma  = new G4Material("H",     1.,  1.0*g/mole, 1.*g/cm3);
  ma  = new G4Material("D",     1.,  2.0*g/mole, 1.*g/cm3);
  ma  = new G4Material("Li",    3.,  6.941*g/mole, 1.*g/cm3);
  ma  = new G4Material("Be",    4.,  9.01*g/mole, 1.848*g/cm3);
  ma  = new G4Material("C",     6.,  12.00*g/mole, 2.0*g/cm3);
        ma  = new G4Material("Graphite",6., 12.00*g/mole, 2.265*g/cm3 );
  ma->SetChemicalFormula("Graphite");
  ma  = new G4Material("Al",    13.,  26.98*g/mole,  2.7 *g/cm3);
        ma  = new G4Material("Si",    14.,  29.055*g/mole, 2.33*g/cm3);
        ma  = new G4Material("LAr",   18.,  39.95*g/mole,  1.393*g/cm3);
        ma  = new G4Material("Zr",    40.,  91.224*g/mole, 4.0*g/cm3);
        ma  = new G4Material("LXe",   54., 131.29*g/mole,  3.02*g/cm3);
        ma  = new G4Material("Fe",    26.,  55.85*g/mole,  7.87*g/cm3);
        ma  = new G4Material("Ni",    29.,  58.6934*g/mole,  8.00*g/cm3);
        ma  = new G4Material("Cu",    29.,  63.55*g/mole,  8.96*g/cm3);
        ma  = new G4Material("Au",    79., 196.97*g/mole, 19.32*g/cm3);
        ma  = new G4Material("Ta",    73., 180.9479*g/mole, 16.67*g/cm3);
        ma  = new G4Material("W",     74., 183.85*g/mole, 19.30*g/cm3);
        ma  = new G4Material("Pb",    82., 207.19*g/mole, 11.35*g/cm3);
        ma  = new G4Material("Bi",    83., 208.98*g/mole, 12.*g/cm3);
        ma  = new G4Material("U",     92., 238.03*g/mole, 18.95*g/cm3);

  G4Element*   H  = new G4Element ("Hydrogen", "H",   1. ,  1.01*g/mole);
  G4Element*   N  = new G4Element ("Nitrigen", "N",   7. , 14.00*g/mole);
  G4Element*   O  = new G4Element ("Oxygen"  , "O",   8. , 16.00*g/mole);
  G4Element*   C  = new G4Element ("Carbon"  , "C",   6. , 12.00*g/mole);
  G4Element*  Cs  = new G4Element ("Cesium"  , "Cs", 55. , 132.905*g/mole);
  G4Element*   I  = new G4Element ("Iodide"  , "I",  53. , 126.9044*g/mole);

  ma = new G4Material("O2", 8., 16.00*g/mole, 1.1*g/cm3);
  ma->SetChemicalFormula("O_2");
  ma = new G4Material ("Water" , 1.*g/cm3, 2);
  ma->AddElement(H,2);
  ma->AddElement(O,1);
  ma->SetChemicalFormula("H_2O");
  ma = new G4Material ("Ethane" , 0.4241*g/cm3, 2);
  ma->AddElement(H,6);
  ma->AddElement(C,2);
  ma->SetChemicalFormula("C_2H_6");
  ma = new G4Material ("CsI" , 4.53*g/cm3, 2);
  ma->AddElement(Cs,1);
  ma->AddElement(I,1);
  ma->SetChemicalFormula("CsI");
  ma = new G4Material("Air"  , 1.290*mg/cm3, 2);
        // use fraction in mass
  ma->AddElement(N, 0.7);
  ma->AddElement(O, 0.3);
  */



}


G4Material* Materials::GetMaterial(const G4String& name)
{ 

  //  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

        G4Material* ma = G4Material::GetMaterial(name);
        
  G4cout << "Material is selected: " << ma->GetName() << G4endl;
  return ma;
}        

  






