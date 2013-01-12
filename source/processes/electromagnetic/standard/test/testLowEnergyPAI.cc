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
//  Test routine based on G4PAIdNdxTest for low energy extension of PAI model
//
// History:
//
// 18.09.09, V. Grichine implementation based on G4PAIdNdxTest

#include "G4ios.hh"
#include <fstream>
#include <cmath>
#include "globals.hh"
#include "Randomize.hh"

#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4SandiaTable.hh"

// #include "G4PAIonisation.hh"
#include "G4PAIxSection.hh"



G4double GetXscSimple(G4Element* el, G4double Tkin)
{
  G4int    shell, NbOfShells, Nshell; 
  G4double tmpxsc, lnt, t, B, xsc  = 0.;
  G4double minust, plust, tplus, pert;

  G4double ryd  = 13.6*eV;  // Rydberg number
  G4double ryd2 = ryd*ryd;

  G4double cof = 4*pi*Bohr_radius*Bohr_radius*ryd2;

  NbOfShells = el->GetNbOfAtomicShells();


  for ( shell = 0; shell < NbOfShells; shell++ )
  { 
    B      = el->GetAtomicShell(shell);            
    Nshell = el->GetNbOfShellElectrons(shell); 
 
    t      = Tkin/B; 
    lnt    = std::log(t);         
    tplus  = t + 1.;
    pert   = 1./t;
    minust = 1. - pert;
    plust  = 1. + pert;

    tmpxsc = 0.5*minust*plust*lnt + minust - lnt/tplus;
    tmpxsc *= Nshell/B/B/(tplus+1); 
    xsc += tmpxsc;
  }
  xsc *= cof;
  return xsc;
}








int main()
{
   std::ofstream outFile("PAIdNdx.out", std::ios::out );
   outFile.setf( std::ios::scientific, std::ios::floatfield );

   std::ofstream fileOut("PAICrenkovPlasmon.out", std::ios::out );
   fileOut.setf( std::ios::scientific, std::ios::floatfield );

   //  std::ifstream fileRead("exp.dat", std::ios::out );
   //  fileRead.setf( std::ios::scientific, std::ios::floatfield );

   std::ofstream fileWrite("exp.dat", std::ios::out );
   fileWrite.setf( std::ios::scientific, std::ios::floatfield );

   std::ofstream fileWrite1("mprrpai.dat", std::ios::out );
   fileWrite1.setf( std::ios::scientific, std::ios::floatfield );

// Create materials  
   

   G4int iz , n,  nel, ncomponents, natoms;
  G4double a, z, ez, density , temperature, pressure, fractionmass;
  G4State state;
  G4String name, symbol;

  // G4Element*   elH = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);

  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", ez=7., a);

  a = 16.00*g/mole;
  // G4Element* elO = new G4Element(name="Oxigen", symbol="O", ez=8., a);

  a = 12.01*g/mole;
  G4Element* elC = new G4Element(name="Carbon",symbol="C", ez=6., a);

  a = 55.85*g/mole;
  G4Element* elFe = new G4Element(name="Iron",symbol="Fe", ez=26., a);

  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxygen",symbol="O", ez=8., a);

  a = 1.01*g/mole;
  G4Isotope* ih1 = new G4Isotope("Hydrogen",iz=1,n=1,a);

  a = 2.01*g/mole;
  G4Isotope* ih2 = new G4Isotope("Deuterium",iz=1,n=2,a);

  G4Element* elH = new G4Element(name="Hydrogen",symbol="H",2);
  elH->AddIsotope(ih1,.999);
  elH->AddIsotope(ih2,.001);

  a = 39.948*g/mole;
  G4Element* elAr = new G4Element(name="Argon", symbol="Ar", z=18., a);

  a = 131.29*g/mole;
  G4Element* elXe = new G4Element(name="Xenon", symbol="Xe", z=54., a);
  
  a = 19.00*g/mole;
  G4Element* elF  = new G4Element(name="Fluorine", symbol="F", z=9., a);

  a = 69.723*g/mole;
  G4Element* elGa  = new G4Element(name="Ga", symbol="Ga", z=31., a);

  a = 74.9216*g/mole;
  G4Element* elAs  = new G4Element(name="As", symbol="As", z=33., a);

 
// G4Isotope::DumpInfo();
// G4Element::DumpInfo();
// G4Material::DumpInfo();

 
  // Materials for PAI model test in low energy region


  a = 9.012*g/mole;
  density = 1.848*g/cm3;
  G4Material* Be = new G4Material(name="Beryllium", z=4. , a, density);

  density = 1.390*g/cm3;
  a = 39.95*g/mole;
  G4Material* lAr = new G4Material(name="liquidArgon", z=18., a, density);

  density = 19.32*g/cm3;
  a =196.97*g/mole;
  G4Material* Au = new G4Material(name="Gold"   , z=79., a, density);

  // Carbon dioxide

  density = 1.977*mg/cm3;
  G4Material* CO2 = new G4Material(name="CO2", density, nel=2,
				       kStateGas,273.15*kelvin,1.*atmosphere);
  CO2->AddElement(elC,1);
  CO2->AddElement(elO,2);

  density = 1.290*mg/cm3;  // old air from elements
  G4Material* air = new G4Material(name="air"  , density, ncomponents=2);
  air->AddElement(elN, fractionmass=0.7);
  air->AddElement(elO, fractionmass=0.3);


  density = 1.25053*mg/cm3;       // STP
  a = 14.01*g/mole;       // get atomic weight !!!
  //  a = 28.016*g/mole;
  G4Material* newN2  = new G4Material(name="newN2", z= 7.,a,density);

  density = 1.25053*mg/cm3;       // STP
  G4Material* anotherN2 = new G4Material(name="anotherN2", density,ncomponents=2);
  anotherN2->AddElement(elN, 1);
  anotherN2->AddElement(elN, 1);

  density = 1.000*g/cm3;
  G4Material* H2O = new G4Material(name="Water", density, ncomponents=2);
  H2O->AddElement(elH, natoms=2);
  H2O->AddElement(elO, natoms=1);



  a = 26.98159*g/mole;
  density = 2.7*g/cm3;
  G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);

  // Silicon as detector material

  density = 2.330*g/cm3;
  a = 28.0855*g/mole;
  G4Material* Si = new G4Material(name="Silicon", z=14., a, density);

 
  density = 7.870*g/cm3;
  a = 55.85*g/mole;
  G4Material* Fe = new G4Material(name="Iron"   , z=26., a, density);

  density = 8.960*g/cm3;
  a = 63.55*g/mole;
  G4Material* Cu = new G4Material(name="Copper"   , z=29., a, density);

  density = 11.35*g/cm3;
  a = 207.19*g/mole;
  G4Material* Pb = new G4Material(name="Lead"     , z=82., a, density);

  // Polypropelene

  G4Material* CH2 = new G4Material ("Polypropelene" , 0.91*g/cm3, 2);
  CH2->AddElement(elH,2);
  CH2->AddElement(elC,1);

  // Kapton (polyimide)

  density = 1.39*g/cm3;
  G4Material* Kapton = new G4Material(name="Kapton", density, nel=3);
  Kapton->AddElement(elO,2);
  Kapton->AddElement(elC,5);
  Kapton->AddElement(elH,4);


  // Germanium as detector material

  density = 5.323*g/cm3;
  a = 72.59*g/mole;
  G4Material* Ge = new G4Material(name="Ge", z=32., a, density);

  // GaAs detectors

  density = 5.32*g/cm3;
  G4Material* GaAs = new G4Material(name="GaAs",density, nel=2);
  GaAs->AddElement(elGa,1);
  GaAs->AddElement(elAs,1);

  // Diamond detectors

  density = 3.5*g/cm3;
  G4Material* Diamond = new G4Material(name="Diamond",density, nel=1);
  Diamond->AddElement(elC,1);

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


  // Dry air (average composition)


  density = 0.0899*mg/cm3;       // STP
  G4Material* matH2 = new G4Material(name="matH2"  , density, ncomponents=1);
  matH2->AddElement(elH, 2);



  density = 1.25053*mg/cm3;       // STP
  G4Material* Nitrogen = new G4Material(name="N2"  , density, ncomponents=1);
  Nitrogen->AddElement(elN, 2);

  density = 1.4289*mg/cm3;       // STP
  G4Material* Oxygen = new G4Material(name="O2"  , density, ncomponents=1);
  Oxygen->AddElement(elO, 2);

  density = 1.7836*mg/cm3;       // STP
  G4Material* Argon = new G4Material(name="Argon"  , density, ncomponents=1);
  Argon->AddElement(elAr, 1);

  density = 1.2928*mg/cm3;       // STP
  G4Material* Air = new G4Material(name="Air"  , density, ncomponents=3);
  Air->AddMaterial( Nitrogen, fractionmass = 0.7557 );
  Air->AddMaterial( Oxygen,   fractionmass = 0.2315 );
  Air->AddMaterial( Argon,    fractionmass = 0.0128 );

  // Xenon as detector gas, STP

  density = 5.858*mg/cm3;
  a = 131.29*g/mole;
  G4Material* Xe  = new G4Material(name="Xenon",z=54., a, density );

  // Helium as detector gas, STP

  density = 0.178*mg/cm3;
  a = 4.0026*g/mole;
  G4Material* He  = new G4Material(name="He",z=2., a, density );


  // Krypton as detector gas, STP

  density = 3.700*mg/cm3;
  a = 83.80*g/mole;
  G4Material* Kr  = new G4Material(name="Kr",z=36., a, density );

  

  // Neon as detector gas, STP

  density = 0.900*mg/cm3;
  a = 20.179*g/mole;
  G4Material* Ne  = new G4Material(name="Ne",z=10., a, density );

  // Carbone dioxide, CO2 STP

  density = 1.977*mg/cm3;
  G4Material* CarbonDioxide = new G4Material(name="CO2", density, nel=2);
  CarbonDioxide->AddElement(elC,1);
  CarbonDioxide->AddElement(elO,2);

   

  // Metane, STP

  density = 0.7174*mg/cm3;
  G4Material* metane = new G4Material(name="CH4",density,nel=2);
  metane->AddElement(elC,1);
  metane->AddElement(elH,4);

  // Propane, STP

  density = 2.005*mg/cm3;
  G4Material* propane = new G4Material(name="C3H8",density,nel=2);
  propane->AddElement(elC,3);
  propane->AddElement(elH,8);

  // iso-Butane (methylpropane), STP

  density = 2.67*mg/cm3;
  G4Material* isobutane = new G4Material(name="isoC4H10",density,nel=2);
  isobutane->AddElement(elC,4);
  isobutane->AddElement(elH,10);

  // 87.5% Xe + 7.5% CH4 + 5% C3H8, 20 C, 1 atm 

  density = 4.9196*mg/cm3;

  G4Material* XeCH4C3H8 = new G4Material(name="XeCH4C3H8"  , density, 
                                                             ncomponents=3);
  XeCH4C3H8->AddMaterial( Xe,       fractionmass = 0.971 );
  XeCH4C3H8->AddMaterial( metane,   fractionmass = 0.010 );
  XeCH4C3H8->AddMaterial( propane,  fractionmass = 0.019 );

  // Propane in MWPC, 2 atm, 20 C

  //  density = 3.758*mg/cm3;
  density = 3.736*mg/cm3;
  G4Material* propaneDet = new G4Material(name="detC3H8",density,nel=2);
  propaneDet->AddElement(elC,3);
  propaneDet->AddElement(elH,8);

  

  // 90% Ne + 10% CO2, STP

  density = 1.0077*mg/cm3;      
  G4Material* Ne10CO2 = new G4Material(name="Ne10CO2"  , density, 

                                                             ncomponents=2);
  Ne10CO2->AddMaterial( Ne,           fractionmass = 0.8038 );
  Ne10CO2->AddMaterial( CarbonDioxide,   fractionmass = 0.1962 );



  // 80% Ar + 20% CO2, STP

  density = 1.8223*mg/cm3;      
  G4Material* Ar20CO2 = new G4Material(name="Ar20CO2"  , density, 
                                                             ncomponents=2);
  Ar20CO2->AddMaterial( Argon,           fractionmass = 0.783 );
  Ar20CO2->AddMaterial( CarbonDioxide,   fractionmass = 0.217 );

  // 93% Ar + 7% CH4, STP

  density = 1.709*mg/cm3;      
  G4Material* Ar7CH4 = new G4Material(name="Ar7CH4"  , density, 
                                                             ncomponents=2);
  Ar7CH4->AddMaterial( Argon,    fractionmass = 0.971 );
  Ar7CH4->AddMaterial( metane,   fractionmass = 0.029 );

  // 80% Xe + 20% CO2, STP

  density = 5.0818*mg/cm3;      
  G4Material* Xe20CO2 = new G4Material(name="Xe20CO2"  , density, 
                                                             ncomponents=2);
  Xe20CO2->AddMaterial( Xe,              fractionmass = 0.922 );
  Xe20CO2->AddMaterial( CarbonDioxide,   fractionmass = 0.078 );

  // 80% Kr + 20% CO2, STP

  density = 3.601*mg/cm3;      
  G4Material* Kr20CO2 = new G4Material(name="Kr20CO2"  , density, 
                                                             ncomponents=2);
  Kr20CO2->AddMaterial( Kr,              fractionmass = 0.89 );
  Kr20CO2->AddMaterial( CarbonDioxide,   fractionmass = 0.11 );

  // 80% He + 20% CO2, STP

  density = 0.5378*mg/cm3;      
  G4Material* He20CO2 = new G4Material(name="He20CO2"  , density, 
                                                             ncomponents=2);
  He20CO2->AddMaterial( He,              fractionmass = 0.265 );
  He20CO2->AddMaterial( CarbonDioxide,   fractionmass = 0.735 );

 


  //  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  G4int i, iMax, j, jMax, k, numOfMaterials, iSan, nbOfElements, sanIndex, row;
  G4double maxEnergyTransfer, kineticEnergy, dEdx, matNat, matA = 0.;
  G4double tau, gamma, bg2, bg, beta2, beBohr, rateMass, Tmax, Tmin, Tkin, ratioTkin, enFermi;

   const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
   numOfMaterials = theMaterialTable->size();




  // G4String testName = "Ne10CO2";
  // G4String testName = "t2kGasMixture";
  // G4String testName = "Ne10CO2T293";
  // G4String testName = "Ne857CO295N2T292";


   // protons

   //  G4String testName = "He";
   // G4String testName = "N2";
   //   G4String testName = "O2";
   // G4String testName = "Ne";
  // G4String testName = "Argon";
  // G4String testName = "Kr";
  // G4String testName = "Xenon";
   // G4String testName = "CO2";


   // electrons mass dE/dx MeV*cm2/g

   // G4String testName = "Water";
   // G4String testName = "matH2";
   // G4String testName = "N2";
   // G4String testName = "O2";
   // G4String testName = "CO2";

   // electrons dE/dx eV/A

   //  G4String testName = "Aluminium";
   // G4String testName = "Gold";
   //  G4String testName = "Copper";
    G4String testName = "Silicon";

   // G4String testName = "Ne";


  // G4cout<<"Enter material name for test : "<<std::flush;
  //  G4cin>>testName;



  for( k = 0; k < numOfMaterials; k++ )
  {
    if((*theMaterialTable)[k]->GetName() != testName) continue;

    // outFile << "Material : " <<(*theMaterialTable)[k]->GetName() << G4endl;
     G4cout << "Material : " <<(*theMaterialTable)[k]->GetName() << G4endl;

     nbOfElements = (*theMaterialTable)[k]->GetNumberOfElements();

     G4cout<<"Sandia cof according old PAI stuff"<<G4endl<<G4endl;
     // outFile<<"Sandia cof according old PAI stuff"<<G4endl<<G4endl;

     G4int* thisMaterialZ = new G4int[nbOfElements];

     matA  = 0.;

     for( iSan = 0; iSan < nbOfElements; iSan++ )
     {
        thisMaterialZ[iSan] = (G4int)(*theMaterialTable)[k]->
                                      GetElement(iSan)->GetZ();
        matA  += (*theMaterialTable)[k]->GetElement(iSan)->GetA(); 
     }
     G4SandiaTable sandia(k);
     // sandia.SetVerbose(1);
     sanIndex = sandia.SandiaIntervals(thisMaterialZ,nbOfElements);    
     sanIndex = sandia.SandiaMixing( thisMaterialZ ,
                             (*theMaterialTable)[k]->GetFractionVector(),
				     nbOfElements,sanIndex);

     for(row = 0; row < sanIndex-1; row++ )
     {
       G4cout<<row+1<<"\t"<<sandia.GetPhotoAbsorpCof(row+1,0)/keV;
       // outFile<<row+1<<"  "<<sandia.GetPhotoAbsorpCof(row+1,0)/keV;

       for( iSan = 1; iSan < 5; iSan++ )
       {
         G4cout<<"\t"<<sandia.GetPhotoAbsorpCof(row+1,iSan);
	 // *(*theMaterialTable)[k]->GetDensity();

         // outFile<<"  "<<sandia.GetPhotoAbsorpCof(row+1,iSan);
	 // *(*theMaterialTable)[k]->GetDensity();
       }
       G4cout<<G4endl;
       // outFile<<G4endl;
     }
     G4cout<<G4endl;
     // outFile<<G4endl;


     // outFile<<G4endl;
     // maxEnergyTransfer = 100*keV;
     maxEnergyTransfer = 10000*keV;
     gamma             = 4.0;
     bg2               = gamma*gamma - 1;

     G4PAIxSection testPAI( k, maxEnergyTransfer, bg2);

     G4cout<<"Interval no."<<"\t"<<"Energy interval"<<G4endl<<G4endl;
     // outFile<<"Interval no."<<"\t"<<"Energy interval"<<G4endl<<G4endl;

     for( j = 1; j <= testPAI.GetIntervalNumber(); j++ )
     {
       G4cout<<j<<"\t\t"<<testPAI.GetEnergyInterval(j)/keV<<G4endl;
       // outFile<<j<<"\t\t"<<testPAI.GetEnergyInterval(j)/keV<<G4endl;
     }
     G4cout<<G4endl;
     // outFile<<G4endl;

     G4cout << "Actual spline size = "<<testPAI.GetSplineSize()<<G4endl;
     G4cout <<"Normalization Cof = "<<testPAI.GetNormalizationCof()<<G4endl;
     G4cout << G4endl;

     // outFile<<"Actual spline size = "<<testPAI.GetSplineSize()<<G4endl;
     // outFile<<"Normalization Cof = "<<testPAI.GetNormalizationCof()<<G4endl;
     // outFile<<G4endl;


     Tmin     = sandia.GetPhotoAbsorpCof(1,0);  // 0.02*keV;
     G4cout<<"Tmin = "<<Tmin/eV<<" eV"<<G4endl;
   
     G4cout 
       //     <<"Tkin, keV"<<"\t"
            << "bg"<<"\t\t"
       //   <<"Max E transfer, kev"<<"\t"
            << "<dN/dxC>, 1/cm"<<"\t"
            << "<dN/dxMM>, 1/cm"<<"\t"
            << "<dN/dxP>, 1/cm"<<"\t"
       //    << "<dN/dxC+dN/dxP>"<<"\t"
            <<"<dN/dx>, 1/cm"<<G4endl<<G4endl;
   
     /*
     outFile
       //     <<"Tkin, keV"<<"\t"
            <<"gamma"<<"\t\t"
       //   <<"Max E transfer, kev"<<"\t"
            <<"<dN/dxC>, 1/cm"<<"\t"
            << "<dN/dxP>, 1/cm"<<"\t"
            <<"<dN/dxC+dN/dxP>"<<"\t"
            <<"<dN/dx>, 1/cm"<<G4endl<<G4endl;
     */
     //   G4PAIxSection testPAIproton(k,maxEnergyTransfer);





     // kineticEnergy = 10.0*keV; // 100.*GeV;    // 10.0*keV;  // 110*MeV; // for proton

     // kineticEnergy = 4*GeV; // for electrons

     kineticEnergy = 15*eV*proton_mass_c2/electron_mass_c2;

     // kineticEnergy = 5.*eV;

     //     for(j=1;j<testPAIproton.GetNumberOfGammas();j++)

     // jMax = 70; // 70;
     jMax = 200; // 70;
     iMax = 2;

     outFile<<jMax<<G4endl;

     if( testName == "He" || testName == "Ne" || 
         testName == "Argon" || testName == "Kr" || testName == "Xe" )
     {
       matNat = Avogadro*(*theMaterialTable)[k]->GetDensity()/matA;
     }
     else matNat = (2.7e19)/cm3;

     // G4cout<<"dNat/dV = "<<matNat*cm3<<G4endl;
     G4cout<<"density = "<<(*theMaterialTable)[k]->GetDensity()*cm3/g<<G4endl;

     for( j = 0; j < jMax; j++ )
     {
       tau      = kineticEnergy/proton_mass_c2;
       // tau      = kineticEnergy/electron_mass_c2;
       gamma    = tau +1.0;
       bg2      = tau*(tau + 2.0);
       bg = std::sqrt(bg2);
       beta2    = bg2/(gamma*gamma);
       beBohr = std::sqrt(beta2)/fine_structure_const;

       // ratioTkin = tau*amu_c2/keV;

       // G4cout<<"bg = "<<bg<<";  beBohr = "<<beBohr<<G4endl<<G4endl;

       rateMass = electron_mass_c2/proton_mass_c2;
       // ratioTkin = kineticEnergy*rateMass/eV;
       enFermi = 12.5*eV;
       ratioTkin = (kineticEnergy*rateMass-enFermi)/eV;
       
       // Tmax = 2.0*electron_mass_c2*bg2/(1.0+2.0*gamma*rateMass+rateMass*rateMass);
       Tmax = 0.5*kineticEnergy*rateMass;

       Tkin = maxEnergyTransfer;

       if ( maxEnergyTransfer > Tmax)         
       {
          Tkin = Tmax;
       }
       if ( Tmax <= Tmin + 0.5*eV )         
       {
          Tkin = Tmin + 0.5*eV;
       }
       // G4cout<<"Tmax = "<<Tmax<<"; Tkin = "<<Tkin<<G4endl;

       G4PAIxSection testPAIproton(k,Tkin,bg2);

       dEdx = testPAIproton.GetIntegralPAIdEdx(1);

       // dEdx /= matNat;
       // dEdx *= (1.0e15)/(cm2*eV);

       // dEdx /= (*theMaterialTable)[k]->GetDensity();
       // dEdx *= g/cm2/MeV;

       dEdx *= angstrom/eV;

       
       G4cout<<"ratioTkin = "<<ratioTkin<<", dEdx = "<<dEdx<<", Tmax = "<<Tmax/keV<<"; Tkin = "<<Tkin/keV<<G4endl;
       outFile<<ratioTkin<<"\t"<<dEdx<<G4endl;


       /*       
       G4cout  
         //      << kineticEnergy/keV<<"\t\t"
         //      << gamma << "\t\t"
               << bg << "\t\t"
	 //    << Tkin/keV<<"\t\t"     
               << testPAIproton.GetIntegralCerenkov(1)*cm << "\t"
               << testPAIproton.GetIntegralMM(1)*cm << "\t"
               << testPAIproton.GetIntegralPlasmon(1)*cm << "\t"
               << testPAIproton.GetIntegralResonance(1)*cm << "\t"
	 // << testPAIproton.GetIntegralCerenkov(1)*cm +
	 //      testPAIproton.GetIntegralPlasmon(1)*cm << "\t"
	 //      << FitALICE(bg) << "\t"
	 //      << FitBichsel(bg) << "\t"
               << testPAIproton.GetIntegralPAIxSection(1)*cm << "\t\t" 
          << G4endl;


       	        
       outFile 
         //      << kineticEnergy/keV<<"\t"
	 //       << gamma << "\t"
               << bg << "\t\t"
         //      << Tkin/keV<<"\t"
               << testPAIproton.GetIntegralCerenkov(1)*cm << "\t"
               << testPAIproton.GetIntegralMM(1)*cm << "\t"
              << testPAIproton.GetIntegralPlasmon(1)*cm << "\t"
               << testPAIproton.GetIntegralResonance(1)*cm << "\t"
         //      << testPAIproton.GetIntegralCerenkov(1)*cm +
         //      testPAIproton.GetIntegralPlasmon(1)*cm << "\t"
	 //      << FitALICE(bg) << "\t"
	 //      << FitBichsel(bg) << "\t"
               << testPAIproton.GetIntegralPAIxSection(1)*cm << "\t" 
               << G4endl;

       //   outFile<<testPAIproton.GetLorentzFactor(j)<<"\t"
       //          <<maxEnergyTransfer/keV<<"\t\t"
       //          <<testPAIproton.GetPAItable(0,j)*cm/keV<<"\t\t"
       //  	      <<testPAIproton.GetPAItable(1,j)*cm<<"\t\t"<<G4endl;

       */  
          
       /*
       outFile<<testPAIproton.GetSplineSize()-1<<G4endl;
       
       // for( i = 1; i < testPAIproton.GetSplineSize(); i++)
       for( i = 1; i < iMax; i++)
       {
       outFile 
               << testPAIproton.GetSplineEnergy(i)/keV       << "\t"
	 //  << testPAIproton.GetIntegralCerenkov(i)*cm    << "\t"
	 //     << testPAIproton.GetIntegralMM(i)*cm    << "\t"
	 //     << testPAIproton.GetIntegralPlasmon(i)*cm     << "\t"
         //      << testPAIproton.GetIntegralResonance(i)*cm   << "\t"
               << testPAIproton.GetIntegralPAIxSection(i)*cm << "\t" 
               << G4endl;

       G4cout 
               << testPAIproton.GetSplineEnergy(i)/keV       << "\t"
               << testPAIproton.GetIntegralCerenkov(i)*cm    << "\t"
               << testPAIproton.GetIntegralMM(i)*cm    << "\t"
               << testPAIproton.GetIntegralPlasmon(i)*cm     << "\t"
               << testPAIproton.GetIntegralResonance(i)*cm   << "\t"
               << testPAIproton.GetIntegralPAIxSection(i)*cm << "\t" 
               << G4endl;

       }
       */       
      

       
       /*
       G4double position, transfer, lambda, range, r2cer=0., r2res=0., r2ruth=0., r2tot=0.;
       G4int nCer = 0, nRes = 0, nRuth = 0, nTot = 0;
       G4double rBin[100], rDistr[100], rTemp, rTemp2, sumDistr = 0., rSum = 0;
       G4double ionBin[100], ionDistr[100], ionMean, ionRand, F = 0.19, ionSum=0., ionSigma;


       for( i = 0; i < 100; i++)
       {
         ionBin[i] = i*1.;
         ionDistr[i] = 0.;
         rBin[i] = i/200.;
         rDistr[i] = 0.;
       }
       for( i = 0; i < 10000; i++)
       {
	 
         position = testPAIproton.GetIntegralPAIxSection(1)*G4UniformRand();

	 if( position < testPAIproton.GetIntegralCerenkov(1) )
	 {
           transfer = testPAIproton.GetCerenkovEnergyTransfer();
           lambda   = testPAIproton.GetPhotonRange(transfer);
           range    = testPAIproton.GetElectronRange(transfer);
           r2cer   += 0.67*(lambda+range)*(lambda+range);
           r2tot   += 0.67*(lambda+range)*(lambda+range);
           rTemp2   = 0.67*(lambda+range)*(lambda+range);
           nCer++;           
	 }
         else if( position < (testPAIproton.GetIntegralCerenkov(1)+
                              testPAIproton.GetIntegralResonance(1) ))
	 {
           transfer = testPAIproton.GetResonanceEnergyTransfer();          
           range    = testPAIproton.GetElectronRange(transfer);
           r2res   += 0.67*range*range;
           r2tot   += 0.67*range*range;           
           rTemp2   = 0.67*range*range;           
           nRes++;           
	 }
	 else
	 {
           transfer = testPAIproton.GetRutherfordEnergyTransfer();          
           range    = testPAIproton.GetElectronRange(transfer);
           r2ruth  += range*range;
           r2tot   += range*range;           
           rTemp2   = range*range;           
           nRuth++;           
	 }
         nTot++;

	 rTemp = std::sqrt(rTemp2);
         rSum += rTemp;

	 // rTemp = rTemp2;

	 for( j = 0; j < 100; j++ )
	 {
           if( rTemp <= rBin[j] )
	   {
             rDistr[j] += 1.;
             break;
	   }
	 }
	 

         transfer = testPAIproton.GetEnergyTransfer();          
	 ionMean  = GetIonisation(transfer);
         ionSigma = std::sqrt(F*ionMean);

         // ionRand  = G4RandGauss::shoot(ionMean, ionSigma);
         ionRand  = ionMean;

         if( ionRand < 0.) ionRand =0.;	

         ionSum += ionRand; 
         nTot++;
 
	 for( j = 0; j < 100; j++ )
	 {
           if( ionRand <= ionBin[j] )
	   {
             ionDistr[j] += 1.;
             break;
	   }
	 }
       }
       
       if(nCer >0)  r2cer  /= nCer;
       if(nRes >0)  r2res  /= nRes;
       if(nRuth >0) r2ruth /= nRuth;
                    r2tot /= nTot;
                    rSum  /= nTot; 
       G4cout<<"nCer = "<<nCer<<"; nRes = "<<nRes<<"; nRuth = "<<nRuth<<G4endl;
       G4cout<<"sum of n = "<<nCer+nRes+nRuth<<"; nTot = "<<nTot<<G4endl;
       G4cout<<"rCer = "<<std::sqrt(r2cer)<<" mm; rRes = "<<std::sqrt(r2res)<<" mm"<<G4endl;
       G4cout<<"rRuth = "<<std::sqrt(r2ruth)<<" mm; rTot = "<<std::sqrt(r2tot)<<" mm"<<G4endl;
       G4cout<<"rSum = "<<rSum<<" mm; "<<G4endl;
       
                    ionSum  /= nTot; 
       G4cout<<"ionSum = "<<ionSum<<" electrons"<<G4endl;

       outFile<<100<<G4endl;

       for( j = 0; j < 100; j++ )
       {
         // outFile<<rBin[j]<<"\t"<<rDistr[j]<<G4endl;
         outFile<<ionBin[j]<<"\t"<<ionDistr[j]<<G4endl;
         sumDistr += rDistr[j];
       }
       G4cout<<"sumDistr = "<<sumDistr<<G4endl;
       */
       


       kineticEnergy *= 1.07;        // was 1.4; 1.5;
     }

     G4cout<<G4endl;
     // outFile<<G4endl;
  }



  return EXIT_SUCCESS;

}










