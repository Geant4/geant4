// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PAIxSectionTest.cc,v 1.4 1999-12-15 14:51:54 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//
//  
//
//  Test routine for G4PAIxSection class code
//
// History:
//
// 21.10.99, V. Grichine implementation based on G4PAIonisationTest

#include "G4ios.hh"
#include "g4std/fstream"
#include <math.h>
#include "globals.hh"
#include "Randomize.hh"

#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4SandiaTable.hh"

#include "G4PAIonisation.hh"
#include "G4PAIxSection.hh"

int main()
{
   G4std::ofstream outFile("PAIdEdx.out", G4std::ios::out ) ;
   outFile.setf( G4std::ios::scientific, G4std::ios::floatfield );

   G4std::ofstream fileOut("PAIdistribution.out", G4std::ios::out ) ;
   fileOut.setf( G4std::ios::scientific, G4std::ios::floatfield );

   //  G4std::ifstream fileRead("exp.dat", G4std::ios::out ) ;
   //  fileRead.setf( G4std::ios::scientific, G4std::ios::floatfield );

   G4std::ofstream fileWrite("exp.dat", G4std::ios::out ) ;
   fileWrite.setf( G4std::ios::scientific, G4std::ios::floatfield );

   G4std::ofstream fileWrite1("mprrpai.dat", G4std::ios::out ) ;
   fileWrite1.setf( G4std::ios::scientific, G4std::ios::floatfield );

// Create materials  
   

  G4int iz , n,  nel, ncomponents ;
  G4double a, z, ez, density , temperature, pressure, fractionmass ;
  G4State state ;
  G4String name, symbol ;

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

  /* ***************************************************************

  a = 9.012*g/mole;
  density = 1.848*g/cm3;
  G4Material* Be = new G4Material(name="Beryllium", z=4. , a, density);

  a = 26.98*g/mole;
  density = 2.7*g/cm3;
  G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);

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

  // Carbon dioxide

  density = 1.977*mg/cm3;
  G4Material* CO2 = new G4Material(name="CO2", density, nel=2,
				       kStateGas,273.15*kelvin,1.*atmosphere);
  CO2->AddElement(elC,1);
  CO2->AddElement(elO,2);

  density = 1.290*mg/cm3;  // old air from elements
  G4Material* air = new G4Material(name="air"  , density, ncomponents=2);
  Air->AddElement(elN, fractionmass=0.7);
  Air->AddElement(elO, fractionmass=0.3);


  density = 1.25053*mg/cm3 ;       // STP
  a = 14.01*g/mole ;       // get atomic weight !!!
  //  a = 28.016*g/mole;
  G4Material* newN2  = new G4Material(name="newN2", z= 7.,a,density) ;

  density = 1.25053*mg/cm3 ;       // STP
  G4Material* anotherN2 = new G4Material(name="anotherN2", density,ncomponents=2);
  anotherN2->AddElement(elN, 1);
  anotherN2->AddElement(elN, 1);

  density = 1.000*g/cm3;
  G4Material* H2O = new G4Material(name="Water", density, ncomponents=2);
  H2O->AddElement(elH, natoms=2);
  H2O->AddElement(elO, natoms=1);

  ***************************************************** */



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

  // Silicon as detector material

  density = 2.330*g/cm3;
  a = 28.09*g/mole;
  G4Material* Si = new G4Material(name="Silicon", z=14., a, density);

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


  density = 1.25053*mg/cm3 ;       // STP
  G4Material* Nitrogen = new G4Material(name="N2"  , density, ncomponents=1);
  Nitrogen->AddElement(elN, 2);

  density = 1.4289*mg/cm3 ;       // STP
  G4Material* Oxygen = new G4Material(name="O2"  , density, ncomponents=1);
  Oxygen->AddElement(elO, 2);

  density = 1.7836*mg/cm3 ;       // STP
  G4Material* Argon = new G4Material(name="Argon"  , density, ncomponents=1);
  Argon->AddElement(elAr, 1);

  density = 1.2928*mg/cm3 ;       // STP
  G4Material* Air = new G4Material(name="Air"  , density, ncomponents=3);
  Air->AddMaterial( Nitrogen, fractionmass = 0.7557 ) ;
  Air->AddMaterial( Oxygen,   fractionmass = 0.2315 ) ;
  Air->AddMaterial( Argon,    fractionmass = 0.0128 ) ;

  // Xenon as detector gas, STP

  density = 5.858*mg/cm3 ;
  a = 131.29*g/mole ;
  G4Material* Xe  = new G4Material(name="Xenon",z=54., a, density );

  // Helium as detector gas, STP

  density = 0.178*mg/cm3 ;
  a = 4.0026*g/mole ;
  G4Material* He  = new G4Material(name="He",z=2., a, density );

  // Neon as detector gas, STP

  density = 0.900*mg/cm3 ;
  a = 20.179*g/mole ;
  G4Material* Ne  = new G4Material(name="Ne",z=10., a, density );

  // Krypton as detector gas, STP

  density = 3.700*mg/cm3 ;
  a = 83.80*g/mole ;
  G4Material* Kr  = new G4Material(name="Kr",z=36., a, density );



  // Carbone dioxide, CO2 STP

  density = 1.977*mg/cm3 ;
  G4Material* CarbonDioxide = new G4Material(name="CO2", density, nel=2) ;
  CarbonDioxide->AddElement(elC,1) ;
  CarbonDioxide->AddElement(elO,2) ;

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

  // 87.5% Xe + 7.5% CH4 + 5% C3H8, 20 C, 1 atm 

  density = 4.9196*mg/cm3 ;

  G4Material* XeCH4C3H8 = new G4Material(name="XeCH4C3H8"  , density, 
                                                             ncomponents=3);
  XeCH4C3H8->AddMaterial( Xe,       fractionmass = 0.971 ) ;
  XeCH4C3H8->AddMaterial( metane,   fractionmass = 0.010 ) ;
  XeCH4C3H8->AddMaterial( propane,  fractionmass = 0.019 ) ;

  // Propane in MWPC, 2 atm, 20 C

  //  density = 3.758*mg/cm3 ;
  density = 3.736*mg/cm3 ;
  G4Material* propaneDet = new G4Material(name="detC3H8",density,nel=2) ;
  propaneDet->AddElement(elC,3) ;
  propaneDet->AddElement(elH,8) ;

  // 80% Ar + 20% CO2, STP

  density = 1.8223*mg/cm3 ;      
  G4Material* Ar20CO2 = new G4Material(name="Ar20CO2"  , density, 
                                                             ncomponents=2);
  Ar20CO2->AddMaterial( Argon,           fractionmass = 0.783 ) ;
  Ar20CO2->AddMaterial( CarbonDioxide,   fractionmass = 0.217 ) ;

  // 93% Ar + 7% CH4, STP

  density = 1.709*mg/cm3 ;      
  G4Material* Ar7CH4 = new G4Material(name="Ar7CH4"  , density, 
                                                             ncomponents=2);
  Ar7CH4->AddMaterial( Argon,    fractionmass = 0.971 ) ;
  Ar7CH4->AddMaterial( metane,   fractionmass = 0.029 ) ;

  // 80% Xe + 20% CO2, STP

  density = 5.0818*mg/cm3 ;      
  G4Material* Xe20CO2 = new G4Material(name="Xe20CO2"  , density, 
                                                             ncomponents=2);
  Xe20CO2->AddMaterial( Xe,              fractionmass = 0.922 ) ;
  Xe20CO2->AddMaterial( CarbonDioxide,   fractionmass = 0.078 ) ;

  // 80% Kr + 20% CO2, STP

  density = 3.601*mg/cm3 ;      
  G4Material* Kr20CO2 = new G4Material(name="Kr20CO2"  , density, 
                                                             ncomponents=2);
  Kr20CO2->AddMaterial( Kr,              fractionmass = 0.89 ) ;
  Kr20CO2->AddMaterial( CarbonDioxide,   fractionmass = 0.11 ) ;

  // 80% He + 20% CO2, STP

  density = 0.5378*mg/cm3 ;      
  G4Material* He20CO2 = new G4Material(name="He20CO2"  , density, 
                                                             ncomponents=2);
  He20CO2->AddMaterial( He,              fractionmass = 0.265 ) ;
  He20CO2->AddMaterial( CarbonDioxide,   fractionmass = 0.735 ) ;




  //  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  //
  //  Create Sandia/PAI tables for given material 
  //

  G4int i, j, k, numOfMaterials, iSan, nbOfElements, sanIndex, row ;
  G4double maxEnergyTransfer, kineticEnergy ;
  G4double tau, gamma, bg2, beta2, rateMass, Tmax ;

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable() ;

  numOfMaterials = theMaterialTable->length();

  G4cout<<"Available materials under test : "<< G4endl<<G4endl ;
  outFile<<"Available materials under test : "<< G4endl<<G4endl ;

  for(k=0;k<numOfMaterials;k++)
  {
  G4cout <<k<<"\t"<< "  Material : " <<(*theMaterialTable)[k]->GetName() << G4endl ;
 outFile <<k<<"\t"<< "  Material : " <<(*theMaterialTable)[k]->GetName() << G4endl ;
  }
  G4String testName ;
  G4cout<<"Enter material name for test : "<<G4std::flush ;
  G4cin>>testName ;

  for(k=0;k<numOfMaterials;k++)
  {
    if((*theMaterialTable)[k]->GetName() != testName) continue ;

     outFile << "Material : " <<(*theMaterialTable)[k]->GetName() << G4endl ;
     G4cout << "Material : " <<(*theMaterialTable)[k]->GetName() << G4endl ;

     nbOfElements = (*theMaterialTable)[k]->GetNumberOfElements() ;

     G4cout<<"Sandia cof according old PAI stuff"<<G4endl<<G4endl ;
     outFile<<"Sandia cof according old PAI stuff"<<G4endl<<G4endl ;

     G4int* thisMaterialZ = new G4int[nbOfElements] ;
     for(iSan=0;iSan<nbOfElements;iSan++)
     {
        thisMaterialZ[iSan] = (G4int)(*theMaterialTable)[k]->
                                      GetElement(iSan)->GetZ() ;
     }
     G4SandiaTable sandia(k) ;
     sanIndex = sandia.SandiaIntervals(thisMaterialZ,nbOfElements) ;    
     sanIndex = sandia.SandiaMixing( thisMaterialZ ,
                             (*theMaterialTable)[k]->GetFractionVector() ,
				     nbOfElements,sanIndex) ;

     for(row=0;row<sanIndex-1;row++)
     {
       G4cout<<row+1<<"\t"<<sandia.GetPhotoAbsorpCof(row+1,0)/keV ;
       outFile<<row+1<<"  "<<sandia.GetPhotoAbsorpCof(row+1,0)/keV ;

       for(iSan=1;iSan<5;iSan++)
       {
         G4cout<<"\t"<<sandia.GetPhotoAbsorpCof(row+1,iSan) ;
	 // *(*theMaterialTable)[k]->GetDensity() ;

         outFile<<"  "<<sandia.GetPhotoAbsorpCof(row+1,iSan) ;
	 // *(*theMaterialTable)[k]->GetDensity() ;
       }
       G4cout<<G4endl ;
       outFile<<G4endl ;
     }
     G4cout<<G4endl ;
     outFile<<G4endl ;


     outFile<<G4endl ;
     maxEnergyTransfer = 100*keV ;
     gamma = 4.0 ;
     bg2 = gamma*gamma - 1 ;

     G4PAIxSection testPAI(k,maxEnergyTransfer,bg2) ;

     G4cout<<"Interval no."<<"\t"<<"Energy interval"<<G4endl<<G4endl ;
     outFile<<"Interval no."<<"\t"<<"Energy interval"<<G4endl<<G4endl ;

     for(j=1;j<=testPAI.GetIntervalNumber();j++)
     {
       G4cout<<j<<"\t\t"<<testPAI.GetEnergyInterval(j)/keV<<G4endl ;
       outFile<<j<<"\t\t"<<testPAI.GetEnergyInterval(j)/keV<<G4endl ;
     }
     G4cout<<G4endl ;
     outFile<<G4endl ;

     outFile<<"Actual spline size = "<<testPAI.GetSplineSize()<<G4endl ;
     outFile<<"Normalization Cof = "<<testPAI.GetNormalizationCof()<<G4endl ;
     outFile<<G4endl ;

     G4cout << "Actual spline size = "<<testPAI.GetSplineSize()<<G4endl ;
     G4cout <<"Normalization Cof = "<<testPAI.GetNormalizationCof()<<G4endl ;
     G4cout << G4endl ;

     outFile<<"Lorentz factor"<<"\t"<<"Max E transfer, kev"<<"\t"
          <<"<dE/dx>, keV/cm"<<"\t\t"<<"<dN/dx>, 1/cm"<<G4endl<<G4endl ;
   
     G4cout << "Lorentz factor"<<"\t"<<"Max E transfer, kev"<<"\t"
            << "<dE/dx>, keV/cm"<<"\t\t"<<"<dN/dx>, 1/cm"<<G4endl<<G4endl ;
   

     //   G4PAIxSection testPAIproton(k,maxEnergyTransfer) ;

     kineticEnergy = 110*MeV ;

     //     for(j=1;j<testPAIproton.GetNumberOfGammas();j++)

     for(j=1;j<30;j++)
     {
       tau      = kineticEnergy/proton_mass_c2 ;
       gamma    = tau +1.0 ;
       bg2      = tau*(tau + 2.0) ;
       beta2    = bg2/(gamma*gamma) ;
       rateMass = electron_mass_c2/proton_mass_c2 ;

       Tmax     = 2.0*electron_mass_c2*bg2
                   /(1.0+2.0*gamma*rateMass+rateMass*rateMass) ;

       if ( maxEnergyTransfer > Tmax)         
       {
          maxEnergyTransfer = Tmax ;
       }
       G4PAIxSection testPAIproton(k,maxEnergyTransfer,bg2) ;
      
       outFile << gamma << "\t"
               << maxEnergyTransfer/keV<<"\t\t"
               << testPAIproton.GetMeanEnergyLoss()*cm/keV << "\t\t"
               << testPAIproton.GetIntegralPAIxSection(1)*cm << "\t\t" << G4endl ;
       G4cout  << gamma << "\t"
               << maxEnergyTransfer/keV<<"\t\t"
               << testPAIproton.GetMeanEnergyLoss()*cm/keV << "\t\t"
               << testPAIproton.GetIntegralPAIxSection(1)*cm << "\t\t" << G4endl ;

       //   outFile<<testPAIproton.GetLorentzFactor(j)<<"\t"
       //          <<maxEnergyTransfer/keV<<"\t\t"
       //          <<testPAIproton.GetPAItable(0,j)*cm/keV<<"\t\t"
       //  	      <<testPAIproton.GetPAItable(1,j)*cm<<"\t\t"<<G4endl ;

       kineticEnergy *= 1.5 ;
     }
     G4cout<<G4endl ;
     outFile<<G4endl ;
  }

  G4String confirm ;
  G4cout<<"Enter 'y' , if you would like to get dE/dx-distribution : "<<G4std::flush ;
  G4cin>>confirm ;
  if(confirm != "y" ) return 1 ;
  G4cout<<G4endl ;

  for(k=0;k<numOfMaterials;k++)
  {
    G4cout <<k<< "  Material : " <<(*theMaterialTable)[k]->GetName() << G4endl ;
  } 
  G4cout<<"Enter material name for dE/dx-distribution : "<<G4std::flush ;
  G4cin>>testName ;
  G4cout<<G4endl ;

  G4int    iLoss, iStat, iStatMax, nGamma ;
  G4double energyLoss[50], Ebin, delta, delta1, delta2, delta3, step, y, pos ;
  G4double intProb[200], colDist, sum, fact, GF, lambda, aaa ;

  G4double alphaCrossTalk = -0.055, betaS = 0.2*0.4*keV ;
  G4int    spectrum[50] ;

  G4cout << " Enter nGamma 1<nGamma<10 : "  <<G4std::flush ;
  G4cin>>nGamma ;
  G4cout<<G4endl ;

  for(k=0;k<numOfMaterials;k++)
  {
     if((*theMaterialTable)[k]->GetName() != testName) continue ;

     G4cout << "Material : " <<(*theMaterialTable)[k]->GetName() << G4endl<<G4endl ;


     G4cout << " Enter Lorentz factor : "  <<G4std::flush ;
     G4cin>>gamma ;
     G4cout<<G4endl ;

     G4cout << " Enter step in mm : " <<G4std::flush ;
     G4cin>>step ;
     G4cout<<G4endl ;
     step *= mm ;

     G4cout << " Enter energy bin in keV : " <<G4std::flush ;
     G4cin>>Ebin ;
     G4cout<<G4endl ;
     Ebin *= keV ;

     G4cout << " Enter number of events : " <<G4std::flush ;
     G4cin>>iStatMax ;

     G4cout<<G4endl<<"Start dE/dx distribution"<<G4endl<<G4endl ;

     maxEnergyTransfer = 100*keV ;
     bg2               = gamma*gamma - 1 ;
     rateMass          = electron_mass_c2/proton_mass_c2 ;

     Tmax              = 2.0*electron_mass_c2*bg2
                          /(1.0+2.0*gamma*rateMass+rateMass*rateMass) ;

     if ( maxEnergyTransfer > Tmax)         maxEnergyTransfer = Tmax ;
       
     G4PAIxSection testPAIenergyLoss(k,maxEnergyTransfer,bg2) ;
 
     for( iLoss = 0 ; iLoss < 50 ; iLoss++ )
     {
        energyLoss[iLoss] = Ebin*iLoss ;
        spectrum[iLoss] = 0 ;
     }
     for(iStat=0;iStat<iStatMax;iStat++)
     {

       //   aaa = (G4double)nGamma ;
       //   lambda = aaa/step ;
       //   colDist = RandGamma::shoot(aaa,lambda) ;

       //  delta = testPAIenergyLoss.GetStepEnergyLoss(colDist) ;

       //  delta = testPAIenergyLoss.GetStepEnergyLoss(step) ;

          delta1 = testPAIenergyLoss.GetStepEnergyLoss(step) ;

          delta = RandGauss::shoot(delta1,0.3*delta1) ;
          if( delta < 0.0 ) delta = 0.0 ;

       //   delta2 = testPAIenergyLoss.GetStepEnergyLoss(step) ;
       //   delta3 = testPAIenergyLoss.GetStepEnergyLoss(step) ;
 
       //   delta = alphaCrossTalk*delta1 + 
       //         delta2 + alphaCrossTalk*delta3 - betaS ;

       for(iLoss=0;iLoss<50;iLoss++)
       {
         if(delta <= energyLoss[iLoss]) break ;
       }
       spectrum[iLoss-1]++ ;
     }
     G4double meanLoss = 0.0 ;

     outFile<<"E, keV"<<"\t\t"<<"Distribution"<<G4endl<<G4endl ;
     G4cout<<"E, keV"<<"\t\t"<<"Distribution"<<G4endl<<G4endl ;
     G4cout<<G4endl ;
     for(iLoss=0;iLoss<50;iLoss++) // with last bin
     {
       fileOut<<energyLoss[iLoss]/keV<<"\t\t"<<spectrum[iLoss]<<G4endl ;
       G4cout<<energyLoss[iLoss]/keV<<"\t\t"<<spectrum[iLoss]<<G4endl ;
       meanLoss +=energyLoss[iLoss]*spectrum[iLoss] ;
     }
     G4cout<<G4endl ;
     G4cout<<"Mean loss over spectrum = "<<meanLoss/keV/iStatMax<<" keV"<<G4endl ;
  }

  G4int exit = 1 ;

  while(exit)
  {
     G4cout<<"Enter 'y' , if you would like to compare with exp. data : "<<G4std::flush ;
     G4cin>>confirm ;
     if(confirm != "y" ) break ;
     G4cout<<G4endl ;

     // Read experimental data file

     G4double delExp[200], distr[200], deltaBin, sumPAI, sumExp ;
     G4int numberOfExpPoints ;

     G4cout<<G4endl ;
     G4cout << " Enter number of experimental points : " <<G4std::flush ;
     G4cin>>numberOfExpPoints ;
     G4cout<<G4endl ;
     G4cout << " Enter energy bin in keV : " <<G4std::flush ;
     G4cin>>deltaBin ;
     G4cout<<G4endl ;
     deltaBin *= keV ;

     G4std::ifstream fileRead ;
     fileRead.open("input.dat") ;
     for(i=0;i<numberOfExpPoints;i++)
     {
       fileRead>>delExp[i]>>distr[i] ;
       delExp[i] *= keV ;
       G4cout<<i<<"\t"<<delExp[i]<<"\t"<<distr[i]<<G4endl ;
     }
     fileRead.close() ;

     // Adjust statistics of experiment to PAI simulation

     sumExp = 0.0 ;
     for(i=0;i<numberOfExpPoints;i++) sumExp +=distr[i] ;
     sumExp *= deltaBin ;

     sumPAI = 0.0 ;
     for(i=0;i<49;i++) sumPAI +=spectrum[i] ;
     sumPAI *= Ebin ;

     for(i=0;i<numberOfExpPoints;i++) distr[i] *= sumPAI/sumExp ;

     for(i=0;i<numberOfExpPoints;i++)
     {
       fileWrite<<delExp[i]/keV<<"\t"<<distr[i]<<G4endl ;
       G4cout<<delExp[i]/keV<<"\t"<<distr[i]<<G4endl ;
     }
     exit = 0 ;
  }

  G4cout<<"Enter 'y' , if you would like to get most probable delta : "<<G4std::flush ;
  G4cin>>confirm ;
  if(confirm != "y" ) return 1 ;
  G4cout<<G4endl ;

  G4int kGamma, iMPLoss, maxSpectrum, iMax ;
  G4double mpDelta[50], meanDelta[50], rrMP[50], rrMean[50] ; 
  G4double mpLoss, tmRatio, mpSum, mpStat ;

  G4double aGamma[33] = 
  {
    4.0, 1.5, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0, // 13
    20., 40.0, 60.0, 80.0, 100.0, 200.0, 400.0, 600.0, 800.0, 1000.0, // 23
    2000.0, 4000.0, 6000.0, 8000.0, 100000.0, 20000.0,                // 29
    40000.0, 60000.0, 80000.0, 100000.0                               // 33
  } ;

  for(k=0;k<numOfMaterials;k++)
  {
    G4cout <<k<< "  Material : " <<(*theMaterialTable)[k]->GetName() << G4endl ;
  } 
  G4cout<<"Enter material name for dE/dx-distribution : "<<G4std::flush ;
  G4cin>>testName ;
  G4cout<<G4endl ;


  for(k=0;k<numOfMaterials;k++)
  {
     if((*theMaterialTable)[k]->GetName() != testName) continue ;

     G4cout << "Material : " <<(*theMaterialTable)[k]->GetName() << G4endl<<G4endl ;

     G4cout << " Enter nGamma 1<nGamma<10 : "  <<G4std::flush ;
     G4cin>>nGamma ;
     G4cout<<G4endl ;


     G4cout << " Enter step in mm : " <<G4std::flush ;
     G4cin>>step ;
     G4cout<<G4endl ;
     step *= mm ;

     G4cout << " Enter energy bin in keV : " <<G4std::flush ;
     G4cin>>Ebin ;
     G4cout<<G4endl ;
     Ebin *= keV ;

     G4cout << " Enter trancated mean ration <1.0 : "  <<G4std::flush ;
     G4cin>>tmRatio ;
     G4cout<<G4endl ;


     G4cout << " Enter number of events : " <<G4std::flush ;
     G4cin>>iStatMax ;
     G4cout<<G4endl ;

     G4cout<<"no."<<"\t"<<"Gamma"<<"\t"<<"Rel. rise"<<"\t"<<"M.P. loss, keV"
           <<"\t"<<"Mean loss, keV"<<G4endl<<G4endl ;
     //   outFile<<"no."<<"\t"<<"Gamma"<<"\t"<<"M.P. loss, keV"
     //      <<"\t"<<"Mean loss, keV"<<G4endl<<G4endl ;
     

     // gamma = 1.1852 ;

     for(kGamma=0;kGamma<33;kGamma++)
     {
       //    G4cout<<G4endl<<"Start dE/dx distribution"<<G4endl<<G4endl ;

       gamma = aGamma[kGamma] ;
       maxEnergyTransfer = 100*keV ;
       bg2               = gamma*gamma - 1 ;
       rateMass          = electron_mass_c2/proton_mass_c2 ;

       Tmax              = 2.0*electron_mass_c2*bg2
                          /(1.0+2.0*gamma*rateMass+rateMass*rateMass) ;

       if ( maxEnergyTransfer > Tmax)         maxEnergyTransfer = Tmax ;
       
       G4PAIxSection testPAIenergyLoss(k,maxEnergyTransfer,bg2) ;
 
       for( iLoss = 0 ; iLoss < 50 ; iLoss++ )
       {
         energyLoss[iLoss] = Ebin*iLoss ;
         spectrum[iLoss] = 0 ;
       }
       for(iStat=0;iStat<iStatMax;iStat++)
       {

         //   aaa = (G4double)nGamma ;
         //   lambda = aaa/step ;
         //   colDist = RandGamma::shoot(aaa,lambda) ;

         //  delta = testPAIenergyLoss.GetStepEnergyLoss(colDist) ;

         delta = testPAIenergyLoss.GetStepEnergyLoss(step) ;

         //   delta1 = testPAIenergyLoss.GetStepEnergyLoss(step) ;
         //   delta2 = testPAIenergyLoss.GetStepEnergyLoss(step) ;
         //   delta3 = testPAIenergyLoss.GetStepEnergyLoss(step) ;
 
         //   delta = alphaCrossTalk*delta1 + 
         //         delta2 + alphaCrossTalk*delta3 - betaS ;

         for(iLoss=0;iLoss<50;iLoss++)
         {
           if(delta <= energyLoss[iLoss]) break ;
         }
         spectrum[iLoss-1]++ ;
       }
       G4int sumStat = 0 ;
       for(iLoss=0;iLoss<49;iLoss++) // without last bin
       {
         sumStat += spectrum[iLoss] ;
         if( sumStat > tmRatio*iStatMax  ) break ;
       }
       if(iLoss == 50) iLoss-- ;
       iMPLoss = iLoss ;
       G4double meanLoss = 0.0 ;
       maxSpectrum = 0 ;

       for(iLoss=0;iLoss<iMPLoss;iLoss++) // without last bin
       {
	 // fileOut<<energyLoss[iLoss]/keV<<"\t\t"<<spectrum[iLoss]<<G4endl ;
	 //  G4cout<<energyLoss[iLoss]/keV<<"\t\t"<<spectrum[iLoss]<<G4endl ;

         meanLoss += energyLoss[iLoss]*spectrum[iLoss] ;

         if( spectrum[iLoss] > maxSpectrum )
	 {
           maxSpectrum = spectrum[iLoss]   ;
           mpLoss      = energyLoss[iLoss] ;
           iMax = iLoss ;
	 }
       }
       mpSum  = 0. ;
       mpStat = 0 ;
       for(iLoss = iMax-5;iLoss<=iMax+5;iLoss++)
       {
         mpSum += energyLoss[iLoss]*spectrum[iLoss] ;
         mpStat += spectrum[iLoss] ;
       }
       mpLoss = mpSum/mpStat ;
       mpLoss /= keV ;
       meanLoss /= keV*sumStat ;
       meanDelta[kGamma] = meanLoss ;
       mpDelta[kGamma] = mpLoss ;

       if(kGamma > 0)
       {
         rrMP[kGamma] = mpLoss/mpDelta[0] ;
         G4cout<<kGamma<<"\t"<<gamma<<"\t"<<rrMP[kGamma]<<"\t"<<mpLoss<<G4endl ;
	 //  outFile<<gamma<<"\t"<<rrMP[kGamma]<<G4endl ;
         fileWrite1<<gamma<<"\t"<<rrMP[kGamma]<<G4endl ;
       }

       //  gamma *= 1.5 ;
    }
    G4cout<<G4endl ;
    outFile<<G4endl ;
  }   

   return EXIT_SUCCESS;

}










