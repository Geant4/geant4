// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PAIxSectionTest.cc,v 1.2 1999-10-27 09:24:43 grichine Exp $
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
#include <fstream.h>
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4SandiaTable.hh"

#include "G4PAIonisation.hh"
#include "G4PAIxSection.hh"

int main()
{
   ofstream outFile("PAIdEdx.out", ios::out ) ;
   outFile.setf( ios::scientific, ios::floatfield );

   ofstream fileOut("PAIdistribution.out", ios::out ) ;
   //  fileOut.setf( ios::scientific, ios::floatfield );

// -------------------------- Create materials -------------------------------  
   

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

  density = 3.758*mg/cm3 ;
  G4Material* propaneDet = new G4Material(name="detC3H8",density,nel=2) ;
  propaneDet->AddElement(elC,3) ;
  propaneDet->AddElement(elH,8) ;

  // 80% Ar + 20% CO2, STP

  density = 1.8223*mg/cm3 ;      
  G4Material* Ar_80CO2_20 = new G4Material(name="Ar20CO2"  , density, 
                                                             ncomponents=2);
  Ar_80CO2_20->AddMaterial( Argon,           fractionmass = 0.783 ) ;
  Ar_80CO2_20->AddMaterial( CarbonDioxide,   fractionmass = 0.217 ) ;

  // 93% Ar + 7% CH4, STP

  density = 1.709*mg/cm3 ;      
  G4Material* Ar7CH4 = new G4Material(name="Ar7CH4"  , density, 
                                                             ncomponents=2);
  Ar7CH4->AddMaterial( Argon,    fractionmass = 0.971 ) ;
  Ar7CH4->AddMaterial( metane,   fractionmass = 0.029 ) ;

  // 80% Xe + 20% CO2, STP

  density = 5.0818*mg/cm3 ;      
  G4Material* Xe_80CO2_20 = new G4Material(name="Xe20CO2"  , density, 
                                                             ncomponents=2);
  Xe_80CO2_20->AddMaterial( Xe,              fractionmass = 0.922 ) ;
  Xe_80CO2_20->AddMaterial( CarbonDioxide,   fractionmass = 0.078 ) ;




  //  G4cout << *(G4Material::GetMaterialTable()) << endl;

  //
  //  Create Sandia/PAI tables for given material 
  //

  G4int i, j, k, numOfMaterials, iSan, nbOfElements, sanIndex, row ;
  G4double maxEnergyTransfer, kineticEnergy ;
  G4double tau, gamma, bg2, beta2, rateMass, Tmax ;

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable() ;

  numOfMaterials = theMaterialTable->length();

  G4cout<<"Available materials under test : "<< endl<<endl ;
  outFile<<"Available materials under test : "<< endl<<endl ;

  for(k=0;k<numOfMaterials;k++)
  {
  G4cout <<k<<"\t"<< "  Material : " <<(*theMaterialTable)[k]->GetName() << endl ;
 outFile <<k<<"\t"<< "  Material : " <<(*theMaterialTable)[k]->GetName() << endl ;
  }
  G4String testName ;
  G4cout<<"Enter material name for test : "<<flush ;
  cin>>testName ;

  for(k=0;k<numOfMaterials;k++)
  {
     if((*theMaterialTable)[k]->GetName() != testName) continue ;

     outFile << "Material : " <<(*theMaterialTable)[k]->GetName() << endl ;
     G4cout << "Material : " <<(*theMaterialTable)[k]->GetName() << endl ;

     nbOfElements = (*theMaterialTable)[k]->GetNumberOfElements() ;

     G4cout<<"Sandia cof according old PAI stuff"<<endl<<endl ;
     outFile<<"Sandia cof according old PAI stuff"<<endl<<endl ;

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
       G4cout<<endl ;
       outFile<<endl ;
     }
     G4cout<<endl ;
     outFile<<endl ;


     outFile<<endl ;
     maxEnergyTransfer = 100*keV ;
     gamma = 4.0 ;
     bg2 = gamma*gamma - 1 ;

     G4PAIxSection testPAI(k,maxEnergyTransfer,bg2) ;

     G4cout<<"Interval no."<<"\t"<<"Energy interval"<<endl<<endl ;
     outFile<<"Interval no."<<"\t"<<"Energy interval"<<endl<<endl ;

     for(j=1;j<=testPAI.GetIntervalNumber();j++)
     {
       G4cout<<j<<"\t\t"<<testPAI.GetEnergyInterval(j)/keV<<endl ;
       outFile<<j<<"\t\t"<<testPAI.GetEnergyInterval(j)/keV<<endl ;
     }
     G4cout<<endl ;
     outFile<<endl ;

     outFile<<"Actual spline size = "<<testPAI.GetSplineSize()<<endl ;
     outFile<<"Normalization Cof = "<<testPAI.GetNormalizationCof()<<endl ;
     outFile<<endl ;

     G4cout << "Actual spline size = "<<testPAI.GetSplineSize()<<endl ;
     G4cout <<"Normalization Cof = "<<testPAI.GetNormalizationCof()<<endl ;
     G4cout << endl ;

     outFile<<"Lorentz factor"<<"\t"<<"Max E transfer, kev"<<"\t"
          <<"<dE/dx>, keV/cm"<<"\t\t"<<"<dN/dx>, 1/cm"<<endl<<endl ;
   
     G4cout << "Lorentz factor"<<"\t"<<"Max E transfer, kev"<<"\t"
            << "<dE/dx>, keV/cm"<<"\t\t"<<"<dN/dx>, 1/cm"<<endl<<endl ;
   

     //   G4PAIxSection testPAIproton(k,maxEnergyTransfer) ;

     kineticEnergy = 100*MeV ;

     //     for(j=1;j<testPAIproton.GetNumberOfGammas();j++)

     for(j=1;j<20;j++)
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
               << testPAIproton.GetIntegralPAIxSection(1)*cm << "\t\t" << endl ;
       G4cout  << gamma << "\t"
               << maxEnergyTransfer/keV<<"\t\t"
               << testPAIproton.GetMeanEnergyLoss()*cm/keV << "\t\t"
               << testPAIproton.GetIntegralPAIxSection(1)*cm << "\t\t" << endl ;

       //   outFile<<testPAIproton.GetLorentzFactor(j)<<"\t"
       //          <<maxEnergyTransfer/keV<<"\t\t"
       //          <<testPAIproton.GetPAItable(0,j)*cm/keV<<"\t\t"
       //  	      <<testPAIproton.GetPAItable(1,j)*cm<<"\t\t"<<endl ;

       kineticEnergy *= 2.0 ;
     }
     G4cout<<endl ;
     outFile<<endl ;
  }

  G4String confirm ;
  G4cout<<"Enter 'y' , if you would like to get dE/dx-distribution : "<<flush ;
  cin>>confirm ;
  if(confirm != "y" ) return 1 ;
  G4cout<<endl ;
  for(k=0;k<numOfMaterials;k++)
  {
    G4cout <<k<< "  Material : " <<(*theMaterialTable)[k]->GetName() << endl ;
  }
 
  G4cout<<"Enter material name for dE/dx-distribution : "<<flush ;
  cin>>testName ;
  G4cout<<endl ;
  for(k=0;k<numOfMaterials;k++)
  {
     if((*theMaterialTable)[k]->GetName() != testName) continue ;

     G4cout << "Material : " <<(*theMaterialTable)[k]->GetName() << endl<<endl ;
     G4cout<<"Start dE/dx distribution"<<endl<<endl ;

     G4int    iLoss, iStat, iStatMax ;
     G4double energyLoss[50], Ebin, delta, step ;
     G4int    spectrum[50] ;

     G4cout << " Enter Lorentz factor : "  <<flush ;
     cin>>gamma ;
     G4cout<<endl ;

     G4cout << " Enter step in mm : " <<flush ;
     cin>>step ;
     G4cout<<endl ;
     step *= mm ;

     G4cout << " Enter energy bin in keV : " <<flush ;
     cin>>Ebin ;
     G4cout<<endl ;
     Ebin *= keV ;

     G4cout << " Enter number of events : " <<flush ;
     cin>>iStatMax ;

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
       delta = testPAIenergyLoss.GetStepEnergyLoss(step) ;

       for(iLoss=0;iLoss<50;iLoss++)
       {
         if(delta <= energyLoss[iLoss]) break ;
       }
       spectrum[iLoss-1]++ ;
     }
     G4double meanLoss = 0.0 ;

     outFile<<"E, keV"<<"\t\t"<<"Distribution"<<endl<<endl ;
     G4cout<<"E, keV"<<"\t\t"<<"Distribution"<<endl<<endl ;
     G4cout<<endl ;
     for(iLoss=0;iLoss<50;iLoss++) // with last bin
     {
       fileOut<<energyLoss[iLoss]/keV<<"\t\t"<<spectrum[iLoss]<<endl ;
       G4cout<<energyLoss[iLoss]/keV<<"\t\t"<<spectrum[iLoss]<<endl ;
       meanLoss +=energyLoss[iLoss]*spectrum[iLoss] ;
     }
     G4cout<<endl ;
     G4cout<<"Mean loss over spectrum = "<<meanLoss/keV/iStatMax<<" keV"<<endl ;
   }

   return EXIT_SUCCESS;

}







