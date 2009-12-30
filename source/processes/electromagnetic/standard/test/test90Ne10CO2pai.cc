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
// $Id: test90Ne10CO2pai.cc,v 1.8 2009-12-30 12:57:41 grichine Exp $
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
// 11.04.08, V. Grichine test of PAI predictions for 90% Ne + 10% CO2
//                       ALICE TPC gas mixture 
// 8.12.09   V. Grichine update for T2K gas mixture

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

// Peter:
// From the AliRoot code in $ALICE_ROOT/TPC/AliTPCv2.cxx:
//  const Float_t kprim = 14.35; // number of primary collisions per 1 cm
//  const Float_t kpoti = 20.77e-9; // first ionization potential for Ne/CO2
//  const Float_t kwIon = 35.97e-9; // energy for the ion-electron pair creation
//
// kprim is the number of primary collisions per cm for a MIP!
// kpoti = I.
// kwIon = W.

G4double FitALICE(G4double bg)
{
  //
  // Bethe-Bloch energy loss formula from ALICE TPC TRD
  //
  const G4double kp1 = 0.76176e-1;
  const G4double kp2 = 10.632;
  const G4double kp3 = 0.13279e-4;
  const G4double kp4 = 1.8631;
  const G4double kp5 = 1.9479;
  const G4double dn1 = 14.35;

  G4double dbg  = (G4double) bg;

  G4double beta = dbg/std::sqrt(1.+dbg*dbg);

  G4double aa   = std::pow(beta,kp4);
  G4double bb   = std::pow(1./dbg,kp5);


  bb  = std::log(kp3 + bb);

  G4double result = ( kp2 - aa - bb)*kp1/aa;

  result *= dn1;

  return result;
}


G4double FitBichsel(G4double bg)
{
  //
  // Primary ionisation from Hans Bichsel fit
  //
  const G4double kp1 = 0.686e-1;
  const G4double kp2 = 11.714;
  const G4double kp3 = 0.218e-4;
  const G4double kp4 = 1.997;
  const G4double kp5 = 2.133;
  const G4double dn1 = 13.32;

  G4double dbg  = (G4double) bg;

  G4double beta = dbg/std::sqrt(1.+dbg*dbg);

  G4double aa   = std::pow(beta,kp4);
  G4double bb   = std::pow(1./dbg,kp5);


  bb=std::log(kp3 + bb);

  G4double result = ( kp2 - aa - bb)*kp1/aa;

  result *= dn1;

  return result;
  
}

////////////////////////////////////////////////////////////


G4double GetIonisation(G4double transfer)
{
  G4double W  = 34.75*eV;
  G4double I1 = 13.62*eV; // first ionisation potential in mixture
  I1 *= 0.9;

  G4double       result = W;

  // result /= 1.-I1/transfer;

  return transfer/result;

}


/////////////////////////////////////////////////



int main()
{
  // std::ofstream outFile("90Ne10CO2pai.dat", std::ios::out );
   std::ofstream outFile("e5GeVt2kPhilippe.dat", std::ios::out );
   outFile.setf( std::ios::scientific, std::ios::floatfield );

   std::ofstream fileOut("PAICerPlasm90Ne10CO2.dat", std::ios::out );
   fileOut.setf( std::ios::scientific, std::ios::floatfield );

   //  std::ifstream fileRead("exp.dat", std::ios::out );
   //  fileRead.setf( std::ios::scientific, std::ios::floatfield );

   std::ofstream fileWrite("exp.dat", std::ios::out );
   fileWrite.setf( std::ios::scientific, std::ios::floatfield );

   std::ofstream fileWrite1("mprrpai.dat", std::ios::out );
   fileWrite1.setf( std::ios::scientific, std::ios::floatfield );

// Create materials  
   

  G4int iz , n,  nel, ncomponents;
  G4double a, z, ez, density , temperature, pressure, fractionmass;
  G4State state;
  G4String name, symbol;


  a = 1.01*g/mole;
  G4Isotope* ih1 = new G4Isotope("Hydrogen",iz=1,n=1,a);

  a = 2.01*g/mole;
  G4Isotope* ih2 = new G4Isotope("Deuterium",iz=1,n=2,a);

  G4Element* elH = new G4Element(name="Hydrogen",symbol="H",2);
  elH->AddIsotope(ih1,.999);
  elH->AddIsotope(ih2,.001);

  a = 12.01*g/mole;
  G4Element* elC = new G4Element(name="Carbon",symbol="C", ez=6., a);

  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", ez=7., a);

  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxygen",symbol="O", ez=8., a);

  
  a = 19.00*g/mole;
  G4Element* elF  = new G4Element(name="Fluorine", symbol="F", z=9., a);
 
  a = 39.948*g/mole;
  G4Element* elAr = new G4Element(name="Argon", symbol="Ar", z=18., a);

  // Neon as detector gas, STP

  density = 0.900*mg/cm3;
  a = 20.179*g/mole;
  G4Material* Ne  = new G4Material(name="Ne",z=10., a, density );

  // Carbone dioxide, CO2 STP

  density = 1.977*mg/cm3;
  G4Material* CarbonDioxide = new G4Material(name="CO2", density, nel=2);
  CarbonDioxide->AddElement(elC,1);
  CarbonDioxide->AddElement(elO,2);

  

  // 90% Ne + 10% CO2, STP

  density = 1.0077*mg/cm3;      
  G4Material* Ne10CO2 = new G4Material(name="Ne10CO2"  , density, 
                          ncomponents=2);
  Ne10CO2->AddMaterial( Ne,           fractionmass = 0.8038 );
  Ne10CO2->AddMaterial( CarbonDioxide,   fractionmass = 0.1962 );

  density *= 273./293.;

  G4Material* Ne10CO2T293 = new G4Material(name="Ne10CO2T293"  , density, 
                           ncomponents=2);
  Ne10CO2T293->AddMaterial( Ne,              fractionmass = 0.8038 );
  Ne10CO2T293->AddMaterial( CarbonDioxide,   fractionmass = 0.1962 );


  density = 1.25053*mg/cm3;       // STP
  G4Material* Nitrogen = new G4Material(name="N2"  , density, ncomponents=1);
  Nitrogen->AddElement(elN, 2);

  // 85.7% Ne + 9.5% CO2 +4.8% N2, STP

  density = 1.0191*mg/cm3;      
  G4Material* Ne857CO295N2 = new G4Material(name="Ne857CO295N2"  , density, 
                             ncomponents=3);
  Ne857CO295N2->AddMaterial( Ne,            fractionmass = 0.7568 );
  Ne857CO295N2->AddMaterial( CarbonDioxide, fractionmass = 0.1843 );
  Ne857CO295N2->AddMaterial( Nitrogen,      fractionmass = 0.0589 );

  density *= 273./292.;
  density *= 0.966/1.01325;

  // G4cout<<"density of Ne857CO295N2T292 = "<<density*cm3/mg<<"  mg/cm3"<<G4endl;

  G4Material* Ne857CO295N2T292 = new G4Material(name="Ne857CO295N2T292"  , density, 
                             ncomponents=3);
  Ne857CO295N2T292->AddMaterial( Ne,            fractionmass = 0.76065 );
  Ne857CO295N2T292->AddMaterial( CarbonDioxide, fractionmass = 0.18140 );
  Ne857CO295N2T292->AddMaterial( Nitrogen,      fractionmass = 0.05795 );


  
  // Ar as detector gas,STP

  density = 1.7836*mg/cm3 ;       // STP
  G4Material* Argon = new G4Material(name="Argon"  , density, ncomponents=1);
  Argon->AddElement(elAr, 1);
  /*
  // iso-Butane (methylpropane), STP

  density = 2.67*mg/cm3 ;
  G4Material* isobutane = new G4Material(name="isoC4H10",density,nel=2) ;
  isobutane->AddElement(elC,4) ;
  isobutane->AddElement(elH,10) ;

  // CF4 from ATLAS TRT estimation

  G4double TRT_CF4_density = 3.9*mg/cm3;
  G4Material* TRT_CF4 = new G4Material(name="TRT_CF4", TRT_CF4_density, nel=2,
                                           kStateGas,293.15*kelvin,1.*atmosphere);
  */

  // Philippe Gros T2K mixture version
  // Argon                                                                                                   
  /*                                          
  density     = 1.66*mg/cm3;
  pressure    = 1*atmosphere;
  temperature = 288.15*kelvin;
  G4Material* Argon = new G4Material(name="Ar", // z=18., a=39.948*g/mole,
				     density, ncomponents=1); // kStateGas,temperature,pressure);
  Argon->AddElement(elAr, 1);
  */
   
  // IsoButane    
                                                                                                    
  density     = 2.51*mg/cm3;
  G4Material* Isobu = new G4Material(name="isoC4H10", z=34.,a=58.123*g/mole, 
                                       density, kStateGas,temperature,pressure);

  // Tetrafluoromethane   
                                                                                            
  density = 3.72*mg/cm3;
  G4Material* FlMet = new
  G4Material(name="CF4",z=42.,a=88.01*g/mole,density,kStateGas,temperature,pressure);

  // Argon + 3% tetrafluoromethane  + 2% iso-butane     
                                                            
  density = 1.748*mg/cm3;
  G4Material* t2kGasMixture = new G4Material(name="t2kGasMixture", density, ncomponents=3);
  t2kGasMixture->AddMaterial(Argon, fractionmass = 90.9*perCent);
  t2kGasMixture->AddMaterial(FlMet, fractionmass = 6.3*perCent);
  t2kGasMixture->AddMaterial(Isobu, fractionmass = 2.8*perCent);


  G4int i, j, jMax, k, numOfMaterials, iSan, nbOfElements, sanIndex, row;
  G4double maxEnergyTransfer, kineticEnergy;
  G4double tau, gamma, bg2, bg, beta2, rateMass, Tmax, Tmin, Tkin;

  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();

  numOfMaterials = theMaterialTable->size();

  G4cout<<"Available materials under test : "<< G4endl<<G4endl;
  // outFile<<"Available materials under test : "<< G4endl<<G4endl;

  for( k = 0; k < numOfMaterials; k++ )
  {
    G4cout <<k<<"\t"<< "  Material : " <<(*theMaterialTable)[k]->GetName() << G4endl;
    // outFile <<k<<"\t"<< "  Material : " <<(*theMaterialTable)[k]->GetName() << G4endl;
  }



  // G4String testName = "N2";
  // G4String testName = "Ne10CO2";
  G4String testName = "t2kGasMixture";
  // G4String testName = "Ar";
  //  G4String testName = "Argon";
  // G4String testName = "Ne10CO2T293";
  // G4String testName = "Ne857CO295N2T292";


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

     for( iSan = 0; iSan < nbOfElements; iSan++ )
     {
        thisMaterialZ[iSan] = (G4int)(*theMaterialTable)[k]->
                                      GetElement(iSan)->GetZ();
     }
     G4SandiaTable sandia(k);
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
     maxEnergyTransfer = 100*keV;
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

     // kineticEnergy = 5*GeV; // for electrons

     kineticEnergy = 5*GeV*proton_mass_c2/electron_mass_c2;

     // kineticEnergy = 5*GeV;

     //     for(j=1;j<testPAIproton.GetNumberOfGammas();j++)

     // jMax = 70; // 70;
     jMax = 1; // 70;

     // outFile<<jMax<<G4endl;

     for( j = 0; j < jMax; j++ )
     {
       tau      = kineticEnergy/proton_mass_c2;
       // tau      = kineticEnergy/electron_mass_c2;
       gamma    = tau +1.0;
       bg2      = tau*(tau + 2.0);
       bg = std::sqrt(bg2);
       beta2    = bg2/(gamma*gamma);
       G4cout<<"bg = "<<bg<<";  b2 = "<<beta2<<G4endl<<G4endl;
       rateMass = electron_mass_c2/proton_mass_c2;
       
       Tmax = 2.0*electron_mass_c2*bg2/(1.0+2.0*gamma*rateMass+rateMass*rateMass);
       // Tmax = 0.5*kineticEnergy;

       Tkin = maxEnergyTransfer;

       if ( maxEnergyTransfer > Tmax)         
       {
          Tkin = Tmax;
       }
       if ( Tmax <= Tmin + 0.5*eV )         
       {
          Tkin = Tmin + 0.5*eV;
       }
       G4PAIxSection testPAIproton(k,Tkin,bg2);
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
       
       outFile<<testPAIproton.GetSplineSize()-1<<G4endl;

       for( i = 1; i < testPAIproton.GetSplineSize(); i++)
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
       


       kineticEnergy *= 1.41;        // was 1.4; 1.5;
     }

     G4cout<<G4endl;
     // outFile<<G4endl;
  }


  return 1;  // end of test















  G4String confirm;
  G4cout<<"Enter 'y' , if you would like to get dE/dx-distribution : "
        <<std::flush;

  G4cin>>confirm;
  if(confirm != "y" ) return 1;
  G4cout<<G4endl;

  for(k=0;k<numOfMaterials;k++)
  {
    G4cout <<k<< "  Material : " <<(*theMaterialTable)[k]->GetName() << G4endl;
  } 
  G4cout<<"Enter material name for dE/dx-distribution : "<<std::flush;
  G4cin>>testName;
  G4cout<<G4endl;

  G4int    iLoss, iStat, iStatMax, nGamma;
  G4double energyLoss[50], Ebin, delta, delta1, delta2, delta3, step, y, pos;
  G4double intProb[200], colDist, sum, fact, GF, lambda, aaa;

  G4double alphaCrossTalk = -0.055, betaS = 0.2*0.4*keV;
  G4int    spectrum[50];

  G4cout << " Enter nGamma 1<nGamma<10 : "  <<std::flush;
  G4cin>>nGamma;
  G4cout<<G4endl;

  for(k=0;k<numOfMaterials;k++)
  {
     if((*theMaterialTable)[k]->GetName() != testName) continue;

     G4cout << "Material : " <<(*theMaterialTable)[k]->GetName() << G4endl<<G4endl;


     G4cout << " Enter Lorentz factor : "  <<std::flush;
     G4cin>>gamma;
     G4cout<<G4endl;

     G4cout << " Enter step in mm : " <<std::flush;
     G4cin>>step;
     G4cout<<G4endl;
     step *= mm;

     G4cout << " Enter energy bin in keV : " <<std::flush;
     G4cin>>Ebin;
     G4cout<<G4endl;
     Ebin *= keV;

     G4cout << " Enter number of events : " <<std::flush;
     G4cin>>iStatMax;

     G4cout<<G4endl<<"Start dE/dx distribution"<<G4endl<<G4endl;

     maxEnergyTransfer = 100*keV;
     bg2               = gamma*gamma - 1;
     rateMass          = electron_mass_c2/proton_mass_c2;

     Tmax              = 2.0*electron_mass_c2*bg2
                          /(1.0+2.0*gamma*rateMass+rateMass*rateMass);

     if ( maxEnergyTransfer > Tmax)         maxEnergyTransfer = Tmax;
       
     G4PAIxSection testPAIenergyLoss(k,maxEnergyTransfer,bg2);
 
     for( iLoss = 0; iLoss < 50; iLoss++ )
     {
        energyLoss[iLoss] = Ebin*iLoss;
        spectrum[iLoss] = 0;
     }
     for(iStat=0;iStat<iStatMax;iStat++)
     {

       //   aaa = (G4double)nGamma;
       //   lambda = aaa/step;
       //   colDist = RandGamma::shoot(aaa,lambda);

       //  delta = testPAIenergyLoss.GetStepEnergyLoss(colDist);

       //  delta = testPAIenergyLoss.GetStepEnergyLoss(step);

          delta1 = testPAIenergyLoss.GetStepEnergyLoss(step);

          delta = G4RandGauss::shoot(delta1,0.3*delta1);
          if( delta < 0.0 ) delta = 0.0;

       //   delta2 = testPAIenergyLoss.GetStepEnergyLoss(step);
       //   delta3 = testPAIenergyLoss.GetStepEnergyLoss(step);
 
       //   delta = alphaCrossTalk*delta1 + 
       //         delta2 + alphaCrossTalk*delta3 - betaS;

       for(iLoss=0;iLoss<50;iLoss++)
       {
         if(delta <= energyLoss[iLoss]) break;
       }
       spectrum[iLoss-1]++;
     }
     G4double meanLoss = 0.0;

     outFile<<"E, keV"<<"\t\t"<<"Distribution"<<G4endl<<G4endl;
     G4cout<<"E, keV"<<"\t\t"<<"Distribution"<<G4endl<<G4endl;
     G4cout<<G4endl;
     for(iLoss=0;iLoss<50;iLoss++) // with last bin
     {
       fileOut<<energyLoss[iLoss]/keV<<"\t\t"<<spectrum[iLoss]<<G4endl;
       G4cout<<energyLoss[iLoss]/keV<<"\t\t"<<spectrum[iLoss]<<G4endl;
       meanLoss +=energyLoss[iLoss]*spectrum[iLoss];
     }
     G4cout<<G4endl;
     G4cout<<"Mean loss over spectrum = "<<meanLoss/keV/iStatMax<<" keV"<<G4endl;
  }

  G4int exit = 1;

  while(exit)
  {
     G4cout<<"Enter 'y' , if you would like to compare with exp. data : "<<std::flush;
     G4cin>>confirm;
     if(confirm != "y" ) break;
     G4cout<<G4endl;

     // Read experimental data file

     G4double delExp[200], distr[200], deltaBin, sumPAI, sumExp;
     G4int numberOfExpPoints;

     G4cout<<G4endl;
     G4cout << " Enter number of experimental points : " <<std::flush;
     G4cin>>numberOfExpPoints;
     G4cout<<G4endl;
     G4cout << " Enter energy bin in keV : " <<std::flush;
     G4cin>>deltaBin;
     G4cout<<G4endl;
     deltaBin *= keV;

     std::ifstream fileRead;
     fileRead.open("input.dat");
     for(i=0;i<numberOfExpPoints;i++)
     {
       fileRead>>delExp[i]>>distr[i];
       delExp[i] *= keV;
       G4cout<<i<<"\t"<<delExp[i]<<"\t"<<distr[i]<<G4endl;
     }
     fileRead.close();

     // Adjust statistics of experiment to PAI simulation

     sumExp = 0.0;
     for(i=0;i<numberOfExpPoints;i++) sumExp +=distr[i];
     sumExp *= deltaBin;

     sumPAI = 0.0;
     for(i=0;i<49;i++) sumPAI +=spectrum[i];
     sumPAI *= Ebin;

     for(i=0;i<numberOfExpPoints;i++) distr[i] *= sumPAI/sumExp;

     for(i=0;i<numberOfExpPoints;i++)
     {
       fileWrite<<delExp[i]/keV<<"\t"<<distr[i]<<G4endl;
       G4cout<<delExp[i]/keV<<"\t"<<distr[i]<<G4endl;
     }
     exit = 0;
  }

  G4cout<<"Enter 'y' , if you would like to get most probable delta : "<<std::flush;
  G4cin>>confirm;
  if(confirm != "y" ) return 1;
  G4cout<<G4endl;

  G4int kGamma, iMPLoss, maxSpectrum, iMax;
  G4double mpDelta[50], meanDelta[50], rrMP[50], rrMean[50]; 
  G4double mpLoss, tmRatio, mpSum, mpStat;

  G4double aGamma[33] = 
  {
    4.0, 1.5, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0, // 13
    20., 40.0, 60.0, 80.0, 100.0, 200.0, 400.0, 600.0, 800.0, 1000.0, // 23
    2000.0, 4000.0, 6000.0, 8000.0, 100000.0, 20000.0,                // 29
    40000.0, 60000.0, 80000.0, 100000.0                               // 33
  };

  for(k=0;k<numOfMaterials;k++)
  {
    G4cout <<k<< "  Material : " <<(*theMaterialTable)[k]->GetName() << G4endl;
  } 
  G4cout<<"Enter material name for dE/dx-distribution : "<<std::flush;
  G4cin>>testName;
  G4cout<<G4endl;


  for(k=0;k<numOfMaterials;k++)
  {
     if((*theMaterialTable)[k]->GetName() != testName) continue;

     G4cout << "Material : " <<(*theMaterialTable)[k]->GetName() << G4endl<<G4endl;

     G4cout << " Enter nGamma 1<nGamma<10 : "  <<std::flush;
     G4cin>>nGamma;
     G4cout<<G4endl;


     G4cout << " Enter step in mm : " <<std::flush;
     G4cin>>step;
     G4cout<<G4endl;
     step *= mm;

     G4cout << " Enter energy bin in keV : " <<std::flush;
     G4cin>>Ebin;
     G4cout<<G4endl;
     Ebin *= keV;

     G4cout << " Enter trancated mean ration <1.0 : "  <<std::flush;
     G4cin>>tmRatio;
     G4cout<<G4endl;


     G4cout << " Enter number of events : " <<std::flush;
     G4cin>>iStatMax;
     G4cout<<G4endl;

     G4cout<<"no."<<"\t"<<"Gamma"<<"\t"<<"Rel. rise"<<"\t"<<"M.P. loss, keV"
           <<"\t"<<"Mean loss, keV"<<G4endl<<G4endl;
     //   outFile<<"no."<<"\t"<<"Gamma"<<"\t"<<"M.P. loss, keV"
     //      <<"\t"<<"Mean loss, keV"<<G4endl<<G4endl;
     

     // gamma = 1.1852;

     for(kGamma=0;kGamma<33;kGamma++)
     {
       //    G4cout<<G4endl<<"Start dE/dx distribution"<<G4endl<<G4endl;

       gamma = aGamma[kGamma];
       maxEnergyTransfer = 100*keV;
       bg2               = gamma*gamma - 1;
       rateMass          = electron_mass_c2/proton_mass_c2;

       Tmax              = 2.0*electron_mass_c2*bg2
                          /(1.0+2.0*gamma*rateMass+rateMass*rateMass);

       if ( maxEnergyTransfer > Tmax)         maxEnergyTransfer = Tmax;
       
       G4PAIxSection testPAIenergyLoss(k,maxEnergyTransfer,bg2);
 
       for( iLoss = 0; iLoss < 50; iLoss++ )
       {
         energyLoss[iLoss] = Ebin*iLoss;
         spectrum[iLoss] = 0;
       }
       for(iStat=0;iStat<iStatMax;iStat++)
       {

         //   aaa = (G4double)nGamma;
         //   lambda = aaa/step;
         //   colDist = RandGamma::shoot(aaa,lambda);

         //  delta = testPAIenergyLoss.GetStepEnergyLoss(colDist);

         delta = testPAIenergyLoss.GetStepEnergyLoss(step);

         //   delta1 = testPAIenergyLoss.GetStepEnergyLoss(step);
         //   delta2 = testPAIenergyLoss.GetStepEnergyLoss(step);
         //   delta3 = testPAIenergyLoss.GetStepEnergyLoss(step);
 
         //   delta = alphaCrossTalk*delta1 + 
         //         delta2 + alphaCrossTalk*delta3 - betaS;

         for(iLoss=0;iLoss<50;iLoss++)
         {
           if(delta <= energyLoss[iLoss]) break;
         }
         spectrum[iLoss-1]++;
       }
       G4int sumStat = 0;
       for(iLoss=0;iLoss<49;iLoss++) // without last bin
       {
         sumStat += spectrum[iLoss];
         if( sumStat > tmRatio*iStatMax  ) break;
       }
       if(iLoss == 50) iLoss--;
       iMPLoss = iLoss;
       G4double meanLoss = 0.0;
       maxSpectrum = 0;

       for(iLoss=0;iLoss<iMPLoss;iLoss++) // without last bin
       {
	 // fileOut<<energyLoss[iLoss]/keV<<"\t\t"<<spectrum[iLoss]<<G4endl;
	 //  G4cout<<energyLoss[iLoss]/keV<<"\t\t"<<spectrum[iLoss]<<G4endl;

         meanLoss += energyLoss[iLoss]*spectrum[iLoss];

         if( spectrum[iLoss] > maxSpectrum )
	 {
           maxSpectrum = spectrum[iLoss]  ;
           mpLoss      = energyLoss[iLoss];
           iMax = iLoss;
	 }
       }
       mpSum  = 0.;
       mpStat = 0;
       for(iLoss = iMax-5;iLoss<=iMax+5;iLoss++)
       {
         mpSum += energyLoss[iLoss]*spectrum[iLoss];
         mpStat += spectrum[iLoss];
       }
       mpLoss = mpSum/mpStat;
       mpLoss /= keV;
       meanLoss /= keV*sumStat;
       meanDelta[kGamma] = meanLoss;
       mpDelta[kGamma] = mpLoss;

       if(kGamma > 0)
       {
         rrMP[kGamma] = mpLoss/mpDelta[0];
         G4cout<<kGamma<<"\t"<<gamma<<"\t"<<rrMP[kGamma]<<"\t"<<mpLoss<<G4endl;
	 //  outFile<<gamma<<"\t"<<rrMP[kGamma]<<G4endl;
         fileWrite1<<gamma<<"\t"<<rrMP[kGamma]<<G4endl;
       }

       //  gamma *= 1.5;
    }
    G4cout<<G4endl;
    outFile<<G4endl;  
  }   

   return EXIT_SUCCESS;

}










