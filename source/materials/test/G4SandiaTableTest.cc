// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SandiaTableTest.cc,v 1.1 1999-11-22 18:16:35 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//
////////////////////////////////////////////////////////////////////////
//
//  This program illustrates the different ways to define photoabsorption 
//  cross section according G4Sandiatable
//
// History:
//
// 15.09.99 V.Grichine, start from G4MaterialTest.cc
//

#include "G4ios.hh"
#include <fstream.h>
#include <iomanip.h>
// #include "g4templates.hh"
#include "globals.hh"
#include "G4Timer.hh"

#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"
#include "G4SandiaTable.hh"

int main() 
{
  // set output format

  G4cout.setf( ios::scientific, ios::floatfield );

  // write results to the file  sandia.out

   ofstream outFile("sandia.out", ios::out ) ;
   outFile.setf( ios::scientific, ios::floatfield );

  G4String name, symbol;             // a=mass of a mole;
  G4double a, z, density;            // z=mean number of protons;  
  G4int iz, n;                       // iz=nb of protons  in an isotope; 
                                     // n=nb of nucleons in an isotope;

  G4int ncomponents, natoms, nel ;
  G4double abundance, fractionmass ;
  G4double temperature, pressure ;

  G4UnitDefinition::BuildUnitsTable();

//
// define Elements
//

  a = 1.01*g/mole;
  G4Element* elH  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a);

  a = 12.01*g/mole;
  G4Element* elC  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);

  a = 14.01*g/mole;
  G4Element* elN  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

  a = 16.00*g/mole;
  G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

  a = 28.09*g/mole;
  G4Element* elSi = new G4Element(name="Silicon",symbol="Si" , z= 14., a);

  a = 55.85*g/mole;
  G4Element* elFe = new G4Element(name="Iron"    ,symbol="Fe", z=26., a);

  a = 131.29*g/mole;
  G4Element* elXe = new G4Element(name="Xenon", symbol="Xe", z=54., a);  

  a = 19.00*g/mole;
  G4Element* elF  = new G4Element(name="Fluorine", symbol="F", z=9., a);

//
// define an Element from isotopes, by relative abundance 
//

  G4Isotope* U5 = new G4Isotope(name="U235", iz=92, n=235, a=235.01*g/mole);
  G4Isotope* U8 = new G4Isotope(name="U238", iz=92, n=238, a=238.03*g/mole);

  G4Element* elU  = new G4Element(name="enriched Uranium", 
                                symbol="U", ncomponents=2);
  elU->AddIsotope(U5, abundance= 90.*perCent);
  elU->AddIsotope(U8, abundance= 10.*perCent);


// G4cout << *(G4Isotope::GetIsotopeTable()) << endl;

// G4cout << *(G4Element::GetElementTable()) << endl;

//
// define simple materials
//

  density = 2.700*g/cm3;
  a = 26.98*g/mole;
  G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);

  density = 1.390*g/cm3;
  a = 39.95*g/mole;
  G4Material* lAr = new G4Material(name="liquidArgon", z=18., a, density);

  density = 8.960*g/cm3;
  a = 63.55*g/mole;
  G4Material* Cu = new G4Material(name="Copper"   , z=29., a, density);

  density = 11.35*g/cm3;
  a = 207.19*g/mole;
  G4Material* Pb = new G4Material(name="Lead "     , z=82., a, density);

//
// define a material from elements.   case 1: chemical molecule
//
 
  density = 1.000*g/cm3;
  G4Material* H2O = new G4Material(name="Water", density, ncomponents=2);
  H2O->AddElement(elH, natoms=2);
  H2O->AddElement(elO, natoms=1);

  density = 1.032*g/cm3;
  G4Material* Sci = new G4Material(name="Scintillator", density, ncomponents=2);
  Sci->AddElement(elC, natoms=9);
  Sci->AddElement(elH, natoms=10);

  density = 2.200*g/cm3;
  G4Material* SiO2 = new G4Material(name="quartz", density, ncomponents=2);
  SiO2->AddElement(elSi, natoms=1);
  SiO2->AddElement(elO , natoms=2);

//
// define a material from elements.   case 2: mixture by fractional mass
//

  density = 1.290*mg/cm3;
  G4Material* Air = new G4Material(name="Air  "  , density, ncomponents=2);
  Air->AddElement(elN, fractionmass=0.7);
  Air->AddElement(elO, fractionmass=0.3);

//
// define a material from elements and/or others materials (mixture of mixtures)
//

  density = 0.200*g/cm3;
  G4Material* Aerog = new G4Material(name="Aerogel", density, ncomponents=3);
  Aerog->AddMaterial(SiO2, fractionmass=0.625);
  Aerog->AddMaterial(H2O , fractionmass=0.374);
  Aerog->AddElement (elC , fractionmass=0.1*perCent);

//
// examples of gas in non STP conditions
//

  density     = 27.*mg/cm3;
  pressure    = 50.*atmosphere;
  temperature = 325.*kelvin;
  G4Material* CO2 = new G4Material(name="Carbonic gas", density, ncomponents=2,
                                     kStateGas,temperature,pressure);
  CO2->AddElement(elC, natoms=1);
  CO2->AddElement(elO, natoms=2);
 
  density     = 0.3*mg/cm3;
  pressure    = 2.*atmosphere;
  temperature = 500.*kelvin;
  G4Material* steam = new G4Material(name="Water steam ", density, ncomponents=1,
                                      kStateGas,temperature,pressure);
  steam->AddMaterial(H2O, fractionmass=1.);

//
// examples of vacuum
//

  density     = universe_mean_density;            //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  new G4Material(name="Galactic", z=1., a=1.01*g/mole, density,
                   kStateGas,temperature,pressure);

  density     = 1.e-5*g/cm3;
  pressure    = 2.e-2*bar;
  temperature = STP_Temperature;                      //from PhysicalConstants.h

  G4Material* beam = new G4Material(name="Beam ", density, ncomponents=1,
                                      kStateGas,temperature,pressure);
  beam->AddMaterial(Air, fractionmass=1.);

  density = 1.39*g/cm3;
  G4Material* Kapton = new G4Material(name="Kapton", density, nel=3);
  Kapton->AddElement(elO,2);
  Kapton->AddElement(elC,5);
  Kapton->AddElement(elH,4);


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

  G4double XeCO2CF4_density = 4.76*mg/cm3;
  G4Material* XeCO2CF4 = new G4Material(name="XeCO2CF4", XeCO2CF4_density,
					ncomponents=3,
					kStateGas,293.15*kelvin,1.*atmosphere);
  XeCO2CF4->AddMaterial(TRT_Xe,0.807);
  XeCO2CF4->AddMaterial(TRT_CO2,0.039);
  XeCO2CF4->AddMaterial(TRT_CF4,0.154);
      
  density = 0.935*g/cm3;
  G4Material* TRT_CH2 = new G4Material(name="TRT_CH2",density, nel=2);
  TRT_CH2->AddElement(elC,1);
  TRT_CH2->AddElement(elH,2);

  density = 0.059*g/cm3;
  G4Material* Radiator = new G4Material(name="Radiator",density, nel=2);
  Radiator->AddElement(elC,1);
  Radiator->AddElement(elH,2);

  density = 0.145*g/cm3;
  G4Material* CarbonFiber = new G4Material(name="CarbonFiber",density, nel=1);
  CarbonFiber->AddElement(elC,1);

//
// Print the table of materials
//

// G4cout << *(G4Material::GetMaterialTable()) << endl;


//
// Checking Sandia table coefficients
//
  G4int numberOfMat, iMat, matIndex, nbOfElements, sanIndex, row, iSan ;
  G4String materialName = "Air" ;
  static const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  numberOfMat = theMaterialTable->length() ;

  for(iMat=0;iMat<numberOfMat;iMat++)
  {
    if(materialName == (*theMaterialTable)[iMat]->GetName() )
    {
      matIndex = (*theMaterialTable)[iMat]->GetIndex() ;
      break ;
    }
  }

  for(iMat=0;iMat<numberOfMat;iMat++)
  {
     G4String matName = (*theMaterialTable)[iMat]->GetName() ;
     matIndex = (*theMaterialTable)[iMat]->GetIndex() ;
     nbOfElements = (*theMaterialTable)[iMat]->GetNumberOfElements() ;

     G4cout<<matIndex<<"\t"<<matName<<endl<<endl ;
     outFile<<matIndex<<"\t"<<matName<<endl<<endl ;
     G4cout<<"Sandia cof according old PAI stuff"<<endl<<endl ;
     outFile<<"Sandia cof according old PAI stuff"<<endl<<endl ;

     G4int* thisMaterialZ = new G4int[nbOfElements] ;
     for(iSan=0;iSan<nbOfElements;iSan++)
     {
        thisMaterialZ[iSan] = (G4int)(*theMaterialTable)[iMat]->
                                      GetElement(iSan)->GetZ() ;
     }
     G4SandiaTable sandia(matIndex) ;
     sanIndex = sandia.SandiaIntervals(thisMaterialZ,nbOfElements) ;    
     sanIndex = sandia.SandiaMixing( thisMaterialZ ,
                             (*theMaterialTable)[iMat]->GetFractionVector() ,
				     nbOfElements,sanIndex) ;

     for(row=0;row<sanIndex-1;row++)
     {
       G4cout<<row+1<<"\t"<<sandia.GetPhotoAbsorpCof(row+1,0)/keV ;
       outFile<<row+1<<"  "<<sandia.GetPhotoAbsorpCof(row+1,0)/keV ;

       for(iSan=1;iSan<5;iSan++)
       {
         G4cout<<"\t"<<sandia.GetPhotoAbsorpCof(row+1,iSan) ;
	 // *(*theMaterialTable)[iMat]->GetDensity() ;

         outFile<<"  "<<sandia.GetPhotoAbsorpCof(row+1,iSan) ;
	 // *(*theMaterialTable)[iMat]->GetDensity() ;
       }
       G4cout<<endl ;
       outFile<<endl ;
     }
     G4cout<<endl ;
     outFile<<endl ;

     G4SandiaTable* sanMatrix = G4Material::GetMaterial(matName)->
                                GetSandiaTable() ;
     sanIndex = sanMatrix->GetMatNbOfIntervals() ;
      
     G4cout<<"Sandia cof according ComputeMixSandiaMatrix()"<<endl<<endl ;
     outFile<<"Sandia cof according ComputeMixSandiaMatrix()"<<endl<<endl ;

     for(row=0;row<sanIndex;row++)
     {
       G4cout<<row+1<<"\t"<<sanMatrix->GetSandiaCofForMaterial(row,0)/keV ;
       outFile<<row+1<<"  "<<sanMatrix->GetSandiaCofForMaterial(row,0)/keV ;

       for(iSan=1;iSan<5;iSan++)
       {
         G4cout<<"\t"<<sanMatrix->GetSandiaCofForMaterial(row,iSan) ;
	 // *(*theMaterialTable)[iMat]->GetDensity() ;
         outFile<<"  "<<sanMatrix->GetSandiaCofForMaterial(row,iSan) ;
	 // *(*theMaterialTable)[iMat]->GetDensity() ;
       }
       G4cout<<endl ;
       outFile<<endl ;
     }      
     G4cout<<endl ;
     outFile<<endl ;
  }




  return EXIT_SUCCESS;
}

//
//
///////////////////////  end of G4SandiaTableTest.cc  /////////////////////////
