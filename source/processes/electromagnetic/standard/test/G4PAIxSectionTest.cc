// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PAIxSectionTest.cc,v 1.1 1999-01-08 16:32:25 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//
//  
//
//  Test routine for G4PAIxSection class code
//

#include "G4ios.hh"
#include <fstream.h>
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4PAIxSection.hh"

int main()
{
   ofstream outFile("EnergyLoss.cc", ios::out ) ;
   outFile.setf( ios::scientific, ios::floatfield );

// -------------------------- Create materials -------------------------------  
   

  G4int iz , n,  nel ;
  G4double a, z, ez, density , temperature, pressure ;
  G4State state ;
  G4String name, symbol ;


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
 
// G4Isotope::DumpInfo();
// G4Element::DumpInfo();

 density = 8.96*g/cm3;

a = 63.55*g/mole;
G4Material* mat = new G4Material(name="Copper", z=29., a, density);

density = 1.00*g/cm3;
mat = new G4Material(name="Water", density, nel=2);
mat->AddElement(elH, 2);
mat->AddElement(elO, 1);

density = 9.43*g/cm3;
mat = new G4Material(name="Ooblium", density, nel=2);
mat->AddElement(elH, .3);
mat->AddElement(elO, .7);

// G4Element::DumpInfo();
// G4Material::DumpInfo();

density = 8.96*g/cm3;
a = 63.55*g/mole;
mat = new G4Material(name="Vacuum1", z=29., a, density, kVacuum);

density = 1.00*g/cm3;
mat = new G4Material(name="Vacuum2", density, nel=2, kVacuum);
mat->AddElement(elH, 2);
mat->AddElement(elO, 1);


density = 5.85e-3*g/cm3 ;
a = 131.29*g/mole ;
mat = new G4Material(name="Xenon", z=54., a, density);

density = 1.782e-3*g/cm3 ;
a = 39.948*g/mole ;
mat = new G4Material(name="Argon", z=18., a, density);


  a = 9.012*g/mole;
  density = 1.848*g/cm3;
  G4Material* Be = new G4Material(name="Beryllium", z=4. , a, density);

  a = 26.98*g/mole;
  density = 2.7*g/cm3;
  G4Material* Al = new G4Material(name="Aluminium", z=13., a, density);

  a = 28.09*g/mole;
  density = 2.33*g/cm3;
  G4Material* Si = new G4Material(name="Silicon", z=14., a, density);

  // G4Element*   elH = new G4Element ("Hydrogen", "H", 1. ,  1.01*g/mole);

  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", ez=7., a);
  a = 16.00*g/mole;
  // G4Element* elO = new G4Element(name="Oxigen", symbol="O", ez=8., a);
  density = 1.29e-03*g/cm3;
  state = kStateGas ;
  temperature = 273.*kelvin ;
  pressure = 1.*atmosphere ;
  G4Material* Air = new G4Material(name="Air", density, nel=2 ,
                                   state ,temperature , pressure ) ;
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);


// G4Element::DumpInfo();
// G4Material::DumpInfo();

// ------------------------ Create PAITable for given material --------------------

G4int i, j, k, numOfMaterials ;
const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable() ;
numOfMaterials = theMaterialTable->length();

for(k=0;k<numOfMaterials;k++)
{
   outFile<<endl ;
   G4double maxEnergyTransfer = 100*keV ;
   G4PAIxSection testPAI(k,maxEnergyTransfer) ;
   outFile << "Material is " <<(*theMaterialTable)[k]->GetName() << endl ;
   outFile<<"Actual spline size = "<<testPAI.GetSplineSize()<<endl ;
   outFile<<"Normalization Cof = "<<testPAI.GetNormalizationCof()<<endl ;
   outFile<<endl ;
   G4cout << "Material is " <<(*theMaterialTable)[k]->GetName() << endl ;
   G4cout << "Actual spline size = "<<testPAI.GetSplineSize()<<endl ;
   G4cout << endl ;
   outFile<<"Lorentz factor"<<"\t"<<"Max E transfer, Mev"<<"\t"
          <<"dE/dx, Mev/mm"<<"\t\t"<<"N, 1/mm"<<endl<<endl ;
   G4double kineticEnergy = 100*MeV ;
      G4PAIxSection testPAIproton(k,maxEnergyTransfer) ;
   // for(j=1;j<20;j++)
   for(j=1;j<testPAIproton.GetNumberOfGammas();j++)
   {
      G4double tau = kineticEnergy/proton_mass_c2 ;
      G4double gamma = tau +1.0 ;
      G4double bg2 = tau*(tau + 2.0) ;
      G4double beta2 = bg2/(gamma*gamma) ;
      G4double rateMass = electron_mass_c2/proton_mass_c2 ;
      G4double Tmax = 2.0*electron_mass_c2*bg2
                   /(1.0+2.0*gamma*rateMass+rateMass*rateMass) ;

      if ( maxEnergyTransfer > Tmax)         
      {
          maxEnergyTransfer = Tmax ;
      }
      // G4PAIxSection testPAIproton(k,maxEnergyTransfer,bg2) ;
      
      //outFile<<gamma<<"\t"
      //       <<maxEnergyTransfer<<"\t\t"
      //       <<testPAIproton.GetMeanEnergyLoss()<<"\t\t"
      //       <<testPAIproton.GetIntegralPAIxSection(1)<<"\t\t"<<endl ;
      outFile<<testPAIproton.GetLorentzFactor(j)<<"\t"
             <<maxEnergyTransfer<<"\t\t"
             <<testPAIproton.GetPAItable(0,j)<<"\t\t"
      	     <<testPAIproton.GetPAItable(1,j)<<"\t\t"<<endl ;
      kineticEnergy *= 2.0 ;
   }
}

return EXIT_SUCCESS;
}


