#include "G4Proton.hh"
#include "G4TripathiCrossSection.hh"
#include "G4ParticleTable.hh"
#include "G4IonConstructor.hh"
#include "G4IonProtonCrossSection.hh"

main()
{
   G4IonConstructor theIons;
   G4IonProtonCrossSection aCrossSection;
   theIons.ConstructParticle();
   G4TripathiCrossSection theIonDataSet;
   cout << "Please select the target"<<endl;
   cout << "1: Aluminum"<<endl;
   cout << "2: Calcium"<<endl;
   cout << "3: Iron"<<endl;
   cout << "4: Zinc"<<endl;
   cout << "5: Cobalt"<<endl;
   cout << "6: Nickel"<<endl;
   cout << "7: Lead"<<endl;
   cout << "8: Silver"<<endl;
   cout << "9: Bismuth"<<endl;
   cout << "10: Uranium"<<endl;
   cout << "11: Silicon"<<endl;
   G4int iEle;
   cin >> iEle;
   G4Element* theElement;
   if(iEle == 1) theElement = new G4Element("Aluminum", "Al", 13, 27.0*g/mole);
   if(iEle == 2) theElement = new G4Element("Calcium", "Ca", 20, 40.0*g/mole);
   if(iEle == 3) theElement = new G4Element("Iron", "Fe", 26, 56.0*g/mole);
   if(iEle == 4) theElement = new G4Element("Zinc", "Zn", 30, 64.0*g/mole);
   if(iEle == 5) theElement = new G4Element("Cobalt", "Co", 27, 59.0*g/mole);
   if(iEle == 6) theElement = new G4Element("Nickel", "Ni", 28, 58.0*g/mole);
   if(iEle == 7) theElement = new G4Element("Lead", "Pb", 82, 208.0*g/mole);
   if(iEle == 8) theElement = new G4Element("Silver", "Ag", 47, 109.0*g/mole);
   if(iEle == 9) theElement = new G4Element("Bismuth", "Bi", 83, 209.0*g/mole);
   if(iEle == 10) theElement = new G4Element("Uranium", "U", 92, 238.0*g/mole);
   if(iEle == 11) theElement = new G4Element("Silicon", "Si", 14, 28.0*g/mole);
   G4double Z = 6;
   G4double A = 12;
   cout << "Please enter the projectile Z, A"<<endl;
   cin >> Z >> A;
   G4ParticleDefinition* theParticleDefinition =  
           G4ParticleTable::GetParticleTable()->FindIon(Z, A, 0, Z);

   G4double ekin = 1*A*MeV;
   G4int it;
   cout << "Full chain (1) or fixed input energy (2)?"<<endl;
   cin >> it;
   if(2==it)
   {
     cout << "Please enter energy/nucleon[AMeV]"<<endl;
     cin >> ekin;
     ekin *= A;
   }
   G4DynamicParticle* theDynamicParticle;
   while(ekin < 0.6*A*GeV)
   {
     ekin *= 1.1;
     theDynamicParticle = new G4DynamicParticle(theParticleDefinition,
                                                 G4ParticleMomentum(1.,0.,0.), ekin);
     if(theIonDataSet.IsApplicable(theDynamicParticle, theElement))
     {
       cout << "Token"<< ekin/GeV/A
            << " " 
            << theIonDataSet.GetCrossSection(theDynamicParticle, theElement)/millibarn
            << G4endl;
     }
     delete theDynamicParticle;
   }
}
