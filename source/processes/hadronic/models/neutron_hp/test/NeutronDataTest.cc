#include "G4Neutron.hh"
#include "G4NeutronHPElasticData.hh"


main()
{
   G4Element* theElement = new G4Element("carbon", "C", 6, 12.0*g/mole);
   G4NeutronHPElasticData anElasticDataSet; 
   G4ParticleDefinition* theParticleDefinition = G4Neutron::NeutronDefinition();

   G4double temp;
   G4cin >> temp;
   temp *= kelvin;
//   G4double ekin = 0.5E+04*eV;
   G4double ekin = 1.E-5*eV;
   G4DynamicParticle* theDynamicParticle;
   G4int counter = -1;
   G4int points = -1; 
   while (++points<30)
   {
//     ekin +=0.002E+04*eV;
     ekin *=3.333333333;
     theDynamicParticle = new G4DynamicParticle(theParticleDefinition,
                                              G4ParticleMomentum(1.,0.,0.), ekin);
     while(++counter<100)
     {
       cout << ekin/MeV 
            << " " 
            << anElasticDataSet.GetCrossSection(theDynamicParticle, theElement, temp)/millibarn
            << G4endl;
     }
   delete theDynamicParticle;
   counter = -1;
   }
}
