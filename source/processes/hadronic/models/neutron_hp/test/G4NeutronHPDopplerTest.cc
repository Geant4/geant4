#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4Alpha.hh"

#include "G4LeptonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4IonConstructor.hh"

main()
{
 G4LeptonConstructor aC1;
 G4BaryonConstructor aC2;
 G4MesonConstructor aC3;
 G4IonConstructor aC4;
 
 aC1.ConstructParticle();
 aC2.ConstructParticle();
 aC3.ConstructParticle();
 aC4.ConstructParticle();

   G4NeutronHPCaptureData genDataSet;

//   G4Element* theElement = new G4Element("copper", "Cu", 29, 63.54*g/mole);
//   G4Element* theElement = new G4Element("copper", "Al", 13, 27.0*g/mole);
//   G4Element* theElement = new G4Element("carbon", "C", 6, 12.0*g/mole);
//   G4Element* theElement = new G4Element("Lead", "Pb", 82, 207.2*g/mole);
//   G4Element* theElement = new G4Element("Hydrogen", "H", 1, 1.01*g/mole);
   G4Element* theElement = new G4Element("Magnesium", "Mg", 12, 24.01*g/mole);
   G4ParticleDefinition * theParticleDefinition;
   theParticleDefinition = G4Neutron::NeutronDefinition();
   genDataSet.BuildPhysicsTable(*theParticleDefinition);
   @@@@@@@@@@@@@@@ Incomplete at the moment @@@@@@@@@@@@@@@@@@

   G4double ekin = 10E-5*eV;
   G4DynamicParticle* theDynamicParticle;
   while(ekin < 20*MeV)
   {
     ekin *= 1.01;
     theDynamicParticle = new G4DynamicParticle(theParticleDefinition,
                                                 G4ParticleMomentum(1.,0.,0.), ekin);
       cout << ekin/MeV 
            << " " 
            << genDataSet.GetCrossSection(theDynamicParticle, theElement)/millibarn
            << G4endl;
     delete theDynamicParticle;
   }
}
