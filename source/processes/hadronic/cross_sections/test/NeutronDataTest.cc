#include "G4Neutron.hh"
#include "G4NeutronInelasticCrossSection.hh"


main()
{
   G4NeutronInelasticCrossSection aDataSet;
   G4Element* theElement = new G4Element("copper", "Cu", 29, 63.54*g/mole);
   G4ParticleDefinition* theParticleDefinition = G4Neutron::NeutronDefinition();

   G4double ekin = 0.1*MeV;
   G4DynamicParticle* theDynamicParticle;
   while(ekin < 20*GeV)
   {
     ekin *= 1.1;
     theDynamicParticle = new G4DynamicParticle(theParticleDefinition,
                                                 G4ParticleMomentum(1.,0.,0.), ekin);
//     if(aDataSet.IsApplicable(theDynamicParticle, theElement))
     {
       cout << ekin/MeV << " " 
            << aDataSet.GetCrossSection(theDynamicParticle, theElement)/millibarn
	    <<endl;
     }
     delete theDynamicParticle;
   }
}
