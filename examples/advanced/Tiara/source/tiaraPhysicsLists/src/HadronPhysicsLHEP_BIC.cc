#include "HadronPhysicsLHEP_BIC.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/iomanip"   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

HadronPhysicsLHEP_BIC::HadronPhysicsLHEP_BIC(const G4String& name)
                    :  G4VPhysicsConstructor(name) 
{
  theNeutrons.RegisterMe(&theLHEPNeutron);
  theNeutrons.RegisterMe(&theBinaryNeutron);
  theLHEPNeutron.SetMinInelasticEnergy(2.8*GeV);
  theBinaryNeutron.SetMaxEnergy(3.2*GeV);

  thePro.RegisterMe(&theLHEPPro);
  thePro.RegisterMe(&theBinaryPro);
  theLHEPPro.SetMinEnergy(2.8*GeV);
  theBinaryPro.SetMaxEnergy(3.2*GeV);
  
  thePiK.RegisterMe(&theLHEPPiK);
  
}

HadronPhysicsLHEP_BIC::~HadronPhysicsLHEP_BIC() {}

void HadronPhysicsLHEP_BIC::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsLHEP_BIC::ConstructProcess()
{
  theNeutrons.Build();
  thePro.Build();
  thePiK.Build();
  theMiscLHEP.Build();
  theStoppingHadron.Build();
  theHadronQED.Build();
}
// 2002 by J.P. Wellisch
