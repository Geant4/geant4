#include "HadronPhysicsLHEP_PRECO.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

HadronPhysicsLHEP_PRECO::HadronPhysicsLHEP_PRECO(const G4String& name)
                    :  G4VPhysicsConstructor(name) 
{
  theNeutrons.RegisterMe(&thePrecoNeutron);
  theNeutrons.RegisterMe(&theLHEPNeutron);
  theLHEPNeutron.SetMinInelasticEnergy(0.15*GeV);

  thePro.RegisterMe(&thePrecoPro);
  thePro.RegisterMe(&theLHEPPro);
  theLHEPPro.SetMinEnergy(0.15*GeV);
  
  thePiK.RegisterMe(&theLHEPPiK);
}

HadronPhysicsLHEP_PRECO::~HadronPhysicsLHEP_PRECO() {}

void HadronPhysicsLHEP_PRECO::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsLHEP_PRECO::ConstructProcess()
{
  theNeutrons.Build();
  thePro.Build();
  thePiK.Build();
  theMiscLHEP.Build();
  theStoppingHadron.Build();
  theHadronQED.Build();
}
// 2002 by J.P. Wellisch
