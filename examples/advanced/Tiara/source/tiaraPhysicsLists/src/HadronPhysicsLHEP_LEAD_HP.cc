#include "HadronPhysicsLHEP_LEAD_HP.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/iomanip"   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

HadronPhysicsLHEP_LEAD_HP::HadronPhysicsLHEP_LEAD_HP(const G4String& name)
                    :  G4VPhysicsConstructor(name) 
{
  theNeutrons.RegisterMe(&theLHEPNeutron);
  theLHEPNeutron.SetMinEnergy(19.9*MeV);
  theLHEPNeutron.SetMinInelasticEnergy(4.99*GeV);
  theNeutrons.RegisterMe(&theLEADNeutron);
  theLEADNeutron.SetMinEnergy(19.9*MeV);
  theNeutrons.RegisterMe(&theHPNeutron);

  thePro.RegisterMe(&theLHEPPro);
  theLHEPPro.SetMinEnergy(4.99*GeV);
  thePro.RegisterMe(&theLEADPro);
  
  thePiK.RegisterMe(&theLHEPPiK);
  theLHEPPiK.SetMinEnergy(4.99*GeV);
  thePiK.RegisterMe(&theLEADPiK);
}

HadronPhysicsLHEP_LEAD_HP::~HadronPhysicsLHEP_LEAD_HP() {}

void HadronPhysicsLHEP_LEAD_HP::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsLHEP_LEAD_HP::ConstructProcess()
{
  theNeutrons.Build();
  thePro.Build();
  thePiK.Build();
  theMiscLHEP.Build();
  theStoppingHadron.Build();
  theHadronQED.Build();
}
// 2002 by J.P. Wellisch
