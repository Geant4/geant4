#include "HadronPhysicsQGSP_BIC.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

HadronPhysicsQGSP_BIC::HadronPhysicsQGSP_BIC(const G4String& name)
                    :  G4VPhysicsConstructor(name) 
{
  theNeutrons.RegisterMe(&theQGSPNeutron);
  theNeutrons.RegisterMe(&theLEPNeutron);
  theNeutrons.RegisterMe(&theBinaryNeutron);
  theLEPNeutron.SetMaxEnergy(25*GeV);
  theLEPNeutron.SetMinInelasticEnergy(2.8*GeV);
  theBinaryNeutron.SetMaxEnergy(3.2*GeV);

  thePro.RegisterMe(&theQGSPPro);
  thePro.RegisterMe(&theLEPPro);
  thePro.RegisterMe(&theBinaryPro);
  theLEPPro.SetMaxEnergy(25*GeV);
  theLEPPro.SetMinEnergy(2.8*GeV);
  theBinaryPro.SetMaxEnergy(3.2*GeV);
  
  thePiK.RegisterMe(&theQGSPPiK);
  thePiK.RegisterMe(&theLEPPiK);
  theLEPPiK.SetMaxEnergy(25*GeV);
  
}

HadronPhysicsQGSP_BIC::~HadronPhysicsQGSP_BIC() {}

void HadronPhysicsQGSP_BIC::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsQGSP_BIC::ConstructProcess()
{
  theNeutrons.Build();
  thePro.Build();
  thePiK.Build();
  theMiscLHEP.Build();
  theStoppingHadron.Build();
  theHadronQED.Build();
}
// 2002 by J.P. Wellisch
