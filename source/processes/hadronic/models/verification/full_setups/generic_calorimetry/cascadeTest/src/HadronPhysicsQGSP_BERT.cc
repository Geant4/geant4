#include "HadronPhysicsQGSP_BERT.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

HadronPhysicsQGSP_BERT::HadronPhysicsQGSP_BERT(const G4String& name)
                    :  G4VPhysicsConstructor(name) 
{
  theNeutrons.RegisterMe(&theQGSPNeutron);
  theNeutrons.RegisterMe(&theLEPNeutron);
  theNeutrons.RegisterMe(&theBertiniNeutron);
  theLEPNeutron.SetMaxEnergy(25*GeV);
  theLEPNeutron.SetMinInelasticEnergy(2.8*GeV);
  theBertiniNeutron.SetMaxEnergy(3.2*GeV);
  theBertiniNeutron.SetMinEnergy(0.0*GeV);

  thePro.RegisterMe(&theQGSPPro);
  thePro.RegisterMe(&theLEPPro);
  thePro.RegisterMe(&theBertiniPro);
  theLEPPro.SetMaxEnergy(25*GeV);
  theLEPPro.SetMinEnergy(2.8*GeV);
  theBertiniPro.SetMaxEnergy(3.2*GeV);
  
  thePiK.RegisterMe(&theQGSPPiK);
  thePiK.RegisterMe(&theLEPPiK);
  thePiK.RegisterMe(&theBertiniPiK);
  theLEPPiK.SetMaxEnergy(25*GeV);
  theLEPPiK.SetMinPionEnergy(2.8*GeV);
  theBertiniPiK.SetMaxEnergy(3.2*GeV);
  
}

HadronPhysicsQGSP_BERT::~HadronPhysicsQGSP_BERT() {}

void HadronPhysicsQGSP_BERT::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsQGSP_BERT::ConstructProcess()
{
  theNeutrons.Build();
  thePro.Build();
  thePiK.Build();
  theMiscLHEP.Build();
  theStoppingHadron.Build();
  theHadronQED.Build();
}
// 2002 by J.P. Wellisch
