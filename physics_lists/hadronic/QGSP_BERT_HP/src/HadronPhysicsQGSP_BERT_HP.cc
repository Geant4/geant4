#include "HadronPhysicsQGSP_BERT_HP.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

HadronPhysicsQGSP_BERT_HP::HadronPhysicsQGSP_BERT_HP(const G4String& name)
                    :  G4VPhysicsConstructor(name) 
{
  theNeutrons.RegisterMe(&theQGSPNeutron);
  theNeutrons.RegisterMe(&theLEPNeutron);
  theNeutrons.RegisterMe(&theBertiniNeutron);
  theNeutrons.RegisterMe(&theHPNeutron);
  theLEPNeutron.SetMaxInelasticEnergy(25*GeV);
  theLEPNeutron.SetMinEnergy(19.9*MeV);
  theLEPNeutron.SetMinInelasticEnergy(9.5*GeV);
  theBertiniNeutron.SetMaxEnergy(9.9*GeV);
  theBertiniNeutron.SetMinEnergy(19.9*MeV);

  thePro.RegisterMe(&theQGSPPro);
  thePro.RegisterMe(&theLEPPro);
  thePro.RegisterMe(&theBertiniPro);
  theLEPPro.SetMaxEnergy(25*GeV);
  theLEPPro.SetMinEnergy(9.5*GeV);
  theBertiniPro.SetMaxEnergy(9.9*GeV);
  
  thePiK.RegisterMe(&theQGSPPiK);
  thePiK.RegisterMe(&theLEPPiK);
  thePiK.RegisterMe(&theBertiniPiK);
  theLEPPiK.SetMaxEnergy(25*GeV);
  theLEPPiK.SetMinPionEnergy(9.5*GeV);
  theBertiniPiK.SetMaxEnergy(9.9*GeV);
  
}

HadronPhysicsQGSP_BERT_HP::~HadronPhysicsQGSP_BERT_HP() {}

void HadronPhysicsQGSP_BERT_HP::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsQGSP_BERT_HP::ConstructProcess()
{
  theNeutrons.Build();
  thePro.Build();
  thePiK.Build();
  theMiscLHEP.Build();
  theStoppingHadron.Build();
  theHadronQED.Build();
}
// 2002 by J.P. Wellisch
