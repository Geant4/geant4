#include "HadronPhysicsLHEP_BIC_HP.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/iomanip"   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

HadronPhysicsLHEP_BIC_HP::HadronPhysicsLHEP_BIC_HP(const G4String& name)
                    :  G4VPhysicsConstructor(name) 
{

  theNeutrons.RegisterMe(&theLHEPNeutron);
  theNeutrons.RegisterMe(&theBinaryNeutron);
  theNeutrons.RegisterMe(&theHPNeutron);
  theLHEPNeutron.SetMinEnergy(19.9*MeV); // added
  theLHEPNeutron.SetMinInelasticEnergy(2.8*GeV);
  theBinaryNeutron.SetMinEnergy(19.9*MeV); // added
  theBinaryNeutron.SetMaxEnergy(3.2*GeV); 

  
  thePro.RegisterMe(&theLHEPPro);
  thePro.RegisterMe(&theBinaryPro);
  theLHEPPro.SetMinEnergy(2.8*GeV); 
  theBinaryPro.SetMaxEnergy(3.2*GeV);
  
  thePiK.RegisterMe(&theLHEPPiK);

}

HadronPhysicsLHEP_BIC_HP::~HadronPhysicsLHEP_BIC_HP() {}

void HadronPhysicsLHEP_BIC_HP::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsLHEP_BIC_HP::ConstructProcess()
{
  theNeutrons.Build();
  thePro.Build();
  thePiK.Build();
  theMiscLHEP.Build();
  theStoppingHadron.Build();
  theHadronQED.Build();
}
// 2002 by J.P. Wellisch
