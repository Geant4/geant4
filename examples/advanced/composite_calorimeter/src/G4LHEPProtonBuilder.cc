#include "G4LHEPProtonBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4LHEPProtonBuilder::
G4LHEPProtonBuilder() 
{
  theMin = 0;
}

G4LHEPProtonBuilder::
~G4LHEPProtonBuilder() {}

void G4LHEPProtonBuilder::
Build(G4HadronElasticProcess & aP)
{
  theElasticModel = new G4LElastic();
  aP.RegisterMe(theElasticModel);
}

void G4LHEPProtonBuilder::
Build(G4ProtonInelasticProcess & aP)
{
  theLEProtonModel = new G4LEProtonInelastic();
  theHEProtonModel = new G4HEProtonInelastic();
  theLEProtonModel->SetMinEnergy(theMin);
  theLEProtonModel->SetMaxEnergy(55*GeV);
  theHEProtonModel->SetMinEnergy(25*GeV);
  aP.RegisterMe(theLEProtonModel);
  aP.RegisterMe(theHEProtonModel);
}

// 2002 by J.P. Wellisch
