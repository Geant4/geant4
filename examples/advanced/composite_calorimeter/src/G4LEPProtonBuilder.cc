#include "G4LEPProtonBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4LEPProtonBuilder::
G4LEPProtonBuilder() 
{
  theMin = 0;
  theMax=55*GeV;
}

G4LEPProtonBuilder::
~G4LEPProtonBuilder() {}

void G4LEPProtonBuilder::
Build(G4HadronElasticProcess & aP)
{
  theElasticModel = new G4LElastic();
  aP.RegisterMe(theElasticModel);
}

void G4LEPProtonBuilder::
Build(G4ProtonInelasticProcess & aP)
{
  theLEProtonModel = new G4LEProtonInelastic();
  theLEProtonModel->SetMinEnergy(theMin);
  theLEProtonModel->SetMaxEnergy(theMax);
  aP.RegisterMe(theLEProtonModel);
}

// 2002 by J.P. Wellisch
