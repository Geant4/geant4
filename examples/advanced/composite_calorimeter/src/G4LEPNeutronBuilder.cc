#include "G4LEPNeutronBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4LEPNeutronBuilder::
G4LEPNeutronBuilder() 
{
  theMin = 0;
  theMax = 55*GeV;
}

G4LEPNeutronBuilder::
~G4LEPNeutronBuilder() {}

void G4LEPNeutronBuilder::
Build(G4HadronElasticProcess & aP)
{
  theElasticModel = new G4LElastic();
  theElasticModel->SetMinEnergy(theMin);
  aP.RegisterMe(theElasticModel);
}

void G4LEPNeutronBuilder::
Build(G4HadronFissionProcess & aP)
{
  theNeutronFissionModel = new G4LFission();
  theNeutronFissionModel->SetMinEnergy(theMin);
  aP.RegisterMe(theNeutronFissionModel);
}

void G4LEPNeutronBuilder::
Build(G4HadronCaptureProcess & aP)
{
  theNeutronCaptureModel = new G4LCapture();
  theNeutronCaptureModel->SetMinEnergy(theMin);
  aP.RegisterMe(theNeutronCaptureModel);
}

void G4LEPNeutronBuilder::
Build(G4NeutronInelasticProcess & aP)
{
  theLENeutronModel = new G4LENeutronInelastic();
  theLENeutronModel->SetMinEnergy(theIMin);
  theLENeutronModel->SetMaxEnergy(theMax);
  aP.RegisterMe(theLENeutronModel);
}

// 2002 by J.P. Wellisch
