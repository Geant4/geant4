#include "G4LHEPNeutronBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4LHEPNeutronBuilder::
G4LHEPNeutronBuilder() 
{
  theMin = 0;
  theIMin = 0;
}

G4LHEPNeutronBuilder::
~G4LHEPNeutronBuilder() {}

void G4LHEPNeutronBuilder::
Build(G4HadronElasticProcess & aP)
{
  theElasticModel = new G4LElastic();
  theElasticModel->SetMinEnergy(theMin);
  aP.RegisterMe(theElasticModel);
}

void G4LHEPNeutronBuilder::
Build(G4HadronFissionProcess & aP)
{
  theNeutronFissionModel = new G4LFission();
  theNeutronFissionModel->SetMinEnergy(theMin);
  aP.RegisterMe(theNeutronFissionModel);
}

void G4LHEPNeutronBuilder::
Build(G4HadronCaptureProcess & aP)
{
  theNeutronCaptureModel = new G4LCapture();
  theNeutronCaptureModel->SetMinEnergy(theMin);
  aP.RegisterMe(theNeutronCaptureModel);
}

void G4LHEPNeutronBuilder::
Build(G4NeutronInelasticProcess & aP)
{
  theLENeutronModel = new G4LENeutronInelastic();
  theHENeutronModel = new G4HENeutronInelastic();
  theLENeutronModel->SetMinEnergy(theIMin);
  theLENeutronModel->SetMaxEnergy(55*GeV);
  theHENeutronModel->SetMinEnergy(25*GeV);
  aP.RegisterMe(theLENeutronModel);
  aP.RegisterMe(theHENeutronModel);
}

// 2002 by J.P. Wellisch
