#include "G4LEADNeutronBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4LEADNeutronBuilder::
G4LEADNeutronBuilder() 
{
  theMin = 0;
  theModel = new G4Mars5GeV;
}

G4LEADNeutronBuilder::
~G4LEADNeutronBuilder() {}

void G4LEADNeutronBuilder::
Build(G4HadronElasticProcess & aP)
{
}

void G4LEADNeutronBuilder::
Build(G4HadronFissionProcess & aP)
{
}

void G4LEADNeutronBuilder::
Build(G4HadronCaptureProcess & aP)
{
}

void G4LEADNeutronBuilder::
Build(G4NeutronInelasticProcess & aP)
{
  theModel->SetMinEnergy(theMin);
  aP.RegisterMe(theModel);
  G4CrossSectionDataStore * theNStore;
  theNStore = aP.GetCrossSectionDataStore();
  theNStore->AddDataSet(&theXSec);  
}

// 2002 by J.P. Wellisch
