#include "G4BertiniNeutronBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4BertiniNeutronBuilder::
G4BertiniNeutronBuilder() 
{
  theMin = 0;
  theMax = 5*GeV;
  theModel = new G4CascadeInterface;
}

G4BertiniNeutronBuilder::
~G4BertiniNeutronBuilder() {}

void G4BertiniNeutronBuilder::
Build(G4HadronElasticProcess & aP)
{
}

void G4BertiniNeutronBuilder::
Build(G4HadronFissionProcess & aP)
{
}

void G4BertiniNeutronBuilder::
Build(G4HadronCaptureProcess & aP)
{
}

void G4BertiniNeutronBuilder::
Build(G4NeutronInelasticProcess & aP)
{
  theModel->SetMinEnergy(0);
  theModel->SetMaxEnergy(theMax);
  aP.RegisterMe(theModel);
  G4CrossSectionDataStore * theNStore;
  theNStore = aP.GetCrossSectionDataStore();
  theNStore->AddDataSet(&theXSec);  
}

// 2002 by J.P. Wellisch
