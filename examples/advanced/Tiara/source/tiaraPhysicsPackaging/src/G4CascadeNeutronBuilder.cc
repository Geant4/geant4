#include "G4CascadeNeutronBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4LCapture.hh"
#include "G4LFission.hh"

G4CascadeNeutronBuilder::
G4CascadeNeutronBuilder() 
{
  theMin = 0;
  theMax = 150.*MeV;
  theModel = new G4BinaryCascade();
}

G4CascadeNeutronBuilder::
~G4CascadeNeutronBuilder() {}

void G4CascadeNeutronBuilder::
Build(G4HadronElasticProcess & aP)
{
}

void G4CascadeNeutronBuilder::
Build(G4HadronFissionProcess & aP)
{
}

void G4CascadeNeutronBuilder::
Build(G4HadronCaptureProcess & aP)
{
}

void G4CascadeNeutronBuilder::
Build(G4NeutronInelasticProcess & aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax);
  aP.RegisterMe(theModel);
  G4CrossSectionDataStore * theNStore;
  theNStore = aP.GetCrossSectionDataStore();
  theNStore->AddDataSet(&theXSec);  
}

// 2002 by J.P. Wellisch
