#include "G4FTFPProtonBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4FTFPProtonBuilder::
G4FTFPProtonBuilder() 
{
  theMin = 15*GeV;
  theModel = new G4TheoFSGenerator;
  theCascade = new G4GeneratorPrecompoundInterface;
  thePreEquilib = new G4PreCompoundModel(&theHandler);
  theCascade->SetDeExcitation(thePreEquilib);  
  theModel->SetTransport(theCascade);
  theModel->SetHighEnergyGenerator(&theStringModel);
  theStringDecay = new G4ExcitedStringDecay(&theFragmentation);
  theStringModel.SetFragmentationModel(theStringDecay);
}

G4FTFPProtonBuilder::
~G4FTFPProtonBuilder() 
{
  delete theStringDecay;
}

void G4FTFPProtonBuilder::
Build(G4HadronElasticProcess & aP)
{
}

void G4FTFPProtonBuilder::
Build(G4ProtonInelasticProcess & aP)
{
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(100*TeV);
  aP.RegisterMe(theModel);
  G4CrossSectionDataStore * thePStore;
  thePStore = aP.GetCrossSectionDataStore();
  thePStore->AddDataSet(&theXSec);  
}

// 2002 by J.P. Wellisch
