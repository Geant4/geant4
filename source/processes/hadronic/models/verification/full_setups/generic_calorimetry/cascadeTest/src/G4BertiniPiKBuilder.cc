#include "G4BertiniPiKBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4BertiniPiKBuilder::
G4BertiniPiKBuilder() 
{
  theMin = 0*GeV;
  theMax = 5*GeV;
  theModel = new G4CascadeInterface;
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax); 
}

G4BertiniPiKBuilder::
~G4BertiniPiKBuilder() {}

void G4BertiniPiKBuilder::
Build(G4HadronElasticProcess & aP) {}

void G4BertiniPiKBuilder::
Build(G4PionPlusInelasticProcess & aP)
{
  aP.RegisterMe(theModel);
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax);
}

void G4BertiniPiKBuilder::
Build(G4PionMinusInelasticProcess & aP)
{
  aP.RegisterMe(theModel);
  aP.GetCrossSectionDataStore()->AddDataSet(&thePiData);
  theModel->SetMinEnergy(theMin);
  theModel->SetMaxEnergy(theMax);
}

void G4BertiniPiKBuilder::
Build(G4KaonPlusInelasticProcess & aP)
{
}

void G4BertiniPiKBuilder::
Build(G4KaonMinusInelasticProcess & aP)
{
}

void G4BertiniPiKBuilder::
Build(G4KaonZeroLInelasticProcess & aP)
{
}

void G4BertiniPiKBuilder::
Build(G4KaonZeroSInelasticProcess & aP)
{
}

// 2002 by J.P. Wellisch
