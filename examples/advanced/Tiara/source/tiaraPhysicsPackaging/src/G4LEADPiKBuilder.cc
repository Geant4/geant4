#include "G4LEADPiKBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4LEADPiKBuilder::
G4LEADPiKBuilder() 
{
  theMin = 0;
  theModel = new G4Mars5GeV;
}

G4LEADPiKBuilder::
~G4LEADPiKBuilder() {}

void G4LEADPiKBuilder::
Build(G4HadronElasticProcess & aP) {}

void G4LEADPiKBuilder::
Build(G4PionPlusInelasticProcess & aP)
{
  aP.RegisterMe(theModel);
  theModel->SetMinEnergy(theMin);
}

void G4LEADPiKBuilder::
Build(G4PionMinusInelasticProcess & aP)
{
  aP.RegisterMe(theModel);
  aP.GetCrossSectionDataStore()->AddDataSet(&thePiData);
  theModel->SetMinEnergy(theMin);
}

void G4LEADPiKBuilder::
Build(G4KaonPlusInelasticProcess & aP)
{
  aP.RegisterMe(theModel);
  aP.GetCrossSectionDataStore()->AddDataSet(&thePiData);
  theModel->SetMinEnergy(theMin);
}

void G4LEADPiKBuilder::
Build(G4KaonMinusInelasticProcess & aP)
{
  aP.RegisterMe(theModel);
  theModel->SetMinEnergy(theMin);
}

void G4LEADPiKBuilder::
Build(G4KaonZeroLInelasticProcess & aP)
{
  aP.RegisterMe(theModel);
  theModel->SetMinEnergy(theMin);
}

void G4LEADPiKBuilder::
Build(G4KaonZeroSInelasticProcess & aP)
{
  aP.RegisterMe(theModel);
  theModel->SetMinEnergy(theMin);
}

// 2002 by J.P. Wellisch
