#include "G4NeutronHPBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

G4NeutronHPBuilder::
G4NeutronHPBuilder() 
{
  theHPElastic = NULL;
  theHPElasticData = NULL;
  
  theHPFission = NULL;
  theHPFissionData = NULL;
  
  theHPCapture = NULL;
  theHPCaptureData = NULL;
  
  theHPInelastic = NULL;
  theHPInelasticData = NULL;
}

G4NeutronHPBuilder::
~G4NeutronHPBuilder() 
{
  delete theHPElasticData;
  delete theHPFissionData;
  delete theHPCaptureData;
  delete theHPInelasticData;
}

void G4NeutronHPBuilder::
Build(G4HadronElasticProcess & aP)
{
  if(theHPElastic==NULL) theHPElastic = new G4NeutronHPElastic;
  if(theHPElasticData == NULL) theHPElasticData = new G4NeutronHPElasticData;
  G4CrossSectionDataStore * theNStore;
  theNStore = aP.GetCrossSectionDataStore();
  theNStore->AddDataSet(theHPElasticData);
  aP.RegisterMe(theHPElastic);
}

void G4NeutronHPBuilder::
Build(G4HadronFissionProcess & aP)
{
  if(theHPFission == NULL) theHPFission = new G4NeutronHPFission;
  if(theHPFissionData==NULL) theHPFissionData=new G4NeutronHPFissionData;
  G4CrossSectionDataStore * theNStore;
  theNStore = aP.GetCrossSectionDataStore();
  theNStore->AddDataSet(theHPFissionData);
  aP.RegisterMe(theHPFission);
}

void G4NeutronHPBuilder::
Build(G4HadronCaptureProcess & aP)
{
  if(theHPCapture==NULL) theHPCapture = new G4NeutronHPCapture;
  if(theHPCaptureData==NULL) theHPCaptureData = new G4NeutronHPCaptureData;
  G4CrossSectionDataStore * theNStore;
  theNStore = aP.GetCrossSectionDataStore();
  theNStore->AddDataSet(theHPCaptureData);
  aP.RegisterMe(theHPCapture);
}

void G4NeutronHPBuilder::
Build(G4NeutronInelasticProcess & aP)
{
  if(theHPInelastic==NULL) theHPInelastic = new G4NeutronHPInelastic;
  if(theHPInelasticData==NULL) theHPInelasticData = new G4NeutronHPInelasticData;
  G4CrossSectionDataStore * theNStore;
  theNStore = aP.GetCrossSectionDataStore();
  theNStore->AddDataSet(theHPInelasticData);
  aP.RegisterMe(theHPInelastic);
}
// 2002 by J.P. Wellisch
