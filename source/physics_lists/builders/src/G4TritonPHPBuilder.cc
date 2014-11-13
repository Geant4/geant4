#include "G4TritonPHPBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleHPInelasticData.hh"

G4TritonPHPBuilder::
G4TritonPHPBuilder() 
{
  theMin = 0;
  theMax=200.*MeV;
}

G4TritonPHPBuilder::
~G4TritonPHPBuilder() 
{
  delete theParticlePHPModel;
}

void G4TritonPHPBuilder::
Build(G4HadronElasticProcess *)
{
  G4cout << "Info - G4TritonPHPBuilder::Build() not adding elastic" << G4endl;
}

void G4TritonPHPBuilder::
Build(G4TritonInelasticProcess * aP)
{
  G4ParticleHPInelasticData* theTritonHPInelasticData=new G4ParticleHPInelasticData(G4Triton::Triton());
  theTritonHPInelasticData->SetMinKinEnergy(theMin);
  theTritonHPInelasticData->SetMaxKinEnergy(theMax);
  aP->AddDataSet(theTritonHPInelasticData);

  theParticlePHPModel = new G4ParticleHPInelastic(G4Triton::Triton(),"ParticleHPInelastic");
  theParticlePHPModel->SetMinEnergy(theMin);
  theParticlePHPModel->SetMaxEnergy(theMax);
  aP->RegisterMe(theParticlePHPModel);

}

