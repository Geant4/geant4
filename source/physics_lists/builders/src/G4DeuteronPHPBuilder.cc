#include "G4DeuteronPHPBuilder.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleHPInelasticData.hh"

G4DeuteronPHPBuilder::
G4DeuteronPHPBuilder() 
{
  theMin = 0;
  theMax=200.*MeV;
}

G4DeuteronPHPBuilder::
~G4DeuteronPHPBuilder() 
{
  delete theParticlePHPModel;
}

void G4DeuteronPHPBuilder::
Build(G4HadronElasticProcess *)
{
  G4cout << "Info - G4DeuteronPHPBuilder::Build() not adding elastic" << G4endl;
}

void G4DeuteronPHPBuilder::
Build(G4DeuteronInelasticProcess * aP)
{
  G4ParticleHPInelasticData* theDeuteronHPInelasticData=new G4ParticleHPInelasticData(G4Deuteron::Deuteron());
  theDeuteronHPInelasticData->SetMinKinEnergy(theMin);
  theDeuteronHPInelasticData->SetMaxKinEnergy(theMax);
  aP->AddDataSet(theDeuteronHPInelasticData);

  theParticlePHPModel = new G4ParticleHPInelastic(G4Deuteron::Deuteron(),"ParticleHPInelastic");
  theParticlePHPModel->SetMinEnergy(theMin);
  theParticlePHPModel->SetMaxEnergy(theMax);
  aP->RegisterMe(theParticlePHPModel);

}

