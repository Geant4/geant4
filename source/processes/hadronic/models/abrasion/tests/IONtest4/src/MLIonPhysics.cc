////////////////////////////////////////////////////////////////////////////////
//
#include "MLIonPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>
////////////////////////////////////////////////////////////////////////////////
//
MLIonPhysics::MLIonPhysics (const G4String& name)
  :G4VPhysicsConstructor(name), mode(name)
{}
////////////////////////////////////////////////////////////////////////////////
//
MLIonPhysics::~MLIonPhysics ()
{}
////////////////////////////////////////////////////////////////////////////////
//
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

// Nuclei
#include "G4IonConstructor.hh"

void MLIonPhysics::ConstructParticle ()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}
////////////////////////////////////////////////////////////////////////////////
//
#include "G4ProcessManager.hh"
#include "G4ionIonisation.hh"
#include "G4hLowEnergyIonisation.hh"
#include "G4WilsonAbrasionModel.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4GeneralSpaceNNCrossSection.hh"
#include "G4ElementTable.hh"
#include "G4MaterialTable.hh"

void MLIonPhysics::ConstructProcess ()
{
  // Elastic Process
  theElasticModel = new G4LElastic();
  theElasticProcess.RegisterMe(theElasticModel);
  
  G4WilsonAbrasionModel *theAM = new G4WilsonAbrasionModel(true);
  theAM->SetVerboseLevel(2);
  theAM->SetMinEnergy(100.0*MeV);
  G4BinaryLightIonReaction * theGenIonBC = new G4BinaryLightIonReaction;
  theGenIonBC->SetMinEnergy(0.0*MeV);
  theGenIonBC->SetMaxEnergy(100.0*MeV);
//  theGenIonBC->SetMaxEnergy(5.0*GeV);
  G4BinaryLightIonReaction * theGenIonBC1 = new G4BinaryLightIonReaction;
  theGenIonBC1->SetMinEnergy(100.0*MeV);
  theGenIonBC1->SetMaxEnergy(5.0*GeV);
  G4ElementTable::iterator iterE;
  G4ElementTable *elementTable =
    const_cast <G4ElementTable*> (G4Element::GetElementTable());
  for (iterE = elementTable->begin(); iterE != elementTable->end(); ++iterE)
    theAM->ActivateFor(*iterE);
  G4MaterialTable::iterator iterM;
  G4MaterialTable *materialTable =
    const_cast <G4MaterialTable*> (G4Material::GetMaterialTable());
  for (iterM = materialTable->begin(); iterM != materialTable->end(); ++iterM)
    theAM->ActivateFor(*iterM);
  
  G4GeneralSpaceNNCrossSection * generalCrossSection =
    new G4GeneralSpaceNNCrossSection;

  theParticleIterator->reset();

  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pManager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     if (   particleName == "alpha"  
               || particleName == "deuteron"  
               || particleName == "tritium"
	       || particleName == "he3"
               || particleName == "GenericIon"  
              )
    {
      pManager->AddDiscreteProcess(&theElasticProcess);
      pManager->AddProcess(new G4MultipleScattering,-1,1,1);
      //
      if (mode == "Ion") {
	G4ionIonisation* iIon = new G4ionIonisation();
	pManager->AddProcess(iIon,-1,2, 2);
      } else {
	G4hLowEnergyIonisation* iIon = new G4hLowEnergyIonisation();
	iIon->SetCutForSecondaryPhotons(100.*eV) ;
	pManager->AddProcess(iIon,-1,2,2);
      }
    }
    if (particleName == "deuteron")
    {
      fDeuteronProcess.AddDataSet(generalCrossSection);
      fDeuteronProcess.RegisterMe(theAM);
//      fDeuteronProcess.RegisterMe(theGenIonBC1);
      fDeuteronModel = new G4LEDeuteronInelastic;
      fDeuteronModel->SetMaxEnergy(100.0*MeV);
      fDeuteronProcess.RegisterMe(fDeuteronModel);
      pManager->AddDiscreteProcess(&fDeuteronProcess);
    }
    else if (particleName == "triton")
    {
      fTritonProcess.AddDataSet(generalCrossSection);
      fTritonProcess.RegisterMe(theAM);
//      fTritonProcess.RegisterMe(theGenIonBC1);
      fTritonModel = new G4LETritonInelastic;
      fTritonModel->SetMaxEnergy(100.0*MeV);
      fTritonProcess.RegisterMe(fTritonModel);
      pManager->AddDiscreteProcess(&fTritonProcess);
    }
    else if (particleName == "alpha")
    {
      fAlphaProcess.AddDataSet(generalCrossSection);
      fAlphaProcess.RegisterMe(theAM);
//      fAlphaProcess.RegisterMe(theGenIonBC1);
      fAlphaModel = new G4LEAlphaInelastic;
      fAlphaModel->SetMaxEnergy(100.0*MeV);
      fAlphaProcess.RegisterMe(fAlphaModel);
      pManager->AddDiscreteProcess(&fTritonProcess);
    }
    else if (particleName == "GenericIon")
    {
      G4HadronInelasticProcess* fGenericIon =
        new G4HadronInelasticProcess("IonInelastic",particle);
      fGenericIon->AddDataSet(generalCrossSection);
      fGenericIon->RegisterMe(theAM);
      fGenericIon->RegisterMe(theGenIonBC);
      pManager->AddDiscreteProcess(fGenericIon);
    }
    else if (particleName == "he3")
    {
      G4HadronInelasticProcess* fhe3Ion =
        new G4HadronInelasticProcess("He3Inelastic",particle);
      fhe3Ion->AddDataSet(generalCrossSection);
      fhe3Ion->RegisterMe(theAM);
      fhe3Ion->RegisterMe(theGenIonBC);
      pManager->AddDiscreteProcess(fhe3Ion);
    }
//##############################################################################
  }
}
////////////////////////////////////////////////////////////////////////////////
