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
#include "G4EMDissociation.hh"
#include "G4EMDissociationCrossSection.hh"

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
  
  G4ElementTable::iterator iter;
  G4ElementTable *elementTable =
    const_cast<G4ElementTable*>(G4Element::GetElementTable());
  for (iter = elementTable->begin(); iter != elementTable->end(); ++iter)
  {
    theAM->ActivateFor(*iter);
  }
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
      fDeuteronModel = new G4LEDeuteronInelastic;
      fDeuteronModel->SetMaxEnergy(100.0*MeV);
      fDeuteronProcess.RegisterMe(fDeuteronModel);
      pManager->AddDiscreteProcess(&fDeuteronProcess);
    }
    else if (particleName == "triton")
    {
      fTritonProcess.AddDataSet(generalCrossSection);
      fTritonProcess.RegisterMe(theAM);
      fTritonModel = new G4LETritonInelastic;
      fTritonModel->SetMaxEnergy(100.0*MeV);
      fTritonProcess.RegisterMe(fTritonModel);
      pManager->AddDiscreteProcess(&fTritonProcess);
    }
    else if (particleName == "alpha")
    {
      fAlphaProcess.AddDataSet(generalCrossSection);
      fAlphaProcess.RegisterMe(theAM);
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
      fGenericIon = new G4HadronInelasticProcess("EMDissociation",particle);
      G4EMDissociation *theEMD = new G4EMDissociation();
      theEMD->SetMinEnergy(100.0*MeV);
      for (iter = elementTable->begin(); iter != elementTable->end(); ++iter)
        theEMD->ActivateFor(*iter);
      theEMD->SetVerboseLevel(2);
      G4EMDissociationCrossSection *EMDCrossSection =
        new G4EMDissociationCrossSection;
      fGenericIon->AddDataSet(EMDCrossSection);
      fGenericIon->RegisterMe(theEMD);
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
