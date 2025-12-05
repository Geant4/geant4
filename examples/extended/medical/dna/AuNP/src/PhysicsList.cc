//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file PhysicsList.cc
/// \brief Implementation of the PhysicsList class

#include "PhysicsList.hh"

#include "DetectorConstruction.hh"
#include "PhysicsListMessenger.hh"

#include "G4PhysicsConstructorRegistry.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

#include "G4ComptonScattering.hh"
#include "G4DNAAttachment.hh"
#include "G4DNABornExcitationModel.hh"
#include "G4DNABornIonisationModel.hh"
#include "G4DNAChampionElasticModel.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNADingfelderChargeDecreaseModel.hh"
#include "G4DNADiracRMatrixExcitationModel.hh"
#include "G4DNAELSEPAElasticModel.hh"
#include "G4DNAElastic.hh"
#include "G4DNAExcitation.hh"
#include "G4DNAGenericIonsManager.hh"
#include "G4DNAIonisation.hh"
#include "G4DNAMeltonAttachmentModel.hh"
#include "G4DNAPlasmonExcitation.hh"
#include "G4DNAQuinnPlasmonExcitationModel.hh"
#include "G4DNARelativisticIonisationModel.hh"
#include "G4DNASancheExcitationModel.hh"
#include "G4DNAVibExcitation.hh"
#include "G4DummyModel.hh"
#include "G4EmConfigurator.hh"
#include "G4Gamma.hh"
#include "G4GammaConversion.hh"
#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4IonFluctuations.hh"
#include "G4LivermoreBremsstrahlungModel.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4LossTableManager.hh"
#include "G4PenelopeBremsstrahlungModel.hh"
#include "G4PenelopeIonisationModel.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"
#include "G4SeltzerBergerModel.hh"
#include "G4StepLimiter.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4UniversalFluctuation.hh"
#include "G4UrbanMscModel.hh"
#include "G4VEmModel.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eIonisation.hh"
#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4ProcessTable.hh"
#include "G4EmDNAChemistry_option3.hh"
#include "G4GenericIon.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PhysicsList::PhysicsList() {
  fpDetector = dynamic_cast<const DetectorConstruction *>(
    G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  defaultCutValue = 0.1 * CLHEP::nanometer;
  fPhysMessenger = std::make_unique<PhysicsListMessenger>(this);
  SetVerboseLevel(1);
  fEmDNAChemistryList = std::make_unique<G4EmDNAChemistry_option3>();
  G4EmParameters *param = G4EmParameters::Instance();
  param->SetMinEnergy(100 * eV);
  param->SetMaxEnergy(1 * GeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructParticle() {
  ConstructBosons();
  ConstructLeptons();
  ConstructBarions();
  if (fEmDNAChemistryList != nullptr) {
    fEmDNAChemistryList->ConstructParticle();
  }
  G4VModularPhysicsList::ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructBosons() {
  G4Gamma::GammaDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructLeptons() {
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructBarions() {
  G4Proton::ProtonDefinition();
  G4GenericIon::GenericIonDefinition();

  auto *genericIonsManager = G4DNAGenericIonsManager::Instance();
  genericIonsManager->GetIon("alpha++");
  genericIonsManager->GetIon("alpha+");
  genericIonsManager->GetIon("helium");
  genericIonsManager->GetIon("hydrogen");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructProcess() {
  ConstructEM();
  ConstructGeneral();
  if (fEmDNAChemistryList != nullptr) {
    fEmDNAChemistryList->ConstructProcess();
  }
  // Contruct processes of the chemistry list
  G4VModularPhysicsList::ConstructProcess();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::ConstructEM() {
  G4bool isGNP = false;
  if (fphysname == "") fphysname = "DNA";

  const G4Material *NPMaterial = fpDetector->GetNPMaterial();
  const G4String npname = NPMaterial->GetName();

  if (npname == "G4_Au") isGNP = true;

  const auto theParticleIterator = GetParticleIterator();
  theParticleIterator->reset();

  while ((*theParticleIterator)()) {
    const G4ParticleDefinition *particle = theParticleIterator->value();

    G4ProcessManager *pm = particle->GetProcessManager();

    G4String particleName = particle->GetParticleName();

    if (particleName == "e-") {
      const auto theDNAElasticProcess = new G4DNAElastic("e-_G4DNAElastic");
      theDNAElasticProcess->SetEmModel(new G4DNAChampionElasticModel(), 1);
      pm->AddDiscreteProcess(theDNAElasticProcess);

      const auto theDNAExcitationProcess = new G4DNAExcitation("e-_G4DNAExcitation");
      theDNAExcitationProcess->SetEmModel(new G4DNABornExcitationModel(), 1);
      pm->AddDiscreteProcess(theDNAExcitationProcess);

      const auto theDNAIonizationProcess = new G4DNAIonisation("e-_G4DNAIonisation");
      theDNAIonizationProcess->SetEmModel(new G4DNABornIonisationModel(), 1);
      pm->AddDiscreteProcess(theDNAIonizationProcess);

      const auto theDNAAttachmentProcess = new G4DNAAttachment("e-_G4DNAAttachment");
      pm->AddDiscreteProcess(theDNAAttachmentProcess);

      const auto theDNAVibExcProcess = new G4DNAVibExcitation("e-_G4DNAVibExcitation");
      pm->AddDiscreteProcess(theDNAVibExcProcess);

      const auto theDNABremProcess = new G4eBremsstrahlung("e-_G4DNABremsstrahlung");
      pm->AddDiscreteProcess(theDNABremProcess);

      if (isGNP) {
        if (fphysname == "DNA") {
          const auto theDNAELSEPAElasticProcess = new G4DNAElastic("e-_G4DNAELSEPAElastic");
          theDNAELSEPAElasticProcess->SetEmModel(new G4DummyModel(), 1);
          pm->AddDiscreteProcess(theDNAELSEPAElasticProcess);
          const auto theDNADRMExcitationProcess =
              new G4DNAExcitation("e-_G4DNADRMExcitation");
          theDNADRMExcitationProcess->SetEmModel(new G4DummyModel(), 1);
          pm->AddDiscreteProcess(theDNADRMExcitationProcess);
          const auto theDNARelativisticIonizationProcess =
              new G4DNAIonisation("e-_G4DNARelativisticIonisation");
          theDNARelativisticIonizationProcess->SetEmModel(new G4DummyModel(), 1);
          pm->AddDiscreteProcess(theDNARelativisticIonizationProcess);
          const auto theDNAPExcitationProcess =
              new G4DNAPlasmonExcitation("e-_G4DNAPlasmonExcitation");
          theDNAPExcitationProcess->SetEmModel(new G4DummyModel(), 1);
          pm->AddDiscreteProcess(theDNAPExcitationProcess);
          const auto theDNAeBremProcess =
              new G4eBremsstrahlung("e-_G4DNABremsstrahlung_GNP");
          theDNAeBremProcess->SetEmModel(new G4DummyModel(), 1);
          pm->AddDiscreteProcess(theDNAeBremProcess);
        }
        if (fphysname == "Livermore") {
          const auto msc = new G4eMultipleScattering();
          msc->SetEmModel(new G4DummyModel(), 1);
          pm->AddDiscreteProcess(msc);
          const auto eIoni = new G4eIonisation("e-_G4LivermoreIoni");
          eIoni->SetEmModel(new G4DummyModel(), 1);
          pm->AddDiscreteProcess(eIoni);
          const auto eBrem = new G4eBremsstrahlung("e-_G4LivermoreBrem");
          eBrem->SetEmModel(new G4DummyModel(), 1);
          pm->AddDiscreteProcess(eBrem);

          const auto steplimit = new G4StepLimiter();
          pm->AddDiscreteProcess(steplimit);
        }
        if (fphysname == "Penelope") {
          const auto msc = new G4eMultipleScattering();
          msc->SetEmModel(new G4DummyModel(), 1);
          pm->AddDiscreteProcess(msc);
          const auto eIoni = new G4eIonisation("e-_G4PenelopeIoni");
          eIoni->SetEmModel(new G4DummyModel(), 1);
          pm->AddDiscreteProcess(eIoni);
          const auto eBrem = new G4eBremsstrahlung("e-_G4PenelopeBrem");
          eBrem->SetEmModel(new G4DummyModel(), 1);
          pm->AddDiscreteProcess(eBrem);
          const auto steplimit = new G4StepLimiter();
          pm->AddDiscreteProcess(steplimit);
        }
      }

      //// Capture of low-energy e-
      // G4ElectronCapture* ecap = new G4ElectronCapture("NP",10.*eV);
      // pm->AddDiscreteProcess(ecap);
    } else if (particleName == "proton") {
    } else if (particleName == "hydrogen") {
    } else if (particleName == "gamma") {
      G4PhotoElectricEffect *thePhotoElectricEffect = new G4PhotoElectricEffect();
      thePhotoElectricEffect->SetEmModel(new G4LivermorePhotoElectricModel(), 1);
      pm->AddDiscreteProcess(thePhotoElectricEffect);

      G4ComptonScattering *theComptonScattering = new G4ComptonScattering();
      theComptonScattering->SetEmModel(new G4LivermoreComptonModel(), 1);
      pm->AddDiscreteProcess(theComptonScattering);

      G4GammaConversion *theGammaConversion = new G4GammaConversion();
      theGammaConversion->SetEmModel(new G4LivermoreGammaConversionModel(), 1);
      pm->AddDiscreteProcess(theGammaConversion);

      G4RayleighScattering *theRayleigh = new G4RayleighScattering();
      pm->AddDiscreteProcess(theRayleigh);
    }
  }

  if (isGNP) {
    G4EmConfigurator *em_config = G4LossTableManager::Instance()->EmConfigurator();
    G4VEmModel *mod;
    mod = new G4DNAChampionElasticModel();
    mod->SetActivationLowEnergyLimit(1 * GeV);
    em_config->SetExtraEmModel("e-", "e-_G4DNAElastic", mod, "NP");
    mod = new G4DNABornExcitationModel();
    mod->SetActivationLowEnergyLimit(1 * GeV);
    em_config->SetExtraEmModel("e-", "e-_G4DNAExcitation", mod, "NP");
    mod = new G4DNABornIonisationModel();
    mod->SetActivationLowEnergyLimit(1 * GeV);
    em_config->SetExtraEmModel("e-", "e-_G4DNAIonisation", mod, "NP");
    mod = new G4DNAMeltonAttachmentModel();
    mod->SetActivationLowEnergyLimit(1 * GeV);
    em_config->SetExtraEmModel("e-", "e-_G4DNAAttachment", mod, "NP");
    mod = new G4DNASancheExcitationModel();
    mod->SetActivationLowEnergyLimit(1 * GeV);
    em_config->SetExtraEmModel("e-", "e-_G4DNAVibExcitation", mod, "NP");
    mod = new G4SeltzerBergerModel();
    mod->SetActivationLowEnergyLimit(1 * GeV);
    em_config->SetExtraEmModel("e-", "e-_G4DNABremsstrahlung", mod, "NP");

    if (fphysname == "DNA") {
      mod = new G4DNAELSEPAElasticModel();
      em_config->SetExtraEmModel("e-", "e-_G4DNAELSEPAElastic", mod, "NP", 10 * eV, 1 * GeV);
      mod = new G4DNADiracRMatrixExcitationModel();
      em_config->SetExtraEmModel("e-", "e-_G4DNADRMExcitation", mod, "NP", 10 * eV, 1 * GeV);
      mod = new G4DNARelativisticIonisationModel();
      em_config->SetExtraEmModel("e-", "e-_G4DNARelativisticIonisation", mod, "NP", 10 * eV,
                                 1 * GeV);
      mod = new G4DNAQuinnPlasmonExcitationModel();
      em_config->SetExtraEmModel("e-", "e-_G4DNAPlasmonExcitation", mod, "NP", 10 * eV, 1 * GeV);
      mod = new G4SeltzerBergerModel();
      em_config->SetExtraEmModel("e-", "e-_G4DNABremsstrahlung_GNP", mod, "NP", 10 * eV, 1 * GeV);
    }
    if (fphysname == "Livermore") {
      mod = new G4UrbanMscModel();
      // mod = new G4GoudsmitSaundersonMscModel();
      em_config->SetExtraEmModel("e-", "msc", mod, "NP", 0 * eV, 100 * MeV,
                                 new G4UniversalFluctuation());
      mod = new G4LivermoreIonisationModel();
      em_config->SetExtraEmModel("e-", "e-_G4LivermoreIoni", mod, "NP", 0 * eV, 1.0 * MeV,
                                 new G4UniversalFluctuation());
      mod = new G4LivermoreBremsstrahlungModel();
      em_config->SetExtraEmModel("e-", "e-_G4LivermoreBrem", mod, "NP", 0 * eV, 1 * GeV,
                                 new G4UniversalFluctuation());
    }
    if (fphysname == "Penelope") {
      mod = new G4UrbanMscModel();
      // mod = new G4GoudsmitSaundersonMscModel();
      em_config->SetExtraEmModel("e-", "msc", mod, "NP", 0 * eV, 100 * MeV,
                                 new G4UniversalFluctuation());
      mod = new G4PenelopeIonisationModel();
      em_config->SetExtraEmModel("e-", "e-_G4PenelopeIoni", mod, "NP", 0 * eV, 1.0 * MeV,
                                 new G4UniversalFluctuation());
      mod = new G4PenelopeBremsstrahlungModel();
      em_config->SetExtraEmModel("e-", "e-_G4PenelopeBrem", mod, "NP", 0 * eV, 1 * GeV,
                                 new G4UniversalFluctuation());
    }
  }

  G4VAtomDeexcitation *de = new G4UAtomicDeexcitation();
  G4LossTableManager::Instance()->SetAtomDeexcitation(de);
  de->SetFluo(true);
  de->SetPIXE(true);
  de->SetAuger(true);
  de->SetAugerCascade(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysicsList::SetCuts() {
  if (verboseLevel > 0) {
    G4cout << "PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue, "Length") << G4endl;
  }

  SetCutValue(fcutForGamma, "gamma");
  SetCutValue(fcutForElectron, "e-");
  SetCutValue(fcutForPositron, "e+");

  if (verboseLevel > 0) {
    DumpCutValuesTable();
  }
}

void PhysicsList::SetPhysics4NP(const G4String &name) {
  fphysname = name;
}