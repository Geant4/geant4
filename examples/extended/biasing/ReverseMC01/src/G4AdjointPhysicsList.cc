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
/// \file G4AdjointPhysicsList.cc
/// \brief Implementation of the G4AdjointPhysicsList class

//  Class Name:        G4AdjointPhysicsList
//        Author:               L. Desorgher
//         Organisation:         SpaceIT GmbH
//        Contract:        ESA contract 21435/08/NL/AT
//         Customer:             ESA/ESTEC
//////////////////////////////////////////////////////////////

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4AdjointPhysicsList.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4AdjointPhysicsMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericIon.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4AdjointGenericIon.hh"
#include "G4AnalysisUtilities.hh"

#include "G4AdjointGamma.hh"
#include "G4AdjointElectron.hh"
#include "G4AdjointProton.hh"
#include "G4AdjointTriton.hh"
#include "G4AdjointDeuteron.hh"
#include "G4AdjointAlpha.hh"
#include "G4AdjointHe3.hh"

#include "G4AdjointAlongStepWeightCorrection.hh"
#include "G4AdjointBremsstrahlungModel.hh"
#include "G4AdjointCSManager.hh"
#include "G4AdjointComptonModel.hh"
#include "G4AdjointForcedInteractionForGamma.hh"
#include "G4AdjointIonIonisationModel.hh"
#include "G4AdjointPhotoElectricModel.hh"
#include "G4AdjointProcessEquivalentToDirectProcess.hh"
#include "G4AdjointSimManager.hh"
#include "G4AdjointeIonisationModel.hh"
#include "G4AdjointhMultipleScattering.hh"
#include "G4AdjointhIonisationModel.hh"
#include "G4ComptonScattering.hh"
#include "G4ContinuousGainOfEnergy.hh"
#include "G4GammaConversion.hh"
#include "G4InversePEEffect.hh"
#include "G4IonInverseIonisation.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UrbanAdjointMscModel.hh"
#include "G4UrbanMscModel.hh"
#include "G4eAdjointMultipleScattering.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eInverseBremsstrahlung.hh"
#include "G4eInverseCompton.hh"
#include "G4eInverseIonisation.hh"
#include "G4eIonisation.hh"
#include "G4eMultipleScattering.hh"
#include "G4eplusAnnihilation.hh"
#include "G4hInverseIonisation.hh"
#include "G4hIonisation.hh"
#include "G4hMultipleScattering.hh"
#include "G4ionIonisation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4AdjointPhysicsList::G4AdjointPhysicsList()
    : G4VUserPhysicsList(),
      fUse_forced_interaction(true),
      fUse_eionisation(true),
      fUse_pionisation(true),
      fUse_brem(true),
      fUse_compton(true),
      fUse_ms(true),
      fUse_egain_fluctuation(true),
      fUse_peeffect(true),
      fEmin_adj_models(1. * keV),
      fEmax_adj_models(1. * MeV),
      fCS_biasing_factor_compton(1.),
      fCS_biasing_factor_brem(1.),
      fCS_biasing_factor_ionisation(1.),
      fCS_biasing_factor_PEeffect(1.)
{
  defaultCutValue = 1.0 * mm;
  SetVerboseLevel(1);
  fPhysicsMessenger = new G4AdjointPhysicsMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4AdjointPhysicsList::~G4AdjointPhysicsList()
{;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4AdjointPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program.
  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();
  G4GenericIon::GenericIonDefinition();
  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();
  ConstructAdjointParticles();
  // Required by MT even if ion physics not used
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4AdjointPhysicsList::ConstructBosons()
{
  G4BosonConstructor pBosonConstructor;
  pBosonConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4AdjointPhysicsList::ConstructLeptons()
{
  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4AdjointPhysicsList::ConstructMesons()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4AdjointPhysicsList::ConstructBaryons()
{
  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4AdjointPhysicsList::ConstructAdjointParticles()
{
  G4AdjointGamma::AdjointGamma();
  G4AdjointElectron::AdjointElectron();
  G4AdjointProton::AdjointProton();
  G4AdjointTriton::Definition();
  G4AdjointDeuteron::Definition();
  G4AdjointAlpha::Definition();
  G4AdjointHe3::Definition();
  G4AdjointGenericIon::Definition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4AdjointPhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4AdjointPhysicsList::ConstructEM()
{
  G4AdjointCSManager *theCSManager =  G4AdjointCSManager::GetAdjointCSManager();
  G4AdjointSimManager *theAdjointSimManager = G4AdjointSimManager::GetInstance();

  theCSManager->RegisterAdjointParticle(G4AdjointElectron::AdjointElectron());

  if (fUse_brem || fUse_peeffect || fUse_compton)
    theCSManager->RegisterAdjointParticle(G4AdjointGamma::AdjointGamma());

  G4eIonisation *fEminusIonisation = nullptr;

  if (fUse_eionisation)
  {
    fEminusIonisation = new G4eIonisation();
    fEminusIonisation->SetLossFluctuations(fUse_egain_fluctuation);
  }
  G4hIonisation *hPIonisation = nullptr;

  if (fUse_pionisation)
  {

    
    hPIonisation = new G4hIonisation();
    hPIonisation->SetLossFluctuations(fUse_egain_fluctuation);
    theCSManager->RegisterAdjointParticle(G4AdjointProton::AdjointProton());
  
  }

  G4eBremsstrahlung *theeminusBremsstrahlung = nullptr;
  if (fUse_brem && fUse_eionisation)
    theeminusBremsstrahlung = new G4eBremsstrahlung();

  G4ComptonScattering *theComptonScattering = nullptr;
  if (fUse_compton)
    theComptonScattering = new G4ComptonScattering();

  G4PhotoElectricEffect *thePEEffect = nullptr;
  if (fUse_peeffect)
    thePEEffect = new G4PhotoElectricEffect();

  G4eMultipleScattering *theeminusMS = nullptr;
  G4hMultipleScattering *thepMS = nullptr;
  G4eAdjointMultipleScattering *theeminusAdjointMS = nullptr;
  if (fUse_ms)
  {
    theeminusMS = new G4eMultipleScattering();
    G4UrbanMscModel *msc1 = new G4UrbanMscModel();
    theeminusMS->SetEmModel(msc1);
    theeminusAdjointMS = new G4eAdjointMultipleScattering();
    G4UrbanAdjointMscModel *msc2 = new G4UrbanAdjointMscModel();
    theeminusAdjointMS->SetEmModel(msc2);
    thepMS = new G4hMultipleScattering();
  }

  G4VProcess *theGammaConversion = 0;
  if (fUse_gamma_conversion)
    theGammaConversion = new G4GammaConversion();
  //Define adjoint e- ionisation
  //-------------------
  G4AdjointeIonisationModel *theeInverseIonisationModel = 0;
  G4eInverseIonisation *theeInverseIonisationProjToProjCase = 0;
  G4eInverseIonisation *theeInverseIonisationProdToProjCase = 0;
  if (fUse_eionisation)
  {
    theeInverseIonisationModel = new G4AdjointeIonisationModel();
    theeInverseIonisationModel->SetHighEnergyLimit(
        fEmax_adj_models);
    theeInverseIonisationModel->SetLowEnergyLimit(
        fEmin_adj_models);
    theeInverseIonisationModel->SetCSBiasingFactor(
        fCS_biasing_factor_ionisation);
    theeInverseIonisationProjToProjCase =
        new G4eInverseIonisation(true, "Inv_eIon",
                                 theeInverseIonisationModel);
    theeInverseIonisationProdToProjCase =
        new G4eInverseIonisation(false, "Inv_eIon1",
                                 theeInverseIonisationModel);
    theAdjointSimManager->ConsiderParticleAsPrimary(G4String("e-"));
  }

  //Define  adjoint Bremsstrahlung
  //-------------------------------
  G4AdjointBremsstrahlungModel *theeInverseBremsstrahlungModel = 0;
  G4eInverseBremsstrahlung *theeInverseBremsstrahlungProjToProjCase = 0;
  G4eInverseBremsstrahlung *theeInverseBremsstrahlungProdToProjCase = 0;
  G4AdjointForcedInteractionForGamma *theForcedInteractionForGamma = 0;

  if (fUse_brem && fUse_eionisation)
  {
    theeInverseBremsstrahlungModel = new G4AdjointBremsstrahlungModel();
    theeInverseBremsstrahlungModel->SetHighEnergyLimit(fEmax_adj_models * 1.01);
    theeInverseBremsstrahlungModel->SetLowEnergyLimit(fEmin_adj_models);
    theeInverseBremsstrahlungModel->SetCSBiasingFactor(
        fCS_biasing_factor_brem);
    theeInverseBremsstrahlungProjToProjCase = new G4eInverseBremsstrahlung(
        true, "Inv_eBrem", theeInverseBremsstrahlungModel);
    if (!fUse_forced_interaction)
      theeInverseBremsstrahlungProdToProjCase =
      new G4eInverseBremsstrahlung(false, G4String("Inv_eBrem1"),
                                   theeInverseBremsstrahlungModel);
    theAdjointSimManager->ConsiderParticleAsPrimary("e-");
    theAdjointSimManager->ConsiderParticleAsPrimary("gamma");

    if (fUse_forced_interaction)
    {
      theForcedInteractionForGamma =
        new G4AdjointForcedInteractionForGamma("ReverseGammaForcedInteraction");
      theForcedInteractionForGamma->RegisterAdjointBremModel(
        theeInverseBremsstrahlungModel);
    }
  }

  //Define  adjoint Compton
  //---------------------
  G4AdjointComptonModel *theeInverseComptonModel = nullptr;
  G4eInverseCompton *theeInverseComptonProjToProjCase = nullptr;
  G4eInverseCompton *theeInverseComptonProdToProjCase = nullptr;

  if (fUse_compton)
  {
    theeInverseComptonModel = new G4AdjointComptonModel();
    theeInverseComptonModel->SetHighEnergyLimit(fEmax_adj_models);
    theeInverseComptonModel->SetLowEnergyLimit(fEmin_adj_models);
    theeInverseComptonModel->SetDirectProcess(theComptonScattering);
    theeInverseComptonModel->SetUseMatrix(false);

    theeInverseComptonModel->SetCSBiasingFactor(fCS_biasing_factor_compton);
    if (!fUse_forced_interaction)
      theeInverseComptonProjToProjCase =
          new G4eInverseCompton(true, "Inv_Compt", theeInverseComptonModel);
    theeInverseComptonProdToProjCase = new G4eInverseCompton(false, "Inv_Compt1",
                                                             theeInverseComptonModel);
    if (fUse_forced_interaction)
    {
      if (!fUse_brem && !fUse_eionisation)
        theForcedInteractionForGamma =
            new G4AdjointForcedInteractionForGamma("ReverseGammaForcedInteraction");
      theForcedInteractionForGamma->RegisterAdjointComptonModel(theeInverseComptonModel);
    }
    theAdjointSimManager->ConsiderParticleAsPrimary(G4String("e-"));
    theAdjointSimManager->ConsiderParticleAsPrimary(G4String("gamma"));
  }

  //Define  adjoint PEEffect
  //---------------------
  G4AdjointPhotoElectricModel *theInversePhotoElectricModel = 0;
  G4InversePEEffect *theInversePhotoElectricProcess = 0;

  if (fUse_peeffect)
  {
    theInversePhotoElectricModel = new G4AdjointPhotoElectricModel();
    theInversePhotoElectricModel->SetHighEnergyLimit(fEmax_adj_models);
    theInversePhotoElectricModel->SetLowEnergyLimit(fEmin_adj_models);
    theInversePhotoElectricModel->SetCSBiasingFactor(
        fCS_biasing_factor_PEeffect);
    theInversePhotoElectricProcess
      = new G4InversePEEffect("Inv_PEEffect", theInversePhotoElectricModel);
    theAdjointSimManager->ConsiderParticleAsPrimary(G4String("e-"));
    theAdjointSimManager->ConsiderParticleAsPrimary(G4String("gamma"));
  }

  //Define  adjoint ionisation for protons
  //---------------------
  G4AdjointhIonisationModel *thepInverseIonisationModel = 0;
  G4hInverseIonisation *thepInverseIonisationProjToProjCase = 0;
  G4hInverseIonisation *thepInverseIonisationProdToProjCase = 0;
  if (fUse_pionisation)
  {
    thepInverseIonisationModel = new G4AdjointhIonisationModel(
        G4Proton::Proton());
    thepInverseIonisationModel->SetHighEnergyLimit(fEmax_adj_models);
    thepInverseIonisationModel->SetLowEnergyLimit(fEmin_adj_models);
    thepInverseIonisationModel->SetUseMatrix(false);
    thepInverseIonisationProjToProjCase
      = new G4hInverseIonisation(true, "Inv_pIon", thepInverseIonisationModel);
    thepInverseIonisationProdToProjCase
      = new G4hInverseIonisation(false, "Inv_pIon1", thepInverseIonisationModel);
    theAdjointSimManager->ConsiderParticleAsPrimary(G4String("e-"));
    theAdjointSimManager->ConsiderParticleAsPrimary(G4String("proton"));
  }

  //Declare the processes active for the different particles
  //--------------------------------------------------------
  auto particleIterator = GetParticleIterator();

  //size_t s1 = theParticleTable->size();
  particleIterator->reset();
  while ((*particleIterator)())
  {
    G4ParticleDefinition *particle = particleIterator->value();
    G4ProcessManager *pmanager = particle->GetProcessManager();
    if (!pmanager)
    {
      pmanager = new G4ProcessManager(particle);
      particle->SetProcessManager(pmanager);
    }

    G4String particleName = particle->GetParticleName();
    if (particleName == "e-")
    {
      if (fUse_ms && fUse_eionisation)
        pmanager->AddProcess(theeminusMS);
      if (fUse_eionisation)
      {
        pmanager->AddProcess(fEminusIonisation);
        G4AdjointCSManager::GetAdjointCSManager()
          ->RegisterEnergyLossProcess(fEminusIonisation, particle);
      }
      if (fUse_brem && fUse_eionisation)
      {
        pmanager->AddProcess(theeminusBremsstrahlung);
        G4AdjointCSManager::GetAdjointCSManager()
          ->RegisterEnergyLossProcess(theeminusBremsstrahlung, particle);
      }
      G4int n_order = 0;
      if (fUse_ms && fUse_eionisation)
      {
        n_order++;
        pmanager->SetProcessOrdering(theeminusMS, idxAlongStep, n_order);
      }
      if (fUse_eionisation)
      {
        n_order++;
        pmanager->SetProcessOrdering(fEminusIonisation, idxAlongStep, n_order);
      }
      if (fUse_brem && fUse_eionisation)
      {
        n_order++;
        pmanager->SetProcessOrdering(theeminusBremsstrahlung,
                                     idxAlongStep, n_order);
      }
      n_order = 0;
      if (fUse_ms && fUse_eionisation)
      {
        n_order++;
        pmanager->SetProcessOrdering(theeminusMS, idxPostStep, n_order);
      }
      if (fUse_eionisation)
      {
        n_order++;
        pmanager->SetProcessOrdering(fEminusIonisation, idxPostStep, n_order);
      }
      if (fUse_brem && fUse_eionisation)
      {
        n_order++;
        pmanager->SetProcessOrdering(theeminusBremsstrahlung, idxPostStep,
                                     n_order);
      }
    }

    if (particleName == "adj_e-")
    {
      G4ContinuousGainOfEnergy *theContinuousGainOfEnergy = 0;
      if (fUse_eionisation)
      {
        theContinuousGainOfEnergy = new G4ContinuousGainOfEnergy();
        theContinuousGainOfEnergy->SetLossFluctuations(
            fUse_egain_fluctuation);
        theContinuousGainOfEnergy->SetDirectEnergyLossProcess(
            fEminusIonisation);
        theContinuousGainOfEnergy->SetDirectParticle(G4Electron::Electron());
        pmanager->AddProcess(theContinuousGainOfEnergy);
      }
      G4int n_order = 0;
      if (fUse_ms)
      {
        n_order++;
        pmanager->AddProcess(theeminusAdjointMS);
        pmanager->SetProcessOrdering(theeminusAdjointMS,
                                     idxAlongStep, n_order);
      }
      n_order++;
      pmanager->SetProcessOrdering(theContinuousGainOfEnergy, idxAlongStep,
                                   n_order);

      n_order++;
      G4AdjointAlongStepWeightCorrection *theAlongStepWeightCorrection =
          new G4AdjointAlongStepWeightCorrection();
      pmanager->AddProcess(theAlongStepWeightCorrection);
      pmanager->SetProcessOrdering(theAlongStepWeightCorrection,
                                   idxAlongStep,
                                   n_order);
      n_order = 0;
      if (fUse_eionisation)
      {
        pmanager->AddProcess(theeInverseIonisationProjToProjCase);
        pmanager->AddProcess(theeInverseIonisationProdToProjCase);
        n_order++;
        pmanager->SetProcessOrdering(theeInverseIonisationProjToProjCase,
                                     idxPostStep, n_order);
        n_order++;
        pmanager->SetProcessOrdering(theeInverseIonisationProdToProjCase,
                                     idxPostStep, n_order);
      }
      if (fUse_brem && fUse_eionisation)
      {
        pmanager->AddProcess(theeInverseBremsstrahlungProjToProjCase);
        n_order++;
        pmanager->SetProcessOrdering(
            theeInverseBremsstrahlungProjToProjCase,
            idxPostStep, n_order);
      }

      if (fUse_compton)
      {
        pmanager->AddProcess(theeInverseComptonProdToProjCase);
        n_order++;
        pmanager->SetProcessOrdering(theeInverseComptonProdToProjCase,
                                     idxPostStep, n_order);
      }
      if (fUse_peeffect)
      {
        pmanager->AddDiscreteProcess(theInversePhotoElectricProcess);
        n_order++;
        pmanager->SetProcessOrdering(theInversePhotoElectricProcess,
                                     idxPostStep, n_order);
      }
      if (fUse_pionisation)
      {
        pmanager->AddProcess(thepInverseIonisationProdToProjCase);
        n_order++;
        pmanager->SetProcessOrdering(thepInverseIonisationProdToProjCase,
                                     idxPostStep, n_order);
      }
      if (fUse_ms && fUse_eionisation)
      {
        n_order++;
        pmanager->SetProcessOrdering(theeminusAdjointMS,
                                     idxPostStep, n_order);
      }
    }

    if (particleName == "adj_gamma")
    {
      G4int n_order = 0;
      if (!fUse_forced_interaction)
      {
        G4AdjointAlongStepWeightCorrection *theAlongStepWeightCorrection =
            new G4AdjointAlongStepWeightCorrection();
        pmanager->AddProcess(theAlongStepWeightCorrection);
        pmanager->SetProcessOrdering(theAlongStepWeightCorrection,
                                     idxAlongStep, 1);

        if (fUse_brem && fUse_eionisation)
        {
          pmanager->AddProcess(theeInverseBremsstrahlungProdToProjCase);
          n_order++;
          pmanager->SetProcessOrdering(
              theeInverseBremsstrahlungProdToProjCase,
              idxPostStep, n_order);
        }
        if (fUse_compton)
        {
          pmanager->AddDiscreteProcess(theeInverseComptonProjToProjCase);
          n_order++;
          pmanager->SetProcessOrdering(theeInverseComptonProjToProjCase,
                                       idxPostStep, n_order);
        }
      }
      else
      {
        //        if (theForcedInteractionForGamma)
        {
          pmanager->AddProcess(theForcedInteractionForGamma);
          n_order++;
          pmanager->SetProcessOrdering(theForcedInteractionForGamma,
                                       idxPostStep, n_order);
          pmanager->SetProcessOrdering(theForcedInteractionForGamma,
                                       idxAlongStep, n_order);
        }
      }
    }

    if (particleName == "gamma")
    {
      if (fUse_compton)
      {
        pmanager->AddDiscreteProcess(theComptonScattering);
        G4AdjointCSManager::GetAdjointCSManager()
          ->RegisterEmProcess(theComptonScattering, particle);
      }
      if (fUse_peeffect)
      {
        pmanager->AddDiscreteProcess(thePEEffect);
        G4AdjointCSManager::GetAdjointCSManager()
          ->RegisterEmProcess(thePEEffect, particle);
      }
      if (fUse_gamma_conversion)
      {
        pmanager->AddDiscreteProcess(theGammaConversion);
      }
    }

    if (particleName == "e+" && fUse_gamma_conversion)
    { //positron
      G4VProcess *theeplusMultipleScattering = new G4eMultipleScattering();
      G4VProcess *theeplusIonisation = new G4eIonisation();
      G4VProcess *theeplusBremsstrahlung = new G4eBremsstrahlung();
      G4VProcess *theeplusAnnihilation = new G4eplusAnnihilation();

      // add processes
      pmanager->AddProcess(theeplusMultipleScattering);
      pmanager->AddProcess(theeplusIonisation);
      pmanager->AddProcess(theeplusBremsstrahlung);
      pmanager->AddProcess(theeplusAnnihilation);

      // set ordering for AtRestDoIt
      pmanager->SetProcessOrderingToFirst(theeplusAnnihilation, idxAtRest);

      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(theeplusMultipleScattering,
                                   idxAlongStep, 1);
      pmanager->SetProcessOrdering(theeplusIonisation, idxAlongStep, 2);
      pmanager->SetProcessOrdering(theeplusBremsstrahlung, idxAlongStep, 3);

      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(theeplusMultipleScattering,
                                   idxPostStep, 1);
      pmanager->SetProcessOrdering(theeplusIonisation, idxPostStep, 2);
      pmanager->SetProcessOrdering(theeplusBremsstrahlung, idxPostStep, 3);
      pmanager->SetProcessOrdering(theeplusAnnihilation, idxPostStep, 4);
    }
    if (particleName == "proton" && fUse_pionisation)
    {
      if (fUse_ms && fUse_pionisation)
        pmanager->AddProcess(thepMS);

      if (fUse_pionisation)
      {
        pmanager->AddProcess(hPIonisation);
        G4AdjointCSManager::GetAdjointCSManager()
          ->RegisterEnergyLossProcess(hPIonisation, particle);
      }

      G4int n_order = 0;
      if (fUse_ms && fUse_pionisation)
      {
        n_order++;
        pmanager->SetProcessOrdering(thepMS, idxAlongStep, n_order);
      }

      if (fUse_pionisation)
      {
        n_order++;
        pmanager->SetProcessOrdering(hPIonisation, idxAlongStep, n_order);
      }

      n_order = 0;
      if (fUse_ms && fUse_pionisation)
      {
        n_order++;
        pmanager->SetProcessOrdering(thepMS, idxPostStep, n_order);
      }

      if (fUse_pionisation)
      {
        n_order++;
        pmanager->SetProcessOrdering(hPIonisation, idxPostStep, n_order);
      }
    }

    if (particleName == "adj_proton" && fUse_pionisation)
    {
      G4ContinuousGainOfEnergy *theContinuousGainOfEnergy = 0;
      if (fUse_pionisation)
      {
        theContinuousGainOfEnergy = new G4ContinuousGainOfEnergy();
        theContinuousGainOfEnergy->SetLossFluctuations(
            fUse_egain_fluctuation);
        theContinuousGainOfEnergy->SetDirectEnergyLossProcess(hPIonisation);
        theContinuousGainOfEnergy->SetDirectParticle(G4Proton::Proton());
        pmanager->AddProcess(theContinuousGainOfEnergy);
      }

      G4int n_order = 0;
      if (fUse_ms)
      {
        n_order++;
        pmanager->AddProcess(thepMS);
        pmanager->SetProcessOrdering(thepMS, idxAlongStep, n_order);
      }

      n_order++;
      pmanager->SetProcessOrdering(theContinuousGainOfEnergy, idxAlongStep,
                                   n_order);

      n_order++;
      G4AdjointAlongStepWeightCorrection *theAlongStepWeightCorrection =
          new G4AdjointAlongStepWeightCorrection();
      pmanager->AddProcess(theAlongStepWeightCorrection);
      pmanager->SetProcessOrdering(theAlongStepWeightCorrection,
                                   idxAlongStep,
                                   n_order);
      n_order = 0;
      if (fUse_pionisation)
      {
        pmanager->AddProcess(thepInverseIonisationProjToProjCase);
        n_order++;
        pmanager->SetProcessOrdering(
            thepInverseIonisationProjToProjCase,
            idxPostStep, n_order);
      }

      if (fUse_ms && fUse_pionisation)
      {
        n_order++;
        pmanager->SetProcessOrdering(thepMS, idxPostStep, n_order);
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Decay.hh"
void G4AdjointPhysicsList::ConstructGeneral()
{

  // Add Decay Process
  G4Decay *theDecayProcess = new G4Decay();

  auto particleIterator = GetParticleIterator();
  particleIterator->reset();
  while ((*particleIterator)())
  {
    G4ParticleDefinition *particle = particleIterator->value();
    G4ProcessManager *pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle))
    {
      pmanager->AddProcess(theDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4AdjointPhysicsList::SetCuts()
{
  if (verboseLevel > 0)
  {
    G4cout << "G4AdjointPhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue, "Length") << G4endl;
  }

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  //
  SetCutValue(defaultCutValue, "gamma");
  SetCutValue(defaultCutValue, "e-");
  SetCutValue(defaultCutValue, "e+");

  if (verboseLevel > 0)
    DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
