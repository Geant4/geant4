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

#include "eRositaPhysicsList.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4PhysicsListHelper.hh"
#include "G4ProductionCutsTable.hh"
#include "G4StepLimiter.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicsConstructor.hh"

// Physics List
#include "G4DecayPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmStandardPhysics.hh"

// Process
#include "G4ComptonScattering.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eIonisation.hh"
#include "G4eMultipleScattering.hh"
#include "G4eplusAnnihilation.hh"
#include "G4GammaConversion.hh"
#include "G4hImpactIonisation.hh"
#include "G4hIonisation.hh"
#include "G4hMultipleScattering.hh"
#include "G4ionIonisation.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"

// Model
#include "G4LivermoreBremsstrahlungModel.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4LivermoreRayleighModel.hh"
#include "G4PenelopeAnnihilationModel.hh"
#include "G4UniversalFluctuation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

eRositaPhysicsList::eRositaPhysicsList()
{
    SetVerboseLevel(1);

    constexpr auto DEFAULT_CUT_VALUE{0.001 * mm};
    SetDefaultCutValue(DEFAULT_CUT_VALUE);

    G4cout << "============================================================================================="
              << G4endl
              << "Geant4 eRosita example, based on a simplified version of eROSITA simulation."
              << G4endl
              << "Further details can be found in:"
              << G4endl
              << " M. G. Pia et al.,"
              << G4endl
              << "  'PIXE Simulation With Geant4',"
              << G4endl
              << "  IEEE Trans. Nucl. Sci., vol. 56, no. 6, pp. 3614-3649, 2009"
              << G4endl
              << " N. Meidinger et al.,"
              << G4endl
              << "  'Development of the focal plane PNCCD camera system for the X-ray space telescope eROSITA',"
              << G4endl
              << "  NIM A 624, 321-329, 2010"
              << G4endl
              << "============================================================================================="
              << G4endl;

    G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

eRositaPhysicsList::~eRositaPhysicsList()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void eRositaPhysicsList::ConstructBosons()
{
    // geantino (pseudo-particle)
    // G4Geantino::GeantinoDefinition();

    // charged geantino (pseudo-particle)
    // G4ChargedGeantino::ChargedGeantinoDefinition();

    // photon (gamma)
    G4Gamma::GammaDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void eRositaPhysicsList::ConstructLeptons()
{
    // leptons

    // e+ / e-
    G4Electron::ElectronDefinition();
    G4Positron::PositronDefinition();

    // mu+ / mu-
    // G4MuonPlus::MuonPlusDefinition();
    // G4MuonMinus::MuonMinusDefinition();

    // nu_e
    // G4NeutrinoE::NeutrinoEDefinition();
    // G4AntiNeutrinoE::AntiNeutrinoEDefinition();

    // nu_mu
    // G4NeutrinoMu::NeutrinoMuDefinition();
    // G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void eRositaPhysicsList::ConstructMesons()
{
    // light mesons

    // pion
    // G4PionPlus::PionPlusDefinition();
    // G4PionMinus::PionMinusDefinition();
    // G4PionZero::PionZeroDefinition();

    // eta
    // G4Eta::EtaDefinition();
    // G4EtaPrime::EtaPrimeDefinition();

    // kaon
    // G4KaonPlus::KaonPlusDefinition();
    // G4KaonMinus::KaonMinusDefinition();
    // G4KaonZero::KaonZeroDefinition();
    // G4AntiKaonZero::AntiKaonZeroDefinition();
    // G4KaonZeroLong::KaonZeroLongDefinition();
    // G4KaonZeroShort::KaonZeroShortDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void eRositaPhysicsList::ConstructBaryons()
{
    // baryons

    // proton
    G4Proton::ProtonDefinition();
    G4AntiProton::AntiProtonDefinition();

    // neutron
    // G4Neutron::NeutronDefinition();
    // G4AntiNeutron::AntiNeutronDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void eRositaPhysicsList::ConstructParticle()
{
    ConstructBosons();
    ConstructLeptons();
    ConstructMesons();
    ConstructBaryons();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void eRositaPhysicsList::ConstructEM()
{
    auto *helper = G4PhysicsListHelper::GetPhysicsListHelper();

    auto particleIterator = GetParticleIterator();
    particleIterator->reset();

    while ((*particleIterator)()) {
        auto *particle = particleIterator->value();
        auto particleName = particle->GetParticleName();

        if (particleName == "gamma") { // photon
            // photoelectric effect
            auto *photoelectricEffect = new G4PhotoElectricEffect();
            // photoelectricEffect->ActivateAuger(true);
            // photoelectricEffect->SetCutForLowEnSecElectrons(0.250 * keV);
            // photoelectricEffect->SetCutForLowEnSecPhotons(0.250 * keV);
            photoelectricEffect->SetEmModel(new G4LivermorePhotoElectricModel()); // Set the model (Livermore)
            helper->RegisterProcess(photoelectricEffect, particle);

            // Compton scattering
            auto *comptonScattering = new G4ComptonScattering();
            comptonScattering->SetEmModel(new G4LivermoreComptonModel()); // Set the model (Livermore)
            helper->RegisterProcess(comptonScattering, particle);

            // gamma conversion
            auto *gammaConversion = new G4GammaConversion();
            gammaConversion->SetEmModel(new G4LivermoreGammaConversionModel()); // Set the model (Livermore)
            helper->RegisterProcess(gammaConversion, particle);

            // Rayleigh scattering
            auto *rayleighScattering = new G4RayleighScattering();
            rayleighScattering->SetEmModel(new G4LivermoreRayleighModel()); // Set the model (Livermore)
            helper->RegisterProcess(rayleighScattering, particle);
        } else if (particleName == "e-") { // electron
            // multiple scattering
            helper->RegisterProcess(new G4eMultipleScattering(), particle);

            // ionization
            auto *ionization = new G4eIonisation();
            ionization->SetEmModel(new G4LivermoreIonisationModel()); // Set the model (Livermore)
            ionization->SetFluctModel(new G4UniversalFluctuation());
            helper->RegisterProcess(ionization, particle);

            // Bremsstrahlung
            auto *bremsstrahlung = new G4eBremsstrahlung();
            bremsstrahlung->SetEmModel(new G4LivermoreBremsstrahlungModel()); // Set the model (Livermore)
            helper->RegisterProcess(bremsstrahlung, particle);
        } else if (particleName == "e+") { // positron
            // multiple scattering
            helper->RegisterProcess(new G4eMultipleScattering(), particle);

            // ionization
            helper->RegisterProcess(new G4eIonisation(), particle);

            // Bremsstrahlung
            auto *bremsstrahlung = new G4eBremsstrahlung();
            bremsstrahlung->SetEmModel(new G4LivermoreBremsstrahlungModel()); // Set the model (Livermore)
            helper->RegisterProcess(bremsstrahlung, particle);

            // annihilation
            auto *annihilation = new G4eplusAnnihilation();
            annihilation->SetEmModel(new G4PenelopeAnnihilationModel()); // Set the model (Penelope)
            helper->RegisterProcess(annihilation, particle);
        // } else if( particleName == "mu+" || particleName == "mu-") {
            // // muon
            // helper->RegisterProcess(new G4MuMultipleScattering, particle);
            // helper->RegisterProcess(new G4MuIonisation, particle);
            // helper->RegisterProcess(new G4MuBremsstrahlung, particle);
            // helper->RegisterProcess(new G4MuPairProduction, particle);
        } else if (particleName == "proton" || particleName == "pi-" || particleName == "pi+") {
            helper->RegisterProcess(new G4hMultipleScattering(), particle);
            helper->RegisterProcess(new G4hIonisation(), particle);
/*
            // proton
            // auto *ionization = new G4hImpactIonisation();
            // ionization->SetPixeCrossSectionK("ecpssr");
            // ionization->SetPixeCrossSectionL("ecpssr");
            // ionization->SetPixeCrossSectionM("ecpssr");
            // ionization->SetPixeProjectileMinEnergy(1.* keV);
            // ionization->SetPixeProjectileMaxEnergy(200. * MeV);
            // ionization->SetCutForSecondaryPhotons(250. * eV);
            // ionization->SetCutForAugerElectrons(250. * eV);

            auto *ionization = new G4hIonisation();
            auto *multipleScattering = new G4hMultipleScattering();

            processManager->AddProcess(multipleScattering, -1, 1, 1);
            processManager->AddProcess(ionization, -1, 2, 2);
*/
        } else if (particleName == "alpha" || particleName == "He3" || particleName == "pi-" || particleName == "pi+" || particleName == "GenericIon") {
            // pion, alpha, ion (should never occur in the current example)
            helper->RegisterProcess(new G4hMultipleScattering, particle);
            helper->RegisterProcess(new G4ionIonisation, particle);
        } else if ((!particle->IsShortLived()) && (particle->GetPDGCharge() != 0.0) && (particle->GetParticleName() != "chargedgeantino")) {
            // every other charged particle, except geantino
            helper->RegisterProcess(new G4hMultipleScattering, particle);
            helper->RegisterProcess(new G4hIonisation, particle);
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void eRositaPhysicsList::ConstructGeneral()
{
    auto *helper = G4PhysicsListHelper::GetPhysicsListHelper();

    // Add decay process
    auto *decay = new G4Decay();

    auto particleIterator = GetParticleIterator();
    particleIterator->reset();

    while ((*particleIterator)()) {
        auto *particle = particleIterator->value();

        if (decay->IsApplicable(*particle)) {
            if (verboseLevel > 1) {
                G4cout << "### Decays for " << particle->GetParticleName() << G4endl;
            }
            helper->RegisterProcess(decay, particle);
/*
            // Set ordering for PostStepDoIt and AtRestDoIt
            processManager->SetProcessOrdering(decay, idxPostStep);
            processManager->SetProcessOrdering(decay, idxAtRest);
*/
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void eRositaPhysicsList::ConstructProcess()
{
    AddTransportation();
    ConstructEM();
    ConstructGeneral();
    // AddStepMax();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void eRositaPhysicsList::SetCuts()
{
    // Set the default cut value for all particle types
    SetCutsWithDefault();

    // Set the secondary production cut lower than 990. eV. Very important for processes at low energies.        
    constexpr auto ENERGY_LOW_LIMIT{250. * eV};
    constexpr auto ENERGY_HIGH_LIMIT{100. * GeV};
    
    G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(ENERGY_LOW_LIMIT, ENERGY_HIGH_LIMIT);

    if (verboseLevel > 0) {
        DumpCutValuesTable();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

/*
#include "G4StepLimiter.hh"
#include "G4UserSpecialCuts.hh"

void eRositaPhysicsList::AddStepMax()
{
    auto *helper = G4PhysicsListHelper::GetPhysicsListHelper();

    // Step limitation seen as a process
    // auto *stepLimiter = new G4StepLimiter();
    // // auto *userCuts = new G4UserSpecialCuts();

    particleIterator->reset();

    while ((*particleIterator)()){
        auto *particle = particleIterator->value();
        // auto *processManager = particle->GetProcessManager();

        if (particle->GetPDGCharge() != 0.0) {
            helper->RegisterProcess(stepLimiter, particle);
            // helper->RegisterProcess(userCuts, particle);
        }
    }
}
*/
