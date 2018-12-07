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
/// \file biasing/ReverseMC01/src/G4AdjointPhysicsList.cc
/// \brief Implementation of the G4AdjointPhysicsList class
//
//
//////////////////////////////////////////////////////////////
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4AdjointPhysicsList::G4AdjointPhysicsList()
 :G4VUserPhysicsList(),
  fEminusIonisation(0),fPIonisation(0),
  fUse_forced_interaction(true),
  fUse_eionisation(true),fUse_pionisation(true),
  fUse_brem(true),fUse_compton(true),fUse_ms(true),
  fUse_egain_fluctuation(true),fUse_peeffect(true),
  fEmin_adj_models(1.*keV), fEmax_adj_models(1.*MeV),
  fCS_biasing_factor_compton(1.),fCS_biasing_factor_brem(1.),
  fCS_biasing_factor_ionisation(1.),fCS_biasing_factor_PEeffect(1.)
{
  defaultCutValue = 1.0*mm;
  SetVerboseLevel(1);
  fPhysicsMessenger = new G4AdjointPhysicsMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4AdjointPhysicsList::~G4AdjointPhysicsList()
{
}
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
  ConstructAdjointParticles();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4AdjointPhysicsList::SetLossFluctuationFlag(bool aBool)
{
 if (fEminusIonisation) fEminusIonisation->SetLossFluctuations(aBool);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4AdjointPhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4AdjointPhysicsList::ConstructLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4AdjointPhysicsList::ConstructMesons()
{
//  mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4Eta::EtaDefinition();
  G4EtaPrime::EtaPrimeDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();
  G4KaonZeroShort::KaonZeroShortDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4AdjointPhysicsList::ConstructBaryons()
{
//  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include"G4AdjointGamma.hh"
#include"G4AdjointElectron.hh"
#include"G4AdjointProton.hh"
void G4AdjointPhysicsList::ConstructAdjointParticles()
{
// adjoint_gammma
  G4AdjointGamma::AdjointGammaDefinition();

// adjoint_electron
  G4AdjointElectron::AdjointElectronDefinition();
  
// adjoint_proton
  G4AdjointProton::AdjointProtonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4AdjointPhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//#include "G4PEEffectFluoModel.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4eMultipleScattering.hh"
#include "G4eAdjointMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
//#include "G4IonParametrisedLossModel.hh"

#include "G4eBremsstrahlung.hh"
#include "G4ContinuousGainOfEnergy.hh"
#include "G4eInverseIonisation.hh"
#include "G4AdjointeIonisationModel.hh"
#include "G4AdjointCSManager.hh"
#include "G4AdjointBremsstrahlungModel.hh"
#include "G4eInverseBremsstrahlung.hh"
#include "G4AdjointComptonModel.hh"
#include "G4eInverseCompton.hh"
#include "G4InversePEEffect.hh"
#include "G4AdjointPhotoElectricModel.hh"
#include "G4AdjointAlongStepWeightCorrection.hh"
#include "G4hInverseIonisation.hh"
#include "G4AdjointhIonisationModel.hh"
#include "G4AdjointhMultipleScattering.hh"
#include "G4IonInverseIonisation.hh"
#include "G4AdjointIonIonisationModel.hh"

#include "G4AdjointSimManager.hh"
#include "G4AdjointProcessEquivalentToDirectProcess.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UrbanMscModel.hh"
#include "G4UrbanAdjointMscModel.hh"
#include "G4UrbanAdjointMscModel.hh"
#include "G4AdjointForcedInteractionForGamma.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4AdjointPhysicsList::ConstructEM()
{
  G4AdjointCSManager* theCSManager =
       G4AdjointCSManager::GetAdjointCSManager();
  G4AdjointSimManager* theAdjointSimManager =
              G4AdjointSimManager::GetInstance();
  
  theCSManager->RegisterAdjointParticle(
            G4AdjointElectron::AdjointElectron());
  
  if (fUse_brem || fUse_peeffect ||fUse_compton)
              theCSManager->RegisterAdjointParticle(
                  G4AdjointGamma::AdjointGamma());

  if (fUse_eionisation) {
          if (!fEminusIonisation) fEminusIonisation  =
                         new G4eIonisation();
        fEminusIonisation->SetLossFluctuations(fUse_egain_fluctuation);
  }
  if (fUse_pionisation) {
          if (!fPIonisation) fPIonisation  = new G4hIonisation();
        fPIonisation->SetLossFluctuations(fUse_egain_fluctuation);
        theCSManager->RegisterAdjointParticle(
                G4AdjointProton::AdjointProton());
  }        
  
  G4eBremsstrahlung* theeminusBremsstrahlung = 0;
  if (fUse_brem && fUse_eionisation)
                  theeminusBremsstrahlung = new G4eBremsstrahlung();
  
  G4ComptonScattering* theComptonScattering =0;
  if (fUse_compton) theComptonScattering = new G4ComptonScattering();
  
  G4PhotoElectricEffect* thePEEffect =0;
  if (fUse_peeffect) thePEEffect = new G4PhotoElectricEffect();
  
  G4eMultipleScattering* theeminusMS = 0;
  G4hMultipleScattering* thepMS= 0;
  G4eAdjointMultipleScattering* theeminusAdjointMS = 0;
  if (fUse_ms) {
         theeminusMS = new G4eMultipleScattering();
         G4UrbanMscModel* msc1 = new G4UrbanMscModel();
         theeminusMS->SetEmModel(msc1);
         theeminusAdjointMS = new G4eAdjointMultipleScattering();
         G4UrbanAdjointMscModel* msc2 = new G4UrbanAdjointMscModel();
         theeminusAdjointMS->SetEmModel(msc2);
         thepMS = new G4hMultipleScattering();
  }        
  
  G4VProcess*  theGammaConversion =0;
  if (fUse_gamma_conversion) theGammaConversion = new G4GammaConversion();
  //Define adjoint e- ionisation
  //-------------------
  G4AdjointeIonisationModel* theeInverseIonisationModel = 0;
  G4eInverseIonisation* theeInverseIonisationProjToProjCase = 0 ;
  G4eInverseIonisation* theeInverseIonisationProdToProjCase = 0;
  if (fUse_eionisation) {
    theeInverseIonisationModel = new G4AdjointeIonisationModel();
    theeInverseIonisationModel->SetHighEnergyLimit(
                                      fEmax_adj_models);
    theeInverseIonisationModel->SetLowEnergyLimit(
                                       fEmin_adj_models);
    theeInverseIonisationModel->SetCSBiasingFactor(
                                    fCS_biasing_factor_ionisation);
    theeInverseIonisationProjToProjCase =
                   new G4eInverseIonisation(true,"Inv_eIon",
                      theeInverseIonisationModel);
    theeInverseIonisationProdToProjCase =
              new G4eInverseIonisation(false,"Inv_eIon1",
                          theeInverseIonisationModel);
    theAdjointSimManager->ConsiderParticleAsPrimary(G4String("e-"));
  }        

  //Define  adjoint Bremsstrahlung
  //-------------------------------
  G4AdjointBremsstrahlungModel* theeInverseBremsstrahlungModel = 0;
  G4eInverseBremsstrahlung* theeInverseBremsstrahlungProjToProjCase = 0;
  G4eInverseBremsstrahlung* theeInverseBremsstrahlungProdToProjCase = 0; 
  G4AdjointForcedInteractionForGamma* theForcedInteractionForGamma = 0;
  if (fUse_brem && fUse_eionisation) {
    theeInverseBremsstrahlungModel = new G4AdjointBremsstrahlungModel();
    theeInverseBremsstrahlungModel->SetHighEnergyLimit(fEmax_adj_models*1.01);
    theeInverseBremsstrahlungModel->SetLowEnergyLimit(fEmin_adj_models);
    theeInverseBremsstrahlungModel->SetCSBiasingFactor(
                                     fCS_biasing_factor_brem);
    theeInverseBremsstrahlungProjToProjCase = new G4eInverseBremsstrahlung(
                             true,"Inv_eBrem",theeInverseBremsstrahlungModel);
    theeInverseBremsstrahlungProdToProjCase = new G4eInverseBremsstrahlung(
        false,"Inv_eBrem1",theeInverseBremsstrahlungModel);
    theAdjointSimManager->ConsiderParticleAsPrimary(G4String("e-"));
    theAdjointSimManager->ConsiderParticleAsPrimary(G4String("gamma"));

    if (!fUse_forced_interaction) theeInverseBremsstrahlungProdToProjCase
       = new G4eInverseBremsstrahlung(false,G4String("Inv_eBrem1"),
                                        theeInverseBremsstrahlungModel);
    theAdjointSimManager->ConsiderParticleAsPrimary(G4String("e-"));
    theAdjointSimManager->ConsiderParticleAsPrimary(G4String("gamma"));
    if (fUse_forced_interaction){
       theForcedInteractionForGamma =
        new G4AdjointForcedInteractionForGamma("ReverseGammaForcedInteraction");
       theForcedInteractionForGamma->RegisterAdjointBremModel(
                                        theeInverseBremsstrahlungModel);
    }
  }

  
  //Define  adjoint Compton
  //---------------------
  
  G4AdjointComptonModel* theeInverseComptonModel = 0;
  G4eInverseCompton* theeInverseComptonProjToProjCase = 0;
  G4eInverseCompton* theeInverseComptonProdToProjCase = 0;

  if (fUse_compton) { 
    theeInverseComptonModel = new G4AdjointComptonModel();
    theeInverseComptonModel->SetHighEnergyLimit(fEmax_adj_models);
    theeInverseComptonModel->SetLowEnergyLimit(fEmin_adj_models);
    theeInverseComptonModel->SetDirectProcess(theComptonScattering);
    theeInverseComptonModel->SetUseMatrix(false);

    theeInverseComptonModel->SetCSBiasingFactor( fCS_biasing_factor_compton);
    if (!fUse_forced_interaction) theeInverseComptonProjToProjCase =
              new G4eInverseCompton(true,"Inv_Compt",theeInverseComptonModel);
    theeInverseComptonProdToProjCase = new G4eInverseCompton(false,"Inv_Compt1",
                                          theeInverseComptonModel);
    if (fUse_forced_interaction){
      if (!theForcedInteractionForGamma )  theForcedInteractionForGamma =
        new G4AdjointForcedInteractionForGamma("ReverseGammaForcedInteraction");
      theForcedInteractionForGamma->
                     RegisterAdjointComptonModel(theeInverseComptonModel);
    }
    theAdjointSimManager->ConsiderParticleAsPrimary(G4String("e-"));
    theAdjointSimManager->ConsiderParticleAsPrimary(G4String("gamma"));
  }

  //Define  adjoint PEEffect
  //---------------------
  G4AdjointPhotoElectricModel* theInversePhotoElectricModel = 0;
  G4InversePEEffect* theInversePhotoElectricProcess = 0;
  
  if (fUse_peeffect) { 
    theInversePhotoElectricModel = new G4AdjointPhotoElectricModel();
    theInversePhotoElectricModel->SetHighEnergyLimit(fEmax_adj_models);
    theInversePhotoElectricModel->SetLowEnergyLimit(fEmin_adj_models);
    theInversePhotoElectricModel->SetCSBiasingFactor(
                                  fCS_biasing_factor_PEeffect);
    theInversePhotoElectricProcess = new G4InversePEEffect("Inv_PEEffect",
                                       theInversePhotoElectricModel);
    theAdjointSimManager->ConsiderParticleAsPrimary(G4String("e-"));
    theAdjointSimManager->ConsiderParticleAsPrimary(G4String("gamma"));
  }
  

  //Define  adjoint ionisation for protons
  //---------------------
   G4AdjointhIonisationModel* thepInverseIonisationModel = 0;
   G4hInverseIonisation* thepInverseIonisationProjToProjCase = 0 ;
   G4hInverseIonisation* thepInverseIonisationProdToProjCase = 0;
   if (fUse_pionisation) {
     thepInverseIonisationModel = new G4AdjointhIonisationModel(
                                           G4Proton::Proton());
     thepInverseIonisationModel->SetHighEnergyLimit(fEmax_adj_models);
     thepInverseIonisationModel->SetLowEnergyLimit(fEmin_adj_models);
     thepInverseIonisationModel->SetUseMatrix(false);
     thepInverseIonisationProjToProjCase = new G4hInverseIonisation(true,
                                     "Inv_pIon",thepInverseIonisationModel);
     thepInverseIonisationProdToProjCase = new G4hInverseIonisation(false,
                                    "Inv_pIon1",thepInverseIonisationModel);
     theAdjointSimManager->ConsiderParticleAsPrimary(G4String("e-"));
     theAdjointSimManager->ConsiderParticleAsPrimary(G4String("proton"));
  }

  //Declare the processes active for the different particles
  //--------------------------------------------------------
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (!pmanager) {
     pmanager = new G4ProcessManager(particle);
     particle->SetProcessManager(pmanager);
    }
        
    G4String particleName = particle->GetParticleName();
    if (particleName == "e-") {
      if (fUse_ms && fUse_eionisation) pmanager->AddProcess(theeminusMS);
      if (fUse_eionisation){
        pmanager->AddProcess(fEminusIonisation);
        G4AdjointCSManager::GetAdjointCSManager()->
                    RegisterEnergyLossProcess(fEminusIonisation,particle);
      }
      if (fUse_brem && fUse_eionisation) {
        pmanager->AddProcess(theeminusBremsstrahlung);
        G4AdjointCSManager::GetAdjointCSManager()->
                   RegisterEnergyLossProcess(theeminusBremsstrahlung,particle);
      }
      G4int n_order=0;
      if (fUse_ms && fUse_eionisation) {
        n_order++;
        pmanager->SetProcessOrdering(theeminusMS, idxAlongStep,n_order);
      }
      if (fUse_eionisation) {
        n_order++;
        pmanager->SetProcessOrdering(fEminusIonisation,idxAlongStep,n_order);
      }
      if (fUse_brem && fUse_eionisation) {
        n_order++;
        pmanager->SetProcessOrdering(theeminusBremsstrahlung,
                                           idxAlongStep,n_order);
      }
      n_order=0;
      if (fUse_ms && fUse_eionisation) {
        n_order++;
        pmanager->SetProcessOrdering(theeminusMS,idxPostStep,n_order);
      }
      if (fUse_eionisation) {
        n_order++;
        pmanager->SetProcessOrdering(fEminusIonisation,idxPostStep,n_order);
      }
      if (fUse_brem && fUse_eionisation) {
        n_order++;
        pmanager->SetProcessOrdering(theeminusBremsstrahlung,idxPostStep,
                                                                     n_order);
          }
      }

      if (particleName == "adj_e-") {
        G4ContinuousGainOfEnergy* theContinuousGainOfEnergy =0;
        if (fUse_eionisation ) {
          theContinuousGainOfEnergy= new G4ContinuousGainOfEnergy();
          theContinuousGainOfEnergy->SetLossFluctuations(
                                   fUse_egain_fluctuation);
          theContinuousGainOfEnergy->SetDirectEnergyLossProcess(
                                              fEminusIonisation);
          theContinuousGainOfEnergy->SetDirectParticle(G4Electron::Electron());
          pmanager->AddProcess(theContinuousGainOfEnergy);
          }
          G4int n_order=0;
          if (fUse_ms) {
            n_order++;
            pmanager->AddProcess(theeminusAdjointMS);
            pmanager->SetProcessOrdering(theeminusAdjointMS,
                                                idxAlongStep,n_order);
          }
          n_order++;
          pmanager->SetProcessOrdering(theContinuousGainOfEnergy,idxAlongStep,
                                                                    n_order);

          n_order++;
          G4AdjointAlongStepWeightCorrection* theAlongStepWeightCorrection =
                                      new G4AdjointAlongStepWeightCorrection();
          pmanager->AddProcess(theAlongStepWeightCorrection);
          pmanager->SetProcessOrdering(theAlongStepWeightCorrection,
                                                            idxAlongStep,
                                                                 n_order);
          n_order=0;
          if (fUse_eionisation) {
            pmanager->AddProcess(theeInverseIonisationProjToProjCase);
            pmanager->AddProcess(theeInverseIonisationProdToProjCase);
            n_order++;
            pmanager->SetProcessOrdering(theeInverseIonisationProjToProjCase,
                                                          idxPostStep,n_order);
            n_order++;
            pmanager->SetProcessOrdering(theeInverseIonisationProdToProjCase,
                                                         idxPostStep,n_order);
          }
          if (fUse_brem && fUse_eionisation) {
            pmanager->AddProcess(theeInverseBremsstrahlungProjToProjCase);
            n_order++;
            pmanager->SetProcessOrdering(
                              theeInverseBremsstrahlungProjToProjCase,
                                                       idxPostStep,n_order);
          }

          if (fUse_compton) {
            pmanager->AddProcess(theeInverseComptonProdToProjCase);
            n_order++;
            pmanager->SetProcessOrdering(theeInverseComptonProdToProjCase,
                                                         idxPostStep,n_order);
          }
          if (fUse_peeffect) {
            pmanager->AddDiscreteProcess(theInversePhotoElectricProcess);
            n_order++;
            pmanager->SetProcessOrdering(theInversePhotoElectricProcess,
                                                           idxPostStep,n_order);
          }
          if (fUse_pionisation) {
            pmanager->AddProcess(thepInverseIonisationProdToProjCase);
            n_order++;
            pmanager->SetProcessOrdering(thepInverseIonisationProdToProjCase,
                                                         idxPostStep,n_order);
          }
          if (fUse_ms && fUse_eionisation) {
            n_order++;
            pmanager->SetProcessOrdering(theeminusAdjointMS,
                                                        idxPostStep,n_order);
          }
        }
        
           
        if(particleName == "adj_gamma") {
          G4int n_order=0;
          if (!fUse_forced_interaction){
           G4AdjointAlongStepWeightCorrection* theAlongStepWeightCorrection =
                                      new G4AdjointAlongStepWeightCorrection();
           pmanager->AddProcess(theAlongStepWeightCorrection);
           pmanager->SetProcessOrdering(theAlongStepWeightCorrection,
                                                       idxAlongStep,1);
                
           if (fUse_brem && fUse_eionisation) {
            pmanager->AddProcess(theeInverseBremsstrahlungProdToProjCase);
            n_order++;
            pmanager->SetProcessOrdering(
                                  theeInverseBremsstrahlungProdToProjCase,
                                                   idxPostStep,n_order);
            }
           if (fUse_compton) {
            pmanager->AddDiscreteProcess(theeInverseComptonProjToProjCase);
            n_order++;
            pmanager->SetProcessOrdering(theeInverseComptonProjToProjCase,
                                                        idxPostStep,n_order);
           }
          }
          else {
           if (theForcedInteractionForGamma) {
            pmanager->AddProcess(theForcedInteractionForGamma);
            n_order++;
            pmanager->SetProcessOrdering(theForcedInteractionForGamma,
                                                       idxPostStep,n_order);
            pmanager->SetProcessOrdering(theForcedInteractionForGamma,
                                                      idxAlongStep,n_order);
           }
          }
        } 
   
        if (particleName == "gamma") {
          if (fUse_compton) {
            pmanager->AddDiscreteProcess(theComptonScattering);
                G4AdjointCSManager::GetAdjointCSManager()->
                             RegisterEmProcess(theComptonScattering,particle);
          }
          if (fUse_peeffect) {
                pmanager->AddDiscreteProcess(thePEEffect);
                G4AdjointCSManager::GetAdjointCSManager()->
                                       RegisterEmProcess(thePEEffect,particle);
          }
      if (fUse_gamma_conversion) {
            pmanager->AddDiscreteProcess(theGammaConversion);
      }
     }
        
     if (particleName == "e+" && fUse_gamma_conversion) {//positron
      G4VProcess* theeplusMultipleScattering = new G4eMultipleScattering();
      G4VProcess* theeplusIonisation         = new G4eIonisation();
      G4VProcess* theeplusBremsstrahlung     = new G4eBremsstrahlung();
      G4VProcess* theeplusAnnihilation       = new G4eplusAnnihilation();

      // add processes
      pmanager->AddProcess(theeplusMultipleScattering);
      pmanager->AddProcess(theeplusIonisation);
      pmanager->AddProcess(theeplusBremsstrahlung);
      pmanager->AddProcess(theeplusAnnihilation);

      // set ordering for AtRestDoIt
      pmanager->SetProcessOrderingToFirst(theeplusAnnihilation, idxAtRest);

      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(theeplusMultipleScattering,
                                                             idxAlongStep,1);
      pmanager->SetProcessOrdering(theeplusIonisation, idxAlongStep,2);
      pmanager->SetProcessOrdering(theeplusBremsstrahlung,idxAlongStep,3);

      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(theeplusMultipleScattering,
                                                             idxPostStep,1);
      pmanager->SetProcessOrdering(theeplusIonisation,idxPostStep,2);
      pmanager->SetProcessOrdering(theeplusBremsstrahlung,idxPostStep,3);
      pmanager->SetProcessOrdering(theeplusAnnihilation,idxPostStep,4);
        }
        if (particleName == "proton" && fUse_pionisation) {
          if (fUse_ms && fUse_pionisation) pmanager->AddProcess(thepMS);

          if (fUse_pionisation){
                pmanager->AddProcess(fPIonisation);
                G4AdjointCSManager::GetAdjointCSManager()->
                          RegisterEnergyLossProcess(fPIonisation,particle);
          }

          G4int n_order=0;
          if (fUse_ms && fUse_pionisation) {
            n_order++;
                pmanager->SetProcessOrdering(thepMS, idxAlongStep,n_order);
          }

          if (fUse_pionisation) {
              n_order++;
              pmanager->SetProcessOrdering(fPIonisation,idxAlongStep,n_order);
          }
                
          n_order=0;
          if (fUse_ms && fUse_pionisation) {
                n_order++;
              pmanager->SetProcessOrdering(thepMS, idxPostStep,n_order);
          }

          if (fUse_pionisation) {
              n_order++;
              pmanager->SetProcessOrdering(fPIonisation,idxPostStep,n_order);
          }
        
        }
        
        if (particleName == "adj_proton" && fUse_pionisation) {
          G4ContinuousGainOfEnergy* theContinuousGainOfEnergy =0;
          if (fUse_pionisation ) {
           theContinuousGainOfEnergy= new G4ContinuousGainOfEnergy();
           theContinuousGainOfEnergy->SetLossFluctuations(
                                   fUse_egain_fluctuation);
           theContinuousGainOfEnergy->SetDirectEnergyLossProcess(fPIonisation);
              theContinuousGainOfEnergy->SetDirectParticle(G4Proton::Proton());
              pmanager->AddProcess(theContinuousGainOfEnergy);
           }

           G4int n_order=0;
           if (fUse_ms) {
             n_order++;
             pmanager->AddProcess(thepMS);
             pmanager->SetProcessOrdering(thepMS, idxAlongStep,n_order);
           }

           n_order++;
           pmanager->SetProcessOrdering(theContinuousGainOfEnergy,idxAlongStep,
                                                                       n_order);
                
           n_order++;
           G4AdjointAlongStepWeightCorrection* theAlongStepWeightCorrection =
                                      new G4AdjointAlongStepWeightCorrection();
          pmanager->AddProcess(theAlongStepWeightCorrection);
          pmanager->SetProcessOrdering(theAlongStepWeightCorrection,
                                                idxAlongStep,
                                                        n_order);
          n_order=0;
          if (fUse_pionisation) {
            pmanager->AddProcess(thepInverseIonisationProjToProjCase);
            n_order++;
                pmanager->SetProcessOrdering(
                                   thepInverseIonisationProjToProjCase,
                                                       idxPostStep,n_order);
          }

          if (fUse_ms && fUse_pionisation) {
            n_order++;
            pmanager->SetProcessOrdering(thepMS,idxPostStep,n_order);
          }
        }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Decay.hh"
void G4AdjointPhysicsList::ConstructGeneral()
{
  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) { 
      pmanager ->AddProcess(theDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4AdjointPhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "G4AdjointPhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  //
  SetCutValue(defaultCutValue, "gamma");
  SetCutValue(defaultCutValue, "e-");
  SetCutValue(defaultCutValue, "e+");

  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
