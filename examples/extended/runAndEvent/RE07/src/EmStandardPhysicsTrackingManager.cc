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
// Implementation of a custom tracking manager for e-/e+ and gamma, using
// the same processes as defined in G4EmStandardPhysics.
//
// Original author: Jonas Hahnfeld, 2021

#include "EmStandardPhysicsTrackingManager.hh"

#include "G4ComptonScattering.hh"
#include "G4CoulombScattering.hh"
#include "G4Electron.hh"
#include "G4EmParameters.hh"
#include "G4Gamma.hh"
#include "G4GammaConversion.hh"
#include "G4KleinNishinaModel.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4LivermorePolarizedRayleighModel.hh"
#include "G4PhotoElectricAngularGeneratorPolarized.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4Positron.hh"
#include "G4RayleighScattering.hh"
#include "G4SystemOfUnits.hh"
#include "G4UrbanMscModel.hh"
#include "G4WentzelVIModel.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eCoulombScatteringModel.hh"
#include "G4eIonisation.hh"
#include "G4eMultipleScattering.hh"
#include "G4eplusAnnihilation.hh"

#include "TrackingManagerHelper.hh"

EmStandardPhysicsTrackingManager* EmStandardPhysicsTrackingManager::fMasterTrackingManager =
  nullptr;

EmStandardPhysicsTrackingManager::EmStandardPhysicsTrackingManager()
{
  G4EmParameters* param = G4EmParameters::Instance();
  G4double highEnergyLimit = param->MscEnergyLimit();
  G4bool polar = param->EnablePolarisation();

  // e-
  {
    G4eMultipleScattering* msc = new G4eMultipleScattering;
    G4UrbanMscModel* msc1 = new G4UrbanMscModel;
    G4WentzelVIModel* msc2 = new G4WentzelVIModel;
    msc1->SetHighEnergyLimit(highEnergyLimit);
    msc2->SetLowEnergyLimit(highEnergyLimit);
    msc->SetEmModel(msc1);
    msc->SetEmModel(msc2);
    fElectronProcs.msc = msc;

    fElectronProcs.ioni = new G4eIonisation;
    fElectronProcs.brems = new G4eBremsstrahlung;

    G4CoulombScattering* ss = new G4CoulombScattering;
    G4eCoulombScatteringModel* ssm = new G4eCoulombScatteringModel;
    ssm->SetLowEnergyLimit(highEnergyLimit);
    ssm->SetActivationLowEnergyLimit(highEnergyLimit);
    ss->SetEmModel(ssm);
    ss->SetMinKinEnergy(highEnergyLimit);
    fElectronProcs.ss = ss;
  }

  // e+
  {
    G4eMultipleScattering* msc = new G4eMultipleScattering;
    G4UrbanMscModel* msc1 = new G4UrbanMscModel;
    G4WentzelVIModel* msc2 = new G4WentzelVIModel;
    msc1->SetHighEnergyLimit(highEnergyLimit);
    msc2->SetLowEnergyLimit(highEnergyLimit);
    msc->SetEmModel(msc1);
    msc->SetEmModel(msc2);
    fPositronProcs.msc = msc;

    fPositronProcs.ioni = new G4eIonisation;
    fPositronProcs.brems = new G4eBremsstrahlung;
    fPositronProcs.annihilation = new G4eplusAnnihilation;

    G4CoulombScattering* ss = new G4CoulombScattering;
    G4eCoulombScatteringModel* ssm = new G4eCoulombScatteringModel;
    ssm->SetLowEnergyLimit(highEnergyLimit);
    ssm->SetActivationLowEnergyLimit(highEnergyLimit);
    ss->SetEmModel(ssm);
    ss->SetMinKinEnergy(highEnergyLimit);
    fPositronProcs.ss = ss;
  }

  {
    G4PhotoElectricEffect* pe = new G4PhotoElectricEffect;
    G4VEmModel* peModel = new G4LivermorePhotoElectricModel;
    if (polar) {
      peModel->SetAngularDistribution(new G4PhotoElectricAngularGeneratorPolarized);
    }
    pe->SetEmModel(peModel);
    fGammaProcs.pe = pe;

    G4ComptonScattering* cs = new G4ComptonScattering;
    if (polar) {
      cs->SetEmModel(new G4KleinNishinaModel);
    }
    fGammaProcs.compton = cs;

    fGammaProcs.conversion = new G4GammaConversion;

    G4RayleighScattering* rl = new G4RayleighScattering;
    if (polar) {
      rl->SetEmModel(new G4LivermorePolarizedRayleighModel);
    }
    fGammaProcs.rayleigh = rl;
  }

  if (fMasterTrackingManager == nullptr) {
    fMasterTrackingManager = this;
  }
  else {
    fElectronProcs.msc->SetMasterProcess(fMasterTrackingManager->fElectronProcs.msc);
    fElectronProcs.ss->SetMasterProcess(fMasterTrackingManager->fElectronProcs.ss);
    fElectronProcs.ioni->SetMasterProcess(fMasterTrackingManager->fElectronProcs.ioni);
    fElectronProcs.brems->SetMasterProcess(fMasterTrackingManager->fElectronProcs.brems);

    fPositronProcs.msc->SetMasterProcess(fMasterTrackingManager->fPositronProcs.msc);
    fPositronProcs.ss->SetMasterProcess(fMasterTrackingManager->fPositronProcs.ss);
    fPositronProcs.ioni->SetMasterProcess(fMasterTrackingManager->fPositronProcs.ioni);
    fPositronProcs.brems->SetMasterProcess(fMasterTrackingManager->fPositronProcs.brems);
    fPositronProcs.annihilation->SetMasterProcess(
      fMasterTrackingManager->fPositronProcs.annihilation);

    fGammaProcs.pe->SetMasterProcess(fMasterTrackingManager->fGammaProcs.pe);
    fGammaProcs.compton->SetMasterProcess(fMasterTrackingManager->fGammaProcs.compton);
    fGammaProcs.conversion->SetMasterProcess(fMasterTrackingManager->fGammaProcs.conversion);
    fGammaProcs.rayleigh->SetMasterProcess(fMasterTrackingManager->fGammaProcs.rayleigh);
  }
}

EmStandardPhysicsTrackingManager::~EmStandardPhysicsTrackingManager()
{
  if (fMasterTrackingManager == this) {
    fMasterTrackingManager = nullptr;
  }
}

void EmStandardPhysicsTrackingManager::BuildPhysicsTable(const G4ParticleDefinition& part)
{
  if (&part == G4Electron::Definition()) {
    fElectronProcs.msc->BuildPhysicsTable(part);
    fElectronProcs.ioni->BuildPhysicsTable(part);
    fElectronProcs.brems->BuildPhysicsTable(part);
    fElectronProcs.ss->BuildPhysicsTable(part);
  }
  else if (&part == G4Positron::Definition()) {
    fPositronProcs.msc->BuildPhysicsTable(part);
    fPositronProcs.ioni->BuildPhysicsTable(part);
    fPositronProcs.brems->BuildPhysicsTable(part);
    fPositronProcs.annihilation->BuildPhysicsTable(part);
    fPositronProcs.ss->BuildPhysicsTable(part);
  }
  else if (&part == G4Gamma::Definition()) {
    fGammaProcs.pe->BuildPhysicsTable(part);
    fGammaProcs.compton->BuildPhysicsTable(part);
    fGammaProcs.conversion->BuildPhysicsTable(part);
    fGammaProcs.rayleigh->BuildPhysicsTable(part);
  }
}

void EmStandardPhysicsTrackingManager::PreparePhysicsTable(const G4ParticleDefinition& part)
{
  if (&part == G4Electron::Definition()) {
    fElectronProcs.msc->PreparePhysicsTable(part);
    fElectronProcs.ioni->PreparePhysicsTable(part);
    fElectronProcs.brems->PreparePhysicsTable(part);
    fElectronProcs.ss->PreparePhysicsTable(part);
  }
  else if (&part == G4Positron::Definition()) {
    fPositronProcs.msc->PreparePhysicsTable(part);
    fPositronProcs.ioni->PreparePhysicsTable(part);
    fPositronProcs.brems->PreparePhysicsTable(part);
    fPositronProcs.annihilation->PreparePhysicsTable(part);
    fPositronProcs.ss->PreparePhysicsTable(part);
  }
  else if (&part == G4Gamma::Definition()) {
    fGammaProcs.pe->PreparePhysicsTable(part);
    fGammaProcs.compton->PreparePhysicsTable(part);
    fGammaProcs.conversion->PreparePhysicsTable(part);
    fGammaProcs.rayleigh->PreparePhysicsTable(part);
  }
}

void EmStandardPhysicsTrackingManager::TrackElectron(G4Track* aTrack)
{
  class ElectronPhysics final : public TrackingManagerHelper::Physics
  {
   public:
    ElectronPhysics(EmStandardPhysicsTrackingManager& mgr) : fMgr(mgr) {}

    void StartTracking(G4Track* aTrack) override
    {
      auto& electronProcs = fMgr.fElectronProcs;

      electronProcs.msc->StartTracking(aTrack);
      electronProcs.ioni->StartTracking(aTrack);
      electronProcs.brems->StartTracking(aTrack);
      electronProcs.ss->StartTracking(aTrack);

      fPreviousStepLength = 0;
    }
    void EndTracking() override
    {
      auto& electronProcs = fMgr.fElectronProcs;

      electronProcs.msc->EndTracking();
      electronProcs.ioni->EndTracking();
      electronProcs.brems->EndTracking();
      electronProcs.ss->EndTracking();
    }

    G4double GetPhysicalInteractionLength(const G4Track& track) override
    {
      auto& electronProcs = fMgr.fElectronProcs;
      G4double physIntLength, proposedSafety = DBL_MAX;
      G4ForceCondition condition;
      G4GPILSelection selection;

      fProposedStep = DBL_MAX;
      fSelected = -1;

      physIntLength = electronProcs.ss->PostStepGPIL(track, fPreviousStepLength, &condition);
      if (physIntLength < fProposedStep) {
        fProposedStep = physIntLength;
        fSelected = 0;
      }

      physIntLength = electronProcs.brems->PostStepGPIL(track, fPreviousStepLength, &condition);
      if (physIntLength < fProposedStep) {
        fProposedStep = physIntLength;
        fSelected = 1;
      }

      physIntLength = electronProcs.ioni->PostStepGPIL(track, fPreviousStepLength, &condition);
      if (physIntLength < fProposedStep) {
        fProposedStep = physIntLength;
        fSelected = 2;
      }

      physIntLength = electronProcs.ioni->AlongStepGPIL(
        track, fPreviousStepLength, fProposedStep, proposedSafety, &selection);
      if (physIntLength < fProposedStep) {
        fProposedStep = physIntLength;
        fSelected = -1;
      }

      physIntLength = electronProcs.msc->AlongStepGPIL(
        track, fPreviousStepLength, fProposedStep, proposedSafety, &selection);
      if (physIntLength < fProposedStep) {
        fProposedStep = physIntLength;
        // Check if MSC actually wants to win, in most cases it only limits the
        // step size.
        if (selection == CandidateForSelection) {
          fSelected = -1;
        }
      }

      return fProposedStep;
    }

    void AlongStepDoIt(G4Track& track, G4Step& step, G4TrackVector&) override
    {
      if (step.GetStepLength() == fProposedStep) {
        step.GetPostStepPoint()->SetStepStatus(fAlongStepDoItProc);
      }
      else {
        // Remember that the step was limited by geometry.
        fSelected = -1;
      }
      auto& electronProcs = fMgr.fElectronProcs;
      G4VParticleChange* particleChange;

      particleChange = electronProcs.msc->AlongStepDoIt(track, step);
      particleChange->UpdateStepForAlongStep(&step);
      track.SetTrackStatus(particleChange->GetTrackStatus());
      particleChange->Clear();

      particleChange = electronProcs.ioni->AlongStepDoIt(track, step);
      particleChange->UpdateStepForAlongStep(&step);
      track.SetTrackStatus(particleChange->GetTrackStatus());
      particleChange->Clear();

      fPreviousStepLength = step.GetStepLength();
    }

    void PostStepDoIt(G4Track& track, G4Step& step, G4TrackVector& secondaries) override
    {
      if (fSelected < 0) {
        return;
      }
      step.GetPostStepPoint()->SetStepStatus(fPostStepDoItProc);

      auto& electronProcs = fMgr.fElectronProcs;
      G4VProcess* process = nullptr;
      G4VParticleChange* particleChange = nullptr;

      switch (fSelected) {
        case 0:
          process = electronProcs.ss;
          particleChange = electronProcs.ss->PostStepDoIt(track, step);
          break;
        case 1:
          process = electronProcs.brems;
          particleChange = electronProcs.brems->PostStepDoIt(track, step);
          break;
        case 2:
          process = electronProcs.ioni;
          particleChange = electronProcs.ioni->PostStepDoIt(track, step);
          break;
      }

      particleChange->UpdateStepForPostStep(&step);
      step.UpdateTrack();

      int numSecondaries = particleChange->GetNumberOfSecondaries();
      for (int i = 0; i < numSecondaries; i++) {
        G4Track* secondary = particleChange->GetSecondary(i);
        secondary->SetParentID(track.GetTrackID());
        secondary->SetCreatorProcess(process);
        secondaries.push_back(secondary);
      }

      track.SetTrackStatus(particleChange->GetTrackStatus());
      particleChange->Clear();
    }

   private:
    EmStandardPhysicsTrackingManager& fMgr;
    G4double fPreviousStepLength;
    G4double fProposedStep;
    G4int fSelected;
  };

  ElectronPhysics physics(*this);
  TrackingManagerHelper::TrackChargedParticle(aTrack, physics);
}

void EmStandardPhysicsTrackingManager::TrackPositron(G4Track* aTrack)
{
  class PositronPhysics final : public TrackingManagerHelper::Physics
  {
   public:
    PositronPhysics(EmStandardPhysicsTrackingManager& mgr) : fMgr(mgr) {}

    void StartTracking(G4Track* aTrack) override
    {
      auto& positronProcs = fMgr.fPositronProcs;

      positronProcs.msc->StartTracking(aTrack);
      positronProcs.ioni->StartTracking(aTrack);
      positronProcs.brems->StartTracking(aTrack);
      positronProcs.annihilation->StartTracking(aTrack);
      positronProcs.ss->StartTracking(aTrack);

      fPreviousStepLength = 0;
    }
    void EndTracking() override
    {
      auto& positronProcs = fMgr.fPositronProcs;

      positronProcs.msc->EndTracking();
      positronProcs.ioni->EndTracking();
      positronProcs.brems->EndTracking();
      positronProcs.annihilation->EndTracking();
      positronProcs.ss->EndTracking();
    }

    G4double GetPhysicalInteractionLength(const G4Track& track) override
    {
      auto& positronProcs = fMgr.fPositronProcs;
      G4double physIntLength, proposedSafety = DBL_MAX;
      G4ForceCondition condition;
      G4GPILSelection selection;

      fProposedStep = DBL_MAX;
      fSelected = -1;

      physIntLength = positronProcs.ss->PostStepGPIL(track, fPreviousStepLength, &condition);
      if (physIntLength < fProposedStep) {
        fProposedStep = physIntLength;
        fSelected = 0;
      }

      physIntLength =
        positronProcs.annihilation->PostStepGPIL(track, fPreviousStepLength, &condition);
      if (physIntLength < fProposedStep) {
        fProposedStep = physIntLength;
        fSelected = 1;
      }

      physIntLength = positronProcs.brems->PostStepGPIL(track, fPreviousStepLength, &condition);
      if (physIntLength < fProposedStep) {
        fProposedStep = physIntLength;
        fSelected = 2;
      }

      physIntLength = positronProcs.ioni->PostStepGPIL(track, fPreviousStepLength, &condition);
      if (physIntLength < fProposedStep) {
        fProposedStep = physIntLength;
        fSelected = 3;
      }

      physIntLength = positronProcs.ioni->AlongStepGPIL(
        track, fPreviousStepLength, fProposedStep, proposedSafety, &selection);
      if (physIntLength < fProposedStep) {
        fProposedStep = physIntLength;
        fSelected = -1;
      }

      physIntLength = positronProcs.msc->AlongStepGPIL(
        track, fPreviousStepLength, fProposedStep, proposedSafety, &selection);
      if (physIntLength < fProposedStep) {
        fProposedStep = physIntLength;
        // Check if MSC actually wants to win, in most cases it only limits the
        // step size.
        if (selection == CandidateForSelection) {
          fSelected = -1;
        }
      }

      return fProposedStep;
    }

    void AlongStepDoIt(G4Track& track, G4Step& step, G4TrackVector&) override
    {
      if (step.GetStepLength() == fProposedStep) {
        step.GetPostStepPoint()->SetStepStatus(fAlongStepDoItProc);
      }
      else {
        // Remember that the step was limited by geometry.
        fSelected = -1;
      }
      auto& positronProcs = fMgr.fPositronProcs;
      G4VParticleChange* particleChange;

      particleChange = positronProcs.msc->AlongStepDoIt(track, step);
      particleChange->UpdateStepForAlongStep(&step);
      track.SetTrackStatus(particleChange->GetTrackStatus());
      particleChange->Clear();

      particleChange = positronProcs.ioni->AlongStepDoIt(track, step);
      particleChange->UpdateStepForAlongStep(&step);
      track.SetTrackStatus(particleChange->GetTrackStatus());
      particleChange->Clear();

      fPreviousStepLength = step.GetStepLength();
    }

    void PostStepDoIt(G4Track& track, G4Step& step, G4TrackVector& secondaries) override
    {
      if (fSelected < 0) {
        return;
      }
      step.GetPostStepPoint()->SetStepStatus(fPostStepDoItProc);

      auto& positronProcs = fMgr.fPositronProcs;
      G4VProcess* process;
      G4VParticleChange* particleChange = nullptr;

      switch (fSelected) {
        case 0:
          process = positronProcs.ss;
          particleChange = positronProcs.ss->PostStepDoIt(track, step);
          break;
        case 1:
          process = positronProcs.annihilation;
          particleChange = positronProcs.annihilation->PostStepDoIt(track, step);
          break;
        case 2:
          process = positronProcs.brems;
          particleChange = positronProcs.brems->PostStepDoIt(track, step);
          break;
        case 3:
          process = positronProcs.ioni;
          particleChange = positronProcs.ioni->PostStepDoIt(track, step);
          break;
      }

      particleChange->UpdateStepForPostStep(&step);
      step.UpdateTrack();

      int numSecondaries = particleChange->GetNumberOfSecondaries();
      for (int i = 0; i < numSecondaries; i++) {
        G4Track* secondary = particleChange->GetSecondary(i);
        secondary->SetParentID(track.GetTrackID());
        secondary->SetCreatorProcess(process);
        secondaries.push_back(secondary);
      }

      track.SetTrackStatus(particleChange->GetTrackStatus());
      particleChange->Clear();
    }

    G4bool HasAtRestProcesses() override { return true; }

    void AtRestDoIt(G4Track& track, G4Step& step, G4TrackVector& secondaries) override
    {
      auto& positronProcs = fMgr.fPositronProcs;
      // Annihilate the positron at rest.
      G4VParticleChange* particleChange = positronProcs.annihilation->AtRestDoIt(track, step);
      particleChange->UpdateStepForAtRest(&step);
      step.UpdateTrack();

      int numSecondaries = particleChange->GetNumberOfSecondaries();
      for (int i = 0; i < numSecondaries; i++) {
        G4Track* secondary = particleChange->GetSecondary(i);
        secondary->SetParentID(track.GetTrackID());
        secondary->SetCreatorProcess(positronProcs.annihilation);
        secondaries.push_back(secondary);
      }

      track.SetTrackStatus(particleChange->GetTrackStatus());
      particleChange->Clear();
    }

   private:
    EmStandardPhysicsTrackingManager& fMgr;
    G4double fPreviousStepLength;
    G4double fProposedStep;
    G4int fSelected;
  };

  PositronPhysics physics(*this);
  TrackingManagerHelper::TrackChargedParticle(aTrack, physics);
}

void EmStandardPhysicsTrackingManager::TrackGamma(G4Track* aTrack)
{
  class GammaPhysics final : public TrackingManagerHelper::Physics
  {
   public:
    GammaPhysics(EmStandardPhysicsTrackingManager& mgr) : fMgr(mgr) {}

    void StartTracking(G4Track* aTrack) override
    {
      auto& gammaProcs = fMgr.fGammaProcs;

      gammaProcs.pe->StartTracking(aTrack);
      gammaProcs.compton->StartTracking(aTrack);
      gammaProcs.conversion->StartTracking(aTrack);
      gammaProcs.rayleigh->StartTracking(aTrack);

      fPreviousStepLength = 0;
    }
    void EndTracking() override
    {
      auto& gammaProcs = fMgr.fGammaProcs;

      gammaProcs.pe->EndTracking();
      gammaProcs.compton->EndTracking();
      gammaProcs.conversion->EndTracking();
      gammaProcs.rayleigh->EndTracking();
    }

    G4double GetPhysicalInteractionLength(const G4Track& track) override
    {
      auto& gammaProcs = fMgr.fGammaProcs;
      G4double physIntLength;
      G4ForceCondition condition;

      fProposedStep = DBL_MAX;
      fSelected = -1;

      physIntLength = gammaProcs.rayleigh->PostStepGPIL(track, fPreviousStepLength, &condition);
      if (physIntLength < fProposedStep) {
        fProposedStep = physIntLength;
        fSelected = 0;
      }

      physIntLength = gammaProcs.conversion->PostStepGPIL(track, fPreviousStepLength, &condition);
      if (physIntLength < fProposedStep) {
        fProposedStep = physIntLength;
        fSelected = 1;
      }

      physIntLength = gammaProcs.compton->PostStepGPIL(track, fPreviousStepLength, &condition);
      if (physIntLength < fProposedStep) {
        fProposedStep = physIntLength;
        fSelected = 2;
      }

      physIntLength = gammaProcs.pe->PostStepGPIL(track, fPreviousStepLength, &condition);
      if (physIntLength < fProposedStep) {
        fProposedStep = physIntLength;
        fSelected = 3;
      }

      return fProposedStep;
    }

    void AlongStepDoIt(G4Track&, G4Step& step, G4TrackVector&) override
    {
      if (step.GetStepLength() == fProposedStep) {
        step.GetPostStepPoint()->SetStepStatus(fAlongStepDoItProc);
      }
      else {
        // Remember that the step was limited by geometry.
        fSelected = -1;
      }
      fPreviousStepLength = step.GetStepLength();
    }

    void PostStepDoIt(G4Track& track, G4Step& step, G4TrackVector& secondaries) override
    {
      if (fSelected < 0) {
        return;
      }
      step.GetPostStepPoint()->SetStepStatus(fPostStepDoItProc);

      auto& gammaProcs = fMgr.fGammaProcs;
      G4VProcess* process = nullptr;
      G4VParticleChange* particleChange = nullptr;

      switch (fSelected) {
        case 0:
          process = gammaProcs.rayleigh;
          particleChange = gammaProcs.rayleigh->PostStepDoIt(track, step);
          break;
        case 1:
          process = gammaProcs.conversion;
          particleChange = gammaProcs.conversion->PostStepDoIt(track, step);
          break;
        case 2:
          process = gammaProcs.compton;
          particleChange = gammaProcs.compton->PostStepDoIt(track, step);
          break;
        case 3:
          process = gammaProcs.pe;
          particleChange = gammaProcs.pe->PostStepDoIt(track, step);
          break;
      }

      particleChange->UpdateStepForPostStep(&step);
      step.UpdateTrack();

      int numSecondaries = particleChange->GetNumberOfSecondaries();
      for (int i = 0; i < numSecondaries; i++) {
        G4Track* secondary = particleChange->GetSecondary(i);
        secondary->SetParentID(track.GetTrackID());
        secondary->SetCreatorProcess(process);
        secondaries.push_back(secondary);
      }

      track.SetTrackStatus(particleChange->GetTrackStatus());
      particleChange->Clear();
    }

   private:
    EmStandardPhysicsTrackingManager& fMgr;
    G4double fPreviousStepLength;
    G4double fProposedStep;
    G4int fSelected;
  };

  GammaPhysics physics(*this);
  TrackingManagerHelper::TrackNeutralParticle(aTrack, physics);
}

void EmStandardPhysicsTrackingManager::HandOverOneTrack(G4Track* aTrack)
{
  const G4ParticleDefinition* part = aTrack->GetParticleDefinition();

  if (part == G4Electron::Definition()) {
    TrackElectron(aTrack);
  }
  else if (part == G4Positron::Definition()) {
    TrackPositron(aTrack);
  }
  else if (part == G4Gamma::Definition()) {
    TrackGamma(aTrack);
  }

  aTrack->SetTrackStatus(fStopAndKill);
  delete aTrack;
}
