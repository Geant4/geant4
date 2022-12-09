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
// G4TransportationWithMsc
//
// Class Description:
//
// It is a generic process of transportation with multiple scattering included
// in the step limitation and propagation.
//
// Original author: Jonas Hahnfeld, 2022

// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4TransportationWithMsc.hh"

#include "G4LossTableBuilder.hh"
#include "G4LossTableManager.hh"
#include "G4EmConfigurator.hh"
#include "G4VMscModel.hh"

#include "G4DynamicParticle.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4StepStatus.hh"
#include "G4Track.hh"

#include "G4Electron.hh"

#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

static constexpr G4double kLowestKinEnergy  = 10 * CLHEP::eV;
static constexpr G4double kGeomMin          = 0.05 * CLHEP::nm;
static constexpr G4double kMinDisplacement2 = kGeomMin * kGeomMin;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4TransportationWithMsc::G4TransportationWithMsc(ScatteringType type,
                                                 G4int verbosity)
  : G4Transportation(verbosity, "TransportationWithMsc")
  , fType(type)
{
  SetVerboseLevel(1);

  fEmManager    = G4LossTableManager::Instance();
  fModelManager = new G4EmModelManager;

  G4ThreeVector zero;
  fSubStepDynamicParticle =
    new G4DynamicParticle(G4Electron::Definition(), zero);
  fSubStepTrack = new G4Track(fSubStepDynamicParticle, 0, zero);
  fSubStep      = new G4Step;
  fSubStepTrack->SetStep(fSubStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4TransportationWithMsc::~G4TransportationWithMsc()
{
  delete fModelManager;
  // fSubStepDynamicParticle is owned and also deleted by fSubStepTrack!
  delete fSubStepTrack;
  delete fSubStep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4TransportationWithMsc::AddMscModel(G4VMscModel* mscModel, G4int order,
                                          const G4Region* region)
{
  if(fType != ScatteringType::MultipleScattering)
  {
    G4Exception("G4TransportationWithMsc::AddMscModel", "em0051",
                FatalException,
                "not allowed unless type == MultipleScattering");
  }

  fModelManager->AddEmModel(order, mscModel, nullptr, region);
  mscModel->SetParticleChange(&fParticleChange);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4TransportationWithMsc::PreparePhysicsTable(
  const G4ParticleDefinition& part)
{
  if(nullptr == fFirstParticle)
  {
    fFirstParticle = &part;
    G4VMultipleScattering* ptr = nullptr;
    auto emConfigurator = fEmManager->EmConfigurator();
    emConfigurator->PrepareModels(&part, ptr, this);
  }

  if(fFirstParticle == &part)
  {
    G4bool master             = fEmManager->IsMaster();
    G4LossTableBuilder* bld   = fEmManager->GetTableBuilder();
    G4bool baseMat            = bld->GetBaseMaterialFlag();
    const auto* theParameters = G4EmParameters::Instance();

    if(master)
    {
      SetVerboseLevel(theParameters->Verbose());
    }
    else
    {
      SetVerboseLevel(theParameters->WorkerVerbose());
    }

    const G4int numberOfModels = fModelManager->NumberOfModels();
    for(G4int i = 0; i < numberOfModels; ++i)
    {
      auto msc = static_cast<G4VMscModel*>(fModelManager->GetModel(i));
      msc->SetMasterThread(master);
      msc->SetPolarAngleLimit(theParameters->MscThetaLimit());
      G4double emax =
        std::min(msc->HighEnergyLimit(), theParameters->MaxKinEnergy());
      msc->SetHighEnergyLimit(emax);
      msc->SetUseBaseMaterials(baseMat);
    }

    fModelManager->Initialise(fFirstParticle, G4Electron::Electron(),
                              verboseLevel);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4TransportationWithMsc::BuildPhysicsTable(
  const G4ParticleDefinition& part)
{
  if(fFirstParticle == &part)
  {
    fEmManager->BuildPhysicsTable(fFirstParticle);

    if(!fEmManager->IsMaster())
    {
      const auto masterProcess =
        static_cast<const G4TransportationWithMsc*>(GetMasterProcess());

      // Initialisation of models.
      const G4int numberOfModels = fModelManager->NumberOfModels();
      for(G4int i = 0; i < numberOfModels; ++i)
      {
        auto msc = static_cast<G4VMscModel*>(fModelManager->GetModel(i));
        auto msc0 =
          static_cast<G4VMscModel*>(masterProcess->fModelManager->GetModel(i));
        msc->SetCrossSectionTable(msc0->GetCrossSectionTable(), false);
        msc->InitialiseLocal(fFirstParticle, msc0);
      }
    }
  }

  if(!G4EmParameters::Instance()->IsPrintLocked() && verboseLevel > 0)
  {
    G4cout << G4endl;
    G4cout << GetProcessName() << ": for " << part.GetParticleName();
    if(fMultipleSteps)
    {
      G4cout << " (multipleSteps: 1)";
    }
    G4cout << G4endl;
    fModelManager->DumpModelList(G4cout, verboseLevel);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4TransportationWithMsc::StartTracking(G4Track* track)
{
  auto* currParticle = track->GetParticleDefinition();
  auto* ionisation   = fEmManager->GetEnergyLossProcess(currParticle);

  fSubStepDynamicParticle->SetDefinition(currParticle);

  const G4int numberOfModels = fModelManager->NumberOfModels();
  for(G4int i = 0; i < numberOfModels; ++i)
  {
    auto msc = static_cast<G4VMscModel*>(fModelManager->GetModel(i));
    msc->StartTracking(track);
    msc->SetIonisation(ionisation, currParticle);
  }

  // Ensure that field propagation state is also cleared / prepared
  G4Transportation::StartTracking(track);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4TransportationWithMsc::AlongStepGetPhysicalInteractionLength(
  const G4Track& track, G4double previousStepSize, G4double currentMinimumStep,
  G4double& proposedSafety, G4GPILSelection* selection)
{
  *selection = NotCandidateForSelection;

  const G4double physStepLimit = currentMinimumStep;

  switch(fType)
  {
    case ScatteringType::MultipleScattering: {
      // Select the MSC model for the current kinetic energy.
      G4VMscModel* mscModel          = nullptr;
      const G4double ekin            = track.GetKineticEnergy();
      const auto* couple             = track.GetMaterialCutsCouple();
      const auto* particleDefinition = track.GetParticleDefinition();
      if(physStepLimit > kGeomMin)
      {
        G4double ekinForSelection = ekin;
        G4double pdgMass          = particleDefinition->GetPDGMass();
        if(pdgMass > CLHEP::GeV)
        {
          ekinForSelection *= proton_mass_c2 / pdgMass;
        }

        if(ekinForSelection >= kLowestKinEnergy)
        {
          mscModel = static_cast<G4VMscModel*>(
            fModelManager->SelectModel(ekinForSelection, couple->GetIndex()));
          if(mscModel == nullptr)
          {
            G4Exception("G4TransportationWithMsc::AlongStepGPIL", "em0052",
                        FatalException, "no MSC model found");
          }
          if(!mscModel->IsActive(ekinForSelection))
          {
            mscModel = nullptr;
          }
        }
      }

      // Call the MSC model to potentially limit the step and convert to
      // geometric path length.
      if(mscModel != nullptr)
      {
        mscModel->SetCurrentCouple(couple);

        // Use the provided track for the first step.
        const G4Track* currentTrackPtr = &track;

        G4double currentSafety = proposedSafety;
        G4double currentEnergy = ekin;

        G4double stepLimitLeft           = physStepLimit;
        G4double totalGeometryStepLength = 0, totalTruePathLength = 0;
        G4bool firstStep = true, continueStepping = fMultipleSteps;

        do
        {
          G4double gPathLength = stepLimitLeft;
          G4double tPathLength =
            mscModel->ComputeTruePathLengthLimit(*currentTrackPtr, gPathLength);
          G4bool mscLimitsStep = (tPathLength < stepLimitLeft);
          if(!fMultipleSteps && mscLimitsStep)
          {
            // MSC limits the step.
            *selection = CandidateForSelection;
          }

          if(!firstStep)
          {
            // Move the navigator to where the previous step ended.
            fLinearNavigator->LocateGlobalPointWithinVolume(
              fTransportEndPosition);
          }

          G4GPILSelection transportSelection;
          G4double geometryStepLength =
            G4Transportation::AlongStepGetPhysicalInteractionLength(
              *currentTrackPtr, previousStepSize, gPathLength, currentSafety,
              &transportSelection);
          if(geometryStepLength < gPathLength)
          {
            // Transportation limits the step, ie the track hit a boundary.
            *selection       = CandidateForSelection;
            continueStepping = false;
          }
          if(fTransportEndKineticEnergy != currentEnergy)
          {
            // Field propagation changed the energy, it's not possible to
            // estimate the continuous energy loss and continue stepping.
            continueStepping = false;
          }

          if(firstStep)
          {
            proposedSafety = currentSafety;
          }
          totalGeometryStepLength += geometryStepLength;

          // Sample MSC direction change and displacement.
          const G4double range =
            mscModel->GetRange(particleDefinition, currentEnergy, couple);

          tPathLength = mscModel->ComputeTrueStepLength(geometryStepLength);

          // Protect against wrong t->g->t conversion.
          tPathLength = std::min(tPathLength, stepLimitLeft);

          totalTruePathLength += tPathLength;
          if(*selection != CandidateForSelection && !mscLimitsStep)
          {
            // If neither MSC nor transportation limits the step, we got the
            // distance we want - make sure we exit the loop.
            continueStepping = false;
          }
          else if(tPathLength >= range)
          {
            // The particle will stop, exit the loop.
            continueStepping = false;
          }
          else
          {
            stepLimitLeft -= tPathLength;
          }

          // Do not sample scattering at the last or at a small step.
          if(tPathLength < range && tPathLength > kGeomMin)
          {
            static constexpr G4double minSafety = 1.20 * CLHEP::nm;
            static constexpr G4double sFact     = 0.99;

            // The call to SampleScattering() *may* directly fill in the changed
            // direction into fParticleChange, so we have to:
            // 1) Make sure the momentum direction is initialized.
            fParticleChange.ProposeMomentumDirection(fTransportEndMomentumDir);
            // 2) Call SampleScattering(), which *may* change it.
            const G4ThreeVector displacement =
              mscModel->SampleScattering(fTransportEndMomentumDir, minSafety);
            // 3) Get the changed direction and inform G4Transportation.
            fMomentumChanged         = true;
            fTransportEndMomentumDir = *fParticleChange.GetMomentumDirection();

            const G4double r2 = displacement.mag2();
            if(r2 > kMinDisplacement2)
            {
              G4bool positionChanged = true;
              G4double dispR         = std::sqrt(r2);
              G4double postSafety    = sFact * fpSafetyHelper->ComputeSafety(
                                                 fTransportEndPosition, dispR);

              // Far away from geometry boundary
              if(postSafety > 0.0 && dispR <= postSafety)
              {
                fTransportEndPosition += displacement;

                // Near the boundary
              }
              else
              {
                // displaced point is definitely within the volume
                if(dispR < postSafety)
                {
                  fTransportEndPosition += displacement;

                  // reduced displacement
                }
                else if(postSafety > kGeomMin)
                {
                  fTransportEndPosition += displacement * (postSafety / dispR);

                  // very small postSafety
                }
                else
                {
                  positionChanged = false;
                }
              }
              if(positionChanged)
              {
                fpSafetyHelper->ReLocateWithinVolume(fTransportEndPosition);
              }
            }
          }

          if(continueStepping)
          {
            // Update safety according to the geometry distance.
            if(currentSafety < fEndPointDistance)
            {
              currentSafety = 0;
            }
            else
            {
              currentSafety -= fEndPointDistance;
            }

            // Update the kinetic energy according to the continuous loss.
            currentEnergy = mscModel->GetEnergy(particleDefinition,
                                                range - tPathLength, couple);

            // From now on, use the track that we can update below.
            currentTrackPtr = fSubStepTrack;

            fSubStepDynamicParticle->SetKineticEnergy(currentEnergy);
            fSubStepDynamicParticle->SetMomentumDirection(
              fTransportEndMomentumDir);
            fSubStepTrack->SetPosition(fTransportEndPosition);

            G4StepPoint& subPreStepPoint = *fSubStep->GetPreStepPoint();
            subPreStepPoint.SetMaterialCutsCouple(couple);
            subPreStepPoint.SetPosition(fTransportEndPosition);
            subPreStepPoint.SetSafety(currentSafety);
            subPreStepPoint.SetStepStatus(fAlongStepDoItProc);
          }
          firstStep = false;
        } while(continueStepping);

        // Note: currentEnergy is only updated if continueStepping is true.
        // In case field propagation changed the energy, this flag is
        // immediately set to false and currentEnergy is still equal to the
        // initial kinetic energy stored in ekin.
        if(currentEnergy != ekin)
        {
          // If field propagation didn't change the energy and we potentially
          // did multiple steps, reset the energy that G4Transportation will
          // propose to not subtract the energy loss twice.
          fTransportEndKineticEnergy = ekin;
          // Also ask for the range again with the initial energy so it is
          // correctly cached in the G4VEnergyLossProcess.
          // FIXME: Asking for a range should never change the cached values!
          (void) mscModel->GetRange(particleDefinition, ekin, couple);
        }

        fParticleChange.ProposeTrueStepLength(totalTruePathLength);

        return totalGeometryStepLength;
      }
    }
  }

  // If we get here, no scattering has happened.
  return G4Transportation::AlongStepGetPhysicalInteractionLength(
    track, previousStepSize, currentMinimumStep, proposedSafety, selection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
