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
//
/// \file optical/OpNovice2/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "HistoManager.hh"
#include "Run.hh"
#include "TrackInformation.hh"

#include "G4Cerenkov.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4Scintillation.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4ProcessManager.hh"
#include "G4Step.hh"
#include "G4SteppingManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::SteppingAction()
  : G4UserSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::~SteppingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SteppingAction::UserSteppingAction(const G4Step* step)
{
  static G4ParticleDefinition* opticalphoton =
    G4OpticalPhoton::OpticalPhotonDefinition();

  G4AnalysisManager* analysisMan = G4AnalysisManager::Instance();
  Run* run =
    static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  G4Track* track          = step->GetTrack();
  G4StepPoint* endPoint   = step->GetPostStepPoint();
  G4StepPoint* startPoint = step->GetPreStepPoint();

  const G4DynamicParticle* theParticle = track->GetDynamicParticle();
  const G4ParticleDefinition* particleDef =
    theParticle->GetParticleDefinition();

  TrackInformation* trackInfo =
    (TrackInformation*) (track->GetUserInformation());

  if(particleDef == opticalphoton)
  {
    const G4VProcess* pds = endPoint->GetProcessDefinedStep();
    G4String procname     = pds->GetProcessName();
    if(procname.compare("OpAbsorption") == 0)
    {
      run->AddOpAbsorption();
      if(trackInfo->GetIsFirstTankX())
      {
        run->AddOpAbsorptionPrior();
      }
    }
    else if(procname.compare("OpRayleigh") == 0)
    {
      run->AddRayleigh();
    }
    else if(procname.compare("OpWLS") == 0)
    {
      G4double en = track->GetKineticEnergy();
      run->AddWLSAbsorption();
      run->AddWLSAbsorptionEnergy(en);
      analysisMan->FillH1(4, en / eV);  // absorption energy
      // loop over secondaries, create statistics
      // const std::vector<const G4Track*>* secondaries =
      auto secondaries = step->GetSecondaryInCurrentStep();
      for(auto sec : *secondaries)
      {
        en = sec->GetKineticEnergy();
        run->AddWLSEmission();
        run->AddWLSEmissionEnergy(en);
        analysisMan->FillH1(5, en / eV);  // emission energy
        G4double time = sec->GetGlobalTime();
        analysisMan->FillH1(6, time / ns);
      }
    }
    else if(procname.compare("OpWLS2") == 0)
    {
      G4double en = track->GetKineticEnergy();
      run->AddWLS2Absorption();
      run->AddWLS2AbsorptionEnergy(en);
      analysisMan->FillH1(7, en / eV);  // absorption energy
      // loop over secondaries, create statistics
      // const std::vector<const G4Track*>* secondaries =
      auto secondaries = step->GetSecondaryInCurrentStep();
      for(auto sec : *secondaries)
      {
        en = sec->GetKineticEnergy();
        run->AddWLS2Emission();
        run->AddWLS2EmissionEnergy(en);
        analysisMan->FillH1(8, en / eV);  // emission energy
        G4double time = sec->GetGlobalTime();
        analysisMan->FillH1(9, time / ns);
      }
    }

    // optical process has endpt on bdry,
    if(endPoint->GetStepStatus() == fGeomBoundary)
    {
      G4ThreeVector p0 = startPoint->GetMomentumDirection();
      G4ThreeVector p1 = endPoint->GetMomentumDirection();

      G4OpBoundaryProcessStatus theStatus = Undefined;

      G4ProcessManager* OpManager = opticalphoton->GetProcessManager();
      G4ProcessVector* postStepDoItVector =
        OpManager->GetPostStepProcessVector(typeDoIt);
      G4int n_proc = postStepDoItVector->entries();

      if(trackInfo->GetIsFirstTankX())
      {
        G4double px1         = p1.x();
        G4double py1         = p1.y();
        G4double pz1         = p1.z();
        // do not count Absorbed or Detected photons here
        if(track->GetTrackStatus() != fStopAndKill)
        {
          if(px1 < 0.)
          {
            analysisMan->FillH1(11, px1);
            analysisMan->FillH1(12, py1);
            analysisMan->FillH1(13, pz1);
          }
          else
          {
            analysisMan->FillH1(14, px1);
            analysisMan->FillH1(15, py1);
            analysisMan->FillH1(16, pz1);
          }
        }

        trackInfo->SetIsFirstTankX(false);
        run->AddTotalSurface();

        for(G4int i = 0; i < n_proc; ++i)
        {
          G4VProcess* currentProcess = (*postStepDoItVector)[i];

          G4OpBoundaryProcess* opProc =
            dynamic_cast<G4OpBoundaryProcess*>(currentProcess);
          if(opProc)
          {
            G4double angle = std::acos(p0.x());
            theStatus      = opProc->GetStatus();
            analysisMan->FillH1(10, theStatus);
            switch(theStatus)
            {
              case Transmission:
                run->AddTransmission();
                analysisMan->FillH1(25, angle / deg);
                break;
              case FresnelRefraction:
                run->AddFresnelRefraction();
                analysisMan->FillH1(17, px1);
                analysisMan->FillH1(18, py1);
                analysisMan->FillH1(19, pz1);
                analysisMan->FillH1(20, angle / deg);
                break;
              case FresnelReflection:
                run->AddFresnelReflection();
                analysisMan->FillH1(21, angle / deg);
                analysisMan->FillH1(23, angle / deg);
                break;
              case TotalInternalReflection:
                run->AddTotalInternalReflection();
                analysisMan->FillH1(22, angle / deg);
                analysisMan->FillH1(23, angle / deg);
                break;
              case LambertianReflection:
                run->AddLambertianReflection();
                break;
              case LobeReflection:
                run->AddLobeReflection();
                break;
              case SpikeReflection:
                run->AddSpikeReflection();
                analysisMan->FillH1(26, angle / deg);
                break;
              case BackScattering:
                run->AddBackScattering();
                break;
              case Absorption:
                analysisMan->FillH1(24, angle / deg);
                run->AddAbsorption();
                break;
              case Detection:
                run->AddDetection();
                break;
              case NotAtBoundary:
                run->AddNotAtBoundary();
                break;
              case SameMaterial:
                run->AddSameMaterial();
                break;
              case StepTooSmall:
                run->AddStepTooSmall();
                break;
              case NoRINDEX:
                run->AddNoRINDEX();
                break;
              case PolishedLumirrorAirReflection:
                run->AddPolishedLumirrorAirReflection();
                break;
              case PolishedLumirrorGlueReflection:
                run->AddPolishedLumirrorGlueReflection();
                break;
              case PolishedAirReflection:
                run->AddPolishedAirReflection();
                break;
              case PolishedTeflonAirReflection:
                run->AddPolishedTeflonAirReflection();
                break;
              case PolishedTiOAirReflection:
                run->AddPolishedTiOAirReflection();
                break;
              case PolishedTyvekAirReflection:
                run->AddPolishedTyvekAirReflection();
                break;
              case PolishedVM2000AirReflection:
                run->AddPolishedVM2000AirReflection();
                break;
              case PolishedVM2000GlueReflection:
                run->AddPolishedVM2000AirReflection();
                break;
              case EtchedLumirrorAirReflection:
                run->AddEtchedLumirrorAirReflection();
                break;
              case EtchedLumirrorGlueReflection:
                run->AddEtchedLumirrorGlueReflection();
                break;
              case EtchedAirReflection:
                run->AddEtchedAirReflection();
                break;
              case EtchedTeflonAirReflection:
                run->AddEtchedTeflonAirReflection();
                break;
              case EtchedTiOAirReflection:
                run->AddEtchedTiOAirReflection();
                break;
              case EtchedTyvekAirReflection:
                run->AddEtchedTyvekAirReflection();
                break;
              case EtchedVM2000AirReflection:
                run->AddEtchedVM2000AirReflection();
                break;
              case EtchedVM2000GlueReflection:
                run->AddEtchedVM2000AirReflection();
                break;
              case GroundLumirrorAirReflection:
                run->AddGroundLumirrorAirReflection();
                break;
              case GroundLumirrorGlueReflection:
                run->AddGroundLumirrorGlueReflection();
                break;
              case GroundAirReflection:
                run->AddGroundAirReflection();
                break;
              case GroundTeflonAirReflection:
                run->AddGroundTeflonAirReflection();
                break;
              case GroundTiOAirReflection:
                run->AddGroundTiOAirReflection();
                break;
              case GroundTyvekAirReflection:
                run->AddGroundTyvekAirReflection();
                break;
              case GroundVM2000AirReflection:
                run->AddGroundVM2000AirReflection();
                break;
              case GroundVM2000GlueReflection:
                run->AddGroundVM2000AirReflection();
                break;
              case Dichroic:
                run->AddDichroic();
                break;
              case CoatedDielectricReflection:
                run->AddCoatedDielectricReflection();
                break;
              case CoatedDielectricRefraction:
                run->AddCoatedDielectricRefraction();
                break;
              case CoatedDielectricFrustratedTransmission:
                run->AddCoatedDielectricFrustratedTransmission();
                break;
              default:
                G4cout << "theStatus: " << theStatus
                       << " was none of the above." << G4endl;
                break;
            }
          }
        }
      }
    }
    // This block serves to test that G4OpBoundaryProcess sets the group
    // velocity correctly. It is not necessary to include in user code.
    // Only steps where pre- and post- are the same material, to avoid
    // incorrect checks (so, in practice, set e.g. OpRayleigh low enough
    // for particles to step in the interior of each volume.
    if (endPoint->GetMaterial() == startPoint->GetMaterial())
    {
      G4double trackVelocity = track->GetVelocity();
      G4double materialVelocity = CLHEP::c_light;
      G4MaterialPropertyVector* velVector = endPoint->GetMaterial()
        ->GetMaterialPropertiesTable()->GetProperty(kGROUPVEL);
      if(velVector)
      {
        materialVelocity =
          velVector->Value(theParticle->GetTotalMomentum(), fIdxVelocity);
      }

      if (std::abs(trackVelocity - materialVelocity) > 1e-9 * CLHEP::c_light)
      {
        G4ExceptionDescription ed;
        ed << "Optical photon group velocity: " << trackVelocity/(cm/ns)
           << " cm/ns is not what is expected from " << G4endl
           << "the material properties, "
           << materialVelocity/(cm/ns) << " cm/ns";
        G4Exception("OpNovice2 SteppingAction", "OpNovice2_1",
                    FatalException, ed);
      }
    }
    // end of group velocity test
  }

  else
  {  // particle != opticalphoton
    // print how many Cerenkov and scint photons produced this step
    // this demonstrates use of GetNumPhotons()
    auto proc_man =
      track->GetDynamicParticle()->GetParticleDefinition()->GetProcessManager();
    G4ProcessVector* proc_vec = proc_man->GetPostStepProcessVector(typeDoIt);
    G4int n_proc              = proc_vec->entries();

    G4int n_scint = 0;
    G4int n_cer   = 0;
    for(G4int i = 0; i < n_proc; ++i)
    {
      G4String proc_name = (*proc_vec)[i]->GetProcessName();
      if(proc_name.compare("Cerenkov") == 0)
      {
        auto cer = (G4Cerenkov*) (*proc_vec)[i];
        n_cer    = cer->GetNumPhotons();
      }
      else if(proc_name.compare("Scintillation") == 0)
      {
        auto scint = (G4Scintillation*) (*proc_vec)[i];
        n_scint    = scint->GetNumPhotons();
      }
    }
    if(fVerbose > 0)
    {
      if(n_cer > 0 || n_scint > 0)
      {
        G4cout << "In this step, " << n_cer << " Cerenkov and " << n_scint
               << " scintillation photons were produced." << G4endl;
      }
    }

    // loop over secondaries, create statistics
    const std::vector<const G4Track*>* secondaries =
      step->GetSecondaryInCurrentStep();

    for(auto sec : *secondaries)
    {
      if(sec->GetDynamicParticle()->GetParticleDefinition() == opticalphoton)
      {
        G4String creator_process = sec->GetCreatorProcess()->GetProcessName();
        if(creator_process.compare("Cerenkov") == 0)
        {
          G4double en = sec->GetKineticEnergy();
          run->AddCerenkovEnergy(en);
          run->AddCerenkov();
          analysisMan->FillH1(1, en / eV);
        }
        else if(creator_process.compare("Scintillation") == 0)
        {
          G4double en = sec->GetKineticEnergy();
          run->AddScintillationEnergy(en);
          run->AddScintillation();
          analysisMan->FillH1(2, en / eV);

          G4double time = sec->GetGlobalTime();
          analysisMan->FillH1(3, time / ns);
        }
      }
    }
  }

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
