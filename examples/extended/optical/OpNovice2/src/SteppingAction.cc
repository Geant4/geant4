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
//#include "EventAction.hh"
#include "HistoManager.hh"
#include "TrackInformation.hh"
#include "Run.hh"

#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4SteppingManager.hh"
#include "G4RunManager.hh"
#include "G4ProcessManager.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::SteppingAction()
: G4UserSteppingAction(),
  fVerbose(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SteppingAction::UserSteppingAction(const G4Step* step)
{
  static G4ParticleDefinition* opticalphoton = 
              G4OpticalPhoton::OpticalPhotonDefinition();
  G4AnalysisManager* analysisMan = G4AnalysisManager::Instance();
  Run* run = static_cast<Run*>(
               G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  G4Track* track = step->GetTrack();
  G4StepPoint* endPoint   = step->GetPostStepPoint();
  G4StepPoint* startPoint = step->GetPreStepPoint();

  G4String particleName = track->GetDynamicParticle()->
                                 GetParticleDefinition()->GetParticleName();

  TrackInformation* trackInfo = 
                        (TrackInformation*)(track->GetUserInformation());

  if (particleName == "opticalphoton") {
    const G4VProcess* pds = endPoint->GetProcessDefinedStep();
    if (pds->GetProcessName() == "OpAbsorption") {
      run->AddOpAbsorption(); 
      if (trackInfo->GetIsFirstTankX()) {
        run->AddOpAbsorptionPrior();
      }
    } 
    else if (pds->GetProcessName() == "OpRayleigh") {
      run->AddRayleigh();
    }

    // optical process has endpt on bdry, 
    if (endPoint->GetStepStatus() == fGeomBoundary) {

      const G4DynamicParticle* theParticle = track->GetDynamicParticle();

      G4ThreeVector oldMomentumDir = theParticle->GetMomentumDirection();

      G4ThreeVector m0 = startPoint->GetMomentumDirection();
      G4ThreeVector m1 = endPoint->GetMomentumDirection();

      G4OpBoundaryProcessStatus theStatus = Undefined;

      G4ProcessManager* OpManager = 
        G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
      G4int MAXofPostStepLoops = 
        OpManager->GetPostStepProcessVector()->entries();
      G4ProcessVector* postStepDoItVector = 
        OpManager->GetPostStepProcessVector(typeDoIt);

      if (trackInfo->GetIsFirstTankX()) {
        G4ThreeVector momdir = endPoint->GetMomentumDirection();
        G4double px1 = momdir.x();
        G4double py1 = momdir.y();
        G4double pz1 = momdir.z();
        if (px1 < 0.) {
          analysisMan->FillH1(4, px1);
          analysisMan->FillH1(5, py1);
          analysisMan->FillH1(6, pz1);
        } else if (px1 >= 0.) {
          analysisMan->FillH1(7, px1);
          analysisMan->FillH1(8, py1);
          analysisMan->FillH1(9, pz1);
        }

        trackInfo->SetIsFirstTankX(false);
        run->AddTotalSurface(); 

        for (G4int i=0; i<MAXofPostStepLoops; ++i) {
          G4VProcess* currentProcess = (*postStepDoItVector)[i];

          G4OpBoundaryProcess* opProc = 
            dynamic_cast<G4OpBoundaryProcess*>(currentProcess);
          if (opProc) {
            theStatus = opProc->GetStatus();
            analysisMan->FillH1(3, theStatus);
            if (theStatus == Transmission) {
              run->AddTransmission();
            }
            else if (theStatus == FresnelRefraction) {
              run->AddFresnelRefraction(); 
              analysisMan->FillH1(10, px1);
              analysisMan->FillH1(11, py1);
              analysisMan->FillH1(12, pz1);
            }
            else if (theStatus == FresnelReflection) { 
              run->AddFresnelReflection(); 
            }
            else if (theStatus == TotalInternalReflection) { 
              run->AddTotalInternalReflection();
            }
            else if (theStatus == LambertianReflection) {
              run->AddLambertianReflection();
            }
            else if (theStatus == LobeReflection) {
              run->AddLobeReflection();
            }
            else if (theStatus == SpikeReflection) {
              run->AddSpikeReflection();
            }
            else if (theStatus == BackScattering) {
              run->AddBackScattering();
            }
            else if (theStatus == Absorption) {
              run->AddAbsorption();
            }
            else if (theStatus == Detection) {
              run->AddDetection();
            }
            else if (theStatus == NotAtBoundary) {
              run->AddNotAtBoundary();
            }
            else if (theStatus == SameMaterial) {
              run->AddSameMaterial();
            }
            else if (theStatus == StepTooSmall) {
              run->AddStepTooSmall();
            }
            else if (theStatus == NoRINDEX) {
              run->AddNoRINDEX();
            }
            else if (theStatus == PolishedLumirrorAirReflection) {
              run->AddPolishedLumirrorAirReflection();
            }
            else if (theStatus == PolishedLumirrorGlueReflection) {
              run->AddPolishedLumirrorGlueReflection();
            }
            else if (theStatus == PolishedAirReflection) {
              run->AddPolishedAirReflection();
            }
            else if (theStatus == PolishedTeflonAirReflection) {
              run->AddPolishedTeflonAirReflection();
            }
            else if (theStatus == PolishedTiOAirReflection) {
              run->AddPolishedTiOAirReflection();
            }
            else if (theStatus == PolishedTyvekAirReflection) {
              run->AddPolishedTyvekAirReflection();
            }
            else if (theStatus == PolishedVM2000AirReflection) {
              run->AddPolishedVM2000AirReflection();
            }
            else if (theStatus == PolishedVM2000GlueReflection) {
              run->AddPolishedVM2000AirReflection();
            }
            else if (theStatus == EtchedLumirrorAirReflection) {
              run->AddEtchedLumirrorAirReflection();
            }
            else if (theStatus == EtchedLumirrorGlueReflection) {
              run->AddEtchedLumirrorGlueReflection();
            }
            else if (theStatus == EtchedAirReflection) {
              run->AddEtchedAirReflection();
            }
            else if (theStatus == EtchedTeflonAirReflection) {
              run->AddEtchedTeflonAirReflection();
            }
            else if (theStatus == EtchedTiOAirReflection) {
              run->AddEtchedTiOAirReflection();
            }
            else if (theStatus == EtchedTyvekAirReflection) {
              run->AddEtchedTyvekAirReflection();
            }
            else if (theStatus == EtchedVM2000AirReflection) {
              run->AddEtchedVM2000AirReflection();
            }
            else if (theStatus == EtchedVM2000GlueReflection) {
              run->AddEtchedVM2000AirReflection();
            }
            else if (theStatus == GroundLumirrorAirReflection) {
              run->AddGroundLumirrorAirReflection();
            }
            else if (theStatus == GroundLumirrorGlueReflection) {
              run->AddGroundLumirrorGlueReflection();
            }
            else if (theStatus == GroundAirReflection) {
              run->AddGroundAirReflection();
            }
            else if (theStatus == GroundTeflonAirReflection) {
              run->AddGroundTeflonAirReflection();
            }
            else if (theStatus == GroundTiOAirReflection) {
              run->AddGroundTiOAirReflection();
            }
            else if (theStatus == GroundTyvekAirReflection) {
              run->AddGroundTyvekAirReflection();
            }
            else if (theStatus == GroundVM2000AirReflection) {
              run->AddGroundVM2000AirReflection();
            }
            else if (theStatus == GroundVM2000GlueReflection) {
              run->AddGroundVM2000AirReflection();
            }
            else if (theStatus == Dichroic) {
              run->AddDichroic();
            }
            
            else {
              G4cout << "theStatus: " << theStatus 
                     << " was none of the above." << G4endl;
            }

          }
        }
      }
    }
  }

  else { // particle != opticalphoton
    // print how many Cerenkov and scint photons produced this step
    // this demonstrates use of GetNumPhotons()
    auto proc_man = track->GetDynamicParticle()->GetParticleDefinition()
                         ->GetProcessManager();
    G4int n_proc = proc_man->GetPostStepProcessVector()->entries();
    G4ProcessVector* proc_vec = proc_man->GetPostStepProcessVector(typeDoIt);

    G4int n_scint = 0;
    G4int n_cer   = 0;
    for (G4int i = 0; i < n_proc; ++i) {
      if ((*proc_vec)[i]->GetProcessName().compare("Cerenkov") == 0) {
        auto cer = (G4Cerenkov*)(*proc_vec)[i];
        n_cer = cer->GetNumPhotons();
      }
      else if ((*proc_vec)[i]->GetProcessName().compare("Scintillation") == 0) {
        auto scint = (G4Scintillation*)(*proc_vec)[i];
        n_scint = scint->GetNumPhotons();
      }
    }
    if (fVerbose > 0) {
      if (n_cer > 0 || n_scint > 0) {
        G4cout << "In this step, " << n_cer
               << " Cerenkov and " << n_scint
               << " scintillation photons were produced." << G4endl;
      }
    }

    // loop over secondaries, create statistics
    const std::vector<const G4Track*>* secondaries =
                                step->GetSecondaryInCurrentStep();

    for (auto sec : *secondaries) {
      if (sec->GetDynamicParticle()->GetParticleDefinition() == opticalphoton){
        if (sec->GetCreatorProcess()->GetProcessName().compare("Cerenkov")==0){
          G4double en = sec->GetKineticEnergy();
          run->AddCerenkovEnergy(en);
          run->AddCerenkov();
          G4AnalysisManager::Instance()->FillH1(1, en/eV);
        }
        else if (sec->GetCreatorProcess()
                    ->GetProcessName().compare("Scintillation") == 0) {
          G4double en = sec->GetKineticEnergy();
          run->AddScintillationEnergy(en);
          run->AddScintillation();
          G4AnalysisManager::Instance()->FillH1(2, en/eV);

          G4double time = sec->GetGlobalTime();
          analysisMan->FillH1(13, time/ns);
        }
      }
    }
  } 

  return;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
