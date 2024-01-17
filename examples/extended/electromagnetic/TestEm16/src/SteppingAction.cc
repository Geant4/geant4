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
/// \file electromagnetic/TestEm16/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "HistoManager.hh"
#include "Run.hh"

#include "G4EmProcessSubType.hh"
#include "G4ParticleTypes.hh"
#include "G4RunManager.hh"
#include "G4SteppingManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4VProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction() : G4UserSteppingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  const G4VProcess* process = aStep->GetPostStepPoint()->GetProcessDefinedStep();
  if (process == nullptr) return;

  G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  G4Track* trk = aStep->GetTrack();

  thread_local G4int iCalled = 0;
  const G4int nprint = 10;  // set to 10 to get debug print for first 10 calls
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  if (fSynchrotronRadiation == process->GetProcessSubType()) {
    const G4StepPoint* PrePoint = aStep->GetPreStepPoint();
    const G4TrackVector* secondary = fpSteppingManager->GetSecondary();
    size_t lp = (*secondary).size();
    if (lp) {
      G4double Egamma = (*secondary)[lp - 1]->GetTotalEnergy();
      run->f_n_gam_sync++;
      run->f_e_gam_sync += Egamma;
      run->f_e_gam_sync2 += Egamma * Egamma;
      if (Egamma > run->f_e_gam_sync_max) run->f_e_gam_sync_max = Egamma;
      run->f_lam_gam_sync += aStep->GetStepLength();
      if (iCalled < nprint) {
        G4double Eelec = PrePoint->GetTotalEnergy();
        G4ThreeVector Pelec = PrePoint->GetMomentum();
        G4ThreeVector PGamma = (*secondary)[lp - 1]->GetMomentum();
        G4bool IsGamma = ((*secondary)[lp - 1]->GetDefinition() == G4Gamma::Gamma());
        const G4int wid = 8;
        if (analysisManager->GetVerboseLevel() > 1 && iCalled < nprint)
          G4cout << __FILE__ << " line " << __LINE__ << " " << __FUNCTION__
                 << " processName=" << process->GetProcessName()
                 << " Step Length=" << std::setw(wid)
                 << G4BestUnit(aStep->GetStepLength(), "Length") << " Eelec=" << std::setw(wid)
                 << G4BestUnit(Eelec, "Energy") << " Pelec=" << G4BestUnit(Pelec, "Energy")
                 << " IsGamma=" << IsGamma << " Egamma=" << std::setw(8)
                 << G4BestUnit(Egamma, "Energy") << " PGamma=" << std::setw(8)
                 << G4BestUnit(PGamma, "Energy") << " #secondaries lp=" << lp << '\n';
      }
      analysisManager->FillH1(1, Egamma);
      analysisManager->FillH1(2, Egamma,
                              Egamma / keV);  // power spectrum : gamma weighted with its energy
      analysisManager->FillH1(3, aStep->GetStepLength());
    }
  }

  thread_local G4ThreeVector IncomingPhotonDirection;
  if (fGammaReflection == process->GetProcessSubType()) {
    const G4TrackVector* secondary = fpSteppingManager->GetSecondary();
    size_t lp = (*secondary).size();
    if (lp) {
      G4double Egamma = (*secondary)[lp - 1]->GetTotalEnergy();
      ++run->f_n_Xray_Refl;
      analysisManager->FillH1(1, Egamma, 1. / run->GetNumberOfEventToBeProcessed());
      if (analysisManager->GetVerboseLevel() > 1 && run->f_n_Xray_Refl < nprint ) 
        G4cout << __FILE__ << " line " << __LINE__ << " " << __FUNCTION__
               << " iCalled=" << std::setw(3) << iCalled
               << " eventID=" << std::setw(3) << eventID
               << " f_n_Xray_Refl=" << run->f_n_Xray_Refl
               << " ProcessName=" << process->GetProcessName()
               << " direction=" << trk->GetMomentumDirection() << " Egamma=" << Egamma / keV
               << " keV" << G4endl;
      IncomingPhotonDirection = trk->GetMomentumDirection();
    }
  }

  const G4VProcess* creator = trk->GetCreatorProcess();
  G4int CreatorProcessSubType=0;
  if (creator != nullptr) CreatorProcessSubType = creator->GetProcessSubType();
  if (CreatorProcessSubType == fGammaReflection) {  // transportation after XrayReflection
    if (IncomingPhotonDirection != G4ThreeVector(0, 0, 0)) {
      G4double cos_angle =
        IncomingPhotonDirection * trk->GetMomentumDirection();  // incoming * outgoing direction
      G4double IncidentAngle = std::acos(cos_angle) / 2.;
      analysisManager->FillH1(4, IncidentAngle, 1. / run->GetNumberOfEventToBeProcessed());
      if (analysisManager->GetVerboseLevel() > 1 && iCalled < nprint)
        G4cout << __FILE__ << " line " << __LINE__ << " " << __FUNCTION__
               << " iCalled=" << std::setw(3) << iCalled
               << " eventID=" << std::setw(3) << eventID
               << " CreatorProcname=" << creator->GetProcessName()
               << " ProcessName=" << process->GetProcessName()
               << " IncomingPhotonDirection=" << IncomingPhotonDirection
               << "  outgoing PhotonDirection=" << trk->GetMomentumDirection()
               << " IncidentAngle=" << IncidentAngle
               << " NumberOfEventToBeProcessed=" << run->GetNumberOfEventToBeProcessed()
               << " RunID=" << run->GetRunID() << G4endl;
    }
  }

  ++iCalled;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
