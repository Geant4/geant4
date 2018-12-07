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
#include "Run.hh"

#include "G4RunManager.hh"
#include "G4SteppingManager.hh"
#include "G4VProcess.hh"
#include "G4ParticleTypes.hh"
#include "G4UnitsTable.hh"
#include "HistoManager.hh"
#include "G4EmProcessSubType.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction()
:G4UserSteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  Run* run = static_cast<Run*>(
             G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  
 const G4VProcess* process = aStep->GetPostStepPoint()->GetProcessDefinedStep();
  if (process == 0) { return; }

  static G4int iCalled=0;
  const  G4int nprint=0; // set to 10 to get debug print for first 10 calls

  if (fSynchrotronRadiation == process->GetProcessSubType())
    {
      ++iCalled;
      const G4StepPoint* PrePoint = aStep->GetPreStepPoint();
      const G4TrackVector* secondary = fpSteppingManager->GetSecondary();

      //(*secondary)[lp-1]  points to the last photon generated
      //   
      size_t lp=(*secondary).size();
      if (lp)
        {
          G4double  Egamma =  (*secondary)[lp-1]->GetTotalEnergy();
          run->f_n_gam_sync++;
          run->f_e_gam_sync += Egamma;
          run->f_e_gam_sync2 += Egamma*Egamma;
          if (Egamma > run->f_e_gam_sync_max) 
            { 
              run->f_e_gam_sync_max = Egamma;
            }
          run->f_lam_gam_sync += aStep->GetStepLength();
          if (iCalled<nprint)
            {
              G4double      Eelec  = PrePoint->GetTotalEnergy();
              G4ThreeVector Pelec  = PrePoint->GetMomentum();
              G4ThreeVector PGamma = (*secondary)[lp-1]->GetMomentum();
              G4bool IsGamma = 
                ((*secondary)[lp-1]->GetDefinition() == G4Gamma::Gamma());
            
              G4cout << "UserSteppingAction processName="
                     << process->GetProcessName()
                     << " Step Length=" << std::setw(6) 
                     << G4BestUnit(aStep->GetStepLength(),"Length")
                     << " Eelec=" << G4BestUnit(Eelec,"Energy")
                     << " Pelec=" << G4BestUnit(Pelec,"Energy")
                     << " IsGamma=" << IsGamma
                     << " Egamma=" << G4BestUnit(Egamma,"Energy")
                     << " PGamma=" << G4BestUnit(PGamma,"Energy")
                     << " #secondaries lp=" << lp
                     << '\n';
            }

          G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
          analysisManager->FillH1(1,Egamma);

          // power spectrum : gamma weighted with its energy
          analysisManager->FillH1(2,Egamma,Egamma/keV); 

          analysisManager->FillH1(3,aStep->GetStepLength());
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
