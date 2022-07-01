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
/// \file optical/OpNovice2/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"

#include "HistoManager.hh"
#include "PrimaryGeneratorAction.hh"
#include "Run.hh"

#include "G4Run.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(PrimaryGeneratorAction* prim)
  : G4UserRunAction()
  , fRun(nullptr)
  , fHistoManager(nullptr)
  , fPrimary(prim)
{
  fHistoManager = new HistoManager();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction() { delete fHistoManager; }

G4Run* RunAction::GenerateRun()
{
  fRun = new Run();
  return fRun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  if(fPrimary)
  {
    G4ParticleDefinition* particle =
      fPrimary->GetParticleGun()->GetParticleDefinition();
    G4double energy       = fPrimary->GetParticleGun()->GetParticleEnergy();
    G4bool polarized      = fPrimary->GetPolarized();
    G4double polarization = fPrimary->GetPolarization();
    fRun->SetPrimary(particle, energy, polarized, polarization);
  }

  // histograms
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if(analysisManager->IsActive())
  {
    analysisManager->OpenFile();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
  if(isMaster)
    fRun->EndOfRun();

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  G4cout << G4endl << " Histogram statistics for the ";
  if(isMaster)
  {
    G4cout << "entire run:" << G4endl << G4endl;
  }
  else
  {
    G4cout << "local thread:" << G4endl << G4endl;
  }

  G4int id = analysisManager->GetH1Id("Cerenkov spectrum");
  if(analysisManager->GetH1Activation(id))
  {
    G4cout << " Cerenkov spectrum: mean = "
           << analysisManager->GetH1(id)->mean()
           << " eV; rms = " << analysisManager->GetH1(id)->rms() << " eV."
           << G4endl;
  }
  id = analysisManager->GetH1Id("Scintillation spectrum");
  if(analysisManager->GetH1Activation(id))
  {
    G4cout << " Scintillation spectrum: mean = "
           << analysisManager->GetH1(id)->mean()
           << " eV; rms = " << analysisManager->GetH1(id)->rms() << " eV."
           << G4endl;
  }
  id = analysisManager->GetH1Id("Scintillation time");
  if(analysisManager->GetH1Activation(id))
  {
    G4cout << " Scintillation time: mean = "
           << analysisManager->GetH1(id)->mean()
           << " ns; rms = " << analysisManager->GetH1(id)->rms() << " ns."
           << G4endl;
  }
  id = analysisManager->GetH1Id("WLS abs");
  if(analysisManager->GetH1Activation(id))
  {
    G4cout << " WLS absorption spectrum: mean = "
           << analysisManager->GetH1(id)->mean()
           << " eV; rms = " << analysisManager->GetH1(id)->rms() << " eV."
           << G4endl;
  }
  id = analysisManager->GetH1Id("WLS em");
  if(analysisManager->GetH1Activation(id))
  {
    G4cout << " WLS emission spectrum: mean = "
           << analysisManager->GetH1(id)->mean()
           << " eV; rms = " << analysisManager->GetH1(id)->rms() << " eV."
           << G4endl;
  }
  id = analysisManager->GetH1Id("WLS time");
  if(analysisManager->GetH1Activation(id))
  {
    G4cout << " WLS emission time: mean = "
           << analysisManager->GetH1(id)->mean()
           << " ns; rms = " << analysisManager->GetH1(id)->rms() << " ns."
           << G4endl;
  }
  id = analysisManager->GetH1Id("WLS2 abs");
  if(analysisManager->GetH1Activation(id))
  {
    G4cout << " WLS emission time: mean = "
           << analysisManager->GetH1(id)->mean()
           << " ns; rms = " << analysisManager->GetH1(id)->rms() << " ns."
           << G4endl;
  }
  id = analysisManager->GetH1Id("WLS2 em");
  if(analysisManager->GetH1Activation(id))
  {
    G4cout << " WLS2 emission spectrum: mean = "
           << analysisManager->GetH1(id)->mean()
           << " eV; rms = " << analysisManager->GetH1(id)->rms() << " eV."
           << G4endl;
  }
  id = analysisManager->GetH1Id("WLS2 time");
  if(analysisManager->GetH1Activation(id))
  {
    G4cout << " WLS2 emission time: mean = "
           << analysisManager->GetH1(id)->mean()
           << " ns; rms = " << analysisManager->GetH1(id)->rms() << " ns."
           << G4endl;
  }
  id = analysisManager->GetH1Id("x_backward");
  if(analysisManager->GetH1Activation(id))
  {
    G4cout << " X momentum dir of backward-going photons: mean = "
           << analysisManager->GetH1(id)->mean()
           << "; rms = " << analysisManager->GetH1(id)->rms() << G4endl;
  }
  id = analysisManager->GetH1Id("y_backward");
  if(analysisManager->GetH1Activation(id))
  {
    G4cout << " Y momentum dir of backward-going photons: mean = "
           << analysisManager->GetH1(id)->mean()
           << "; rms = " << analysisManager->GetH1(id)->rms() << G4endl;
  }
  id = analysisManager->GetH1Id("z_backward");
  if(analysisManager->GetH1Activation(id))
  {
    G4cout << " Z momentum dir of backward-going photons: mean = "
           << analysisManager->GetH1(id)->mean()
           << "; rms = " << analysisManager->GetH1(id)->rms() << G4endl;
  }
  id = analysisManager->GetH1Id("x_forward");
  if(analysisManager->GetH1Activation(id))
  {
    G4cout << " X momentum dir of forward-going photons: mean = "
           << analysisManager->GetH1(id)->mean()
           << "; rms = " << analysisManager->GetH1(id)->rms() << G4endl;
  }
  id = analysisManager->GetH1Id("y_forward");
  if(analysisManager->GetH1Activation(id))
  {
    G4cout << " Y momentum dir of forward-going photons: mean = "
           << analysisManager->GetH1(id)->mean()
           << "; rms = " << analysisManager->GetH1(id)->rms() << G4endl;
  }
  id = analysisManager->GetH1Id("z_forward");
  if(analysisManager->GetH1Activation(id))
  {
    G4cout << " Z momentum dir of forward-going photons: mean = "
           << analysisManager->GetH1(id)->mean()
           << "; rms = " << analysisManager->GetH1(id)->rms() << G4endl;
  }
  id = analysisManager->GetH1Id("x_fresnel");
  if(analysisManager->GetH1Activation(id))
  {
    G4cout << " X momentum dir of Fresnel-refracted photons: mean = "
           << analysisManager->GetH1(id)->mean()
           << "; rms = " << analysisManager->GetH1(id)->rms() << G4endl;
  }
  id = analysisManager->GetH1Id("y_fresnel");
  if(analysisManager->GetH1Activation(id))
  {
    G4cout << " Y momentum dir of Fresnel-refracted photons: mean = "
           << analysisManager->GetH1(id)->mean()
           << "; rms = " << analysisManager->GetH1(id)->rms() << G4endl;
  }
  id = analysisManager->GetH1Id("z_fresnel");
  if(analysisManager->GetH1Activation(id))
  {
    G4cout << " Z momentum dir of Fresnel-refracted photons: mean = "
           << analysisManager->GetH1(id)->mean()
           << "; rms = " << analysisManager->GetH1(id)->rms() << G4endl;
  }
  id = analysisManager->GetH1Id("Transmitted");
  if(analysisManager->GetH1Activation(id))
  {
    G4cout << " Angle of transmitted photons: mean = "
           << analysisManager->GetH1(id)->mean()
           << "; rms = " << analysisManager->GetH1(id)->rms() << G4endl;
  }
  id = analysisManager->GetH1Id("Fresnel reflection");
  if(analysisManager->GetH1Activation(id))
  {
    G4cout << " Angle of Fresnel-reflected photons: mean = "
           << analysisManager->GetH1(id)->mean()
           << "; rms = " << analysisManager->GetH1(id)->rms() << G4endl;
  }
  id = analysisManager->GetH1Id("Total internal reflection");
  if(analysisManager->GetH1Activation(id))
  {
    G4cout << " Angle of total internal reflected photons: mean = "
           << analysisManager->GetH1(id)->mean()
           << "; rms = " << analysisManager->GetH1(id)->rms() << G4endl;
  }

  id = analysisManager->GetH1Id("Fresnel refraction");
  if(analysisManager->GetH1Activation(id))
  {
    G4cout << " Angle of Fresnel-refracted photons: mean = "
           << analysisManager->GetH1(id)->mean()
           << "; rms = " << analysisManager->GetH1(id)->rms() << G4endl;
  }

  id = analysisManager->GetH1Id("Absorption");
  if(analysisManager->GetH1Activation(id))
  {
    G4cout << " Angle of absorbed photons: mean = "
           << analysisManager->GetH1(id)->mean()
           << "; rms = " << analysisManager->GetH1(id)->rms() << G4endl;
  }

  G4cout << G4endl;

  if(analysisManager->IsActive())
  {
    analysisManager->Write();
    analysisManager->CloseFile();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
