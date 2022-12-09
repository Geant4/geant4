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
/// \file src/Run.cc
/// \brief Implementation of the Run class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"

#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Positron.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(DetectorConstruction* det)
  : G4Run(),
    fDetector(det),
    fParticle(nullptr),
    fEkin(0.),
    fChargedStep(0),
    fNeutralStep(0),
    fN_gamma(0),
    fN_elec(0),
    fN_pos(0)
{
  // initialize cumulative quantities
  //
  for (G4int k = 0; k < kMaxAbsor; k++) {
    fSumEAbs[k] = fSum2EAbs[k] = fSumLAbs[k] = fSum2LAbs[k] = 0.;
    fEnergyDeposit[k].clear();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetPrimary(G4ParticleDefinition* particle, G4double energy)
{
  fParticle = particle;
  fEkin = energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::FillPerEvent(G4int kAbs, G4double EAbs, G4double LAbs)
{
  // accumulate statistic with restriction
  //
  fEnergyDeposit[kAbs].push_back(EAbs);
  fSumEAbs[kAbs] += EAbs;
  fSum2EAbs[kAbs] += EAbs * EAbs;
  fSumLAbs[kAbs] += LAbs;
  fSum2LAbs[kAbs] += LAbs * LAbs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::AddChargedStep() { fChargedStep += 1.0; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::AddNeutralStep() { fNeutralStep += 1.0; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::AddSecondaryTrack(const G4Track* track)
{
  const G4ParticleDefinition* d = track->GetDefinition();
  if (d == G4Gamma::Gamma()) {
    ++fN_gamma;
  }
  else if (d == G4Electron::Electron()) {
    ++fN_elec;
  }
  else if (d == G4Positron::Positron()) {
    ++fN_pos;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);

  // pass information about primary particle
  fParticle = localRun->fParticle;
  fEkin = localRun->fEkin;

  // accumulate sums
  //
  for (G4int k = 0; k < kMaxAbsor; k++) {
    fSumEAbs[k] += localRun->fSumEAbs[k];
    fSum2EAbs[k] += localRun->fSum2EAbs[k];
    fSumLAbs[k] += localRun->fSumLAbs[k];
    fSum2LAbs[k] += localRun->fSum2LAbs[k];
  }

  fChargedStep += localRun->fChargedStep;
  fNeutralStep += localRun->fNeutralStep;

  fN_gamma += localRun->fN_gamma;
  fN_elec += localRun->fN_elec;
  fN_pos += localRun->fN_pos;

  G4Run::Merge(run);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun()
{
  G4int nEvt = numberOfEvent;
  G4double norm = G4double(nEvt);
  if (norm > 0) norm = 1. / norm;
  G4double qnorm = std::sqrt(norm);

  fChargedStep *= norm;
  fNeutralStep *= norm;

  // compute and print statistic
  //
  G4double beamEnergy = fEkin;
  G4double sqbeam = std::sqrt(beamEnergy / GeV);

  G4double MeanEAbs, MeanEAbs2, rmsEAbs, resolution, rmsres;
  G4double MeanLAbs, MeanLAbs2, rmsLAbs;

  std::ios::fmtflags mode = G4cout.flags();
  G4int prec = G4cout.precision(2);
  G4cout << "\n------------------------------------------------------------\n";
  G4cout << std::setw(14) << "material" << std::setw(17) << "Edep       RMS" << std::setw(33)
         << "sqrt(E0(GeV))*rmsE/Emean" << std::setw(23) << "total tracklen \n \n";

  for (G4int k = 1; k <= fDetector->GetNbOfAbsor(); k++) {
    MeanEAbs = fSumEAbs[k] * norm;
    MeanEAbs2 = fSum2EAbs[k] * norm;
    rmsEAbs = std::sqrt(std::abs(MeanEAbs2 - MeanEAbs * MeanEAbs));

    resolution = 100. * sqbeam * rmsEAbs / MeanEAbs;
    rmsres = resolution * qnorm;

    // Save mean and RMS
    fSumEAbs[k] = MeanEAbs;
    fSum2EAbs[k] = rmsEAbs;

    MeanLAbs = fSumLAbs[k] * norm;
    MeanLAbs2 = fSum2LAbs[k] * norm;
    rmsLAbs = std::sqrt(std::abs(MeanLAbs2 - MeanLAbs * MeanLAbs));

    // print
    //
    G4cout << std::setw(14) << fDetector->GetAbsorMaterial(k)->GetName() << ": "
           << std::setprecision(5) << std::setw(6) << G4BestUnit(MeanEAbs, "Energy") << " :  "
           << std::setprecision(4) << std::setw(5) << G4BestUnit(rmsEAbs, "Energy") << std::setw(10)
           << resolution << " +- " << std::setw(5) << rmsres << " %" << std::setprecision(3)
           << std::setw(10) << G4BestUnit(MeanLAbs, "Length") << " +- " << std::setw(4)
           << G4BestUnit(rmsLAbs, "Length") << G4endl;
  }
  G4cout << "\n------------------------------------------------------------\n";

  G4cout << " Beam particle " << fParticle->GetParticleName()
         << "  E = " << G4BestUnit(beamEnergy, "Energy") << G4endl;
  G4cout << " Mean number of gamma       " << (G4double)fN_gamma * norm << G4endl;
  G4cout << " Mean number of e-          " << (G4double)fN_elec * norm << G4endl;
  G4cout << " Mean number of e+          " << (G4double)fN_pos * norm << G4endl;
  G4cout << std::setprecision(6) << " Mean number of charged steps  " << fChargedStep << G4endl;
  G4cout << " Mean number of neutral steps  " << fNeutralStep << G4endl;
  G4cout << "------------------------------------------------------------\n" << G4endl;

  G4cout.setf(mode, std::ios::floatfield);
  G4cout.precision(prec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
