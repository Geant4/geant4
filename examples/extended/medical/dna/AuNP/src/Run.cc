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
/// \file Run.cc
/// \brief Implementation of the Run class

#include "Run.hh"

#include "DetectorConstruction.hh"
#include "HistoManager.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4Material.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::RecordEvent(const G4Event * /*aEvent*/) {
  numberOfEvent++; // This is an original line.
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4THitsMap<G4double> *Run::GetHitsMap(const G4String &detName, const G4String &colName) const {
  const G4String fullName = detName + "/" + colName;
  return GetHitsMap(fullName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4THitsMap<G4double> *Run::GetHitsMap(const G4String &fullName) const {
  if (fCollName == fullName) {
    return fRunMap;
  }
  return nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::DumpAllScorer() const {
  // - Number of HitsMap in this RUN.
  // - GetHitsMap and dump values.

  G4THitsMap<G4double> *RunMap = GetHitsMap();
  if (RunMap) {
    G4cout << " PrimitiveScorer RUN " << RunMap->GetSDname() << "," << RunMap->GetName() << G4endl;
    G4cout << " Number of entries " << RunMap->entries() << G4endl;
    std::map<G4int, G4double *>::iterator itr = RunMap->GetMap()->begin();
    for (; itr != RunMap->GetMap()->end(); itr++) {
      G4cout << "  copy no.: " << itr->first << "  Run Value : " << *(itr->second) << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetPrimary(G4ParticleDefinition *particle, G4double energy) {
  fParticle = particle;
  fEkin = energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run *run) {
  const auto localRun = static_cast<const Run *>(run);

  fParticle = localRun->fParticle;
  fEkin = localRun->fEkin;

  G4Run::Merge(run);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::EndOfRun() {
  const std::ios::fmtflags mode = G4cout.flags();
  G4cout.setf(std::ios::fixed, std::ios::floatfield);
  const G4int prec = G4cout.precision(2);

  G4cout.setf(mode, std::ios::floatfield);
  G4cout.precision(prec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......