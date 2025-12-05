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
/// \file GB03Run.cc
/// \brief Implementation of the GB03Run class

#include "GB03Run.hh"

#include <algorithm>
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ParticlesCount::Count(const G4String& partName, G4double weight, G4double E)
{
  // G4cout << "Count " << partName << " " << weight << G4endl;

  const auto& [it, inserted] = fPartCntr.try_emplace(
    partName, std::tuple<G4double, G4int, G4double, G4double>(weight, 1, E, E));

  // G4cout << "Count " << inserted << " " << G4endl;

  // return;

  if (!inserted) {
    // Exists already
    std::get<0>(it->second) += weight;
    std::get<1>(it->second) += 1;
    std::get<2>(it->second) = std::min(std::get<2>(it->second), E);
    std::get<3>(it->second) = std::max(std::get<3>(it->second), E);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ParticlesCount::Print() const
{
  for (const auto& it : fPartCntr) {
    G4String partName = it.first;
    auto [weight, count, Emin, Emax] = it.second;
    G4cout << " " << std::setw(20) << partName << " " << std::setw(7) << count << " "
           << std::setw(12) << weight << " Emin " << std::setw(12) << Emin << " Emax "
           << std::setw(12) << Emax << " (MeV) " << G4endl;
  }
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ParticlesCount::Merge(const ParticlesCount& worker)
{
  for (const auto& part : worker.fPartCntr) {
    const auto& [it, inserted] = fPartCntr.emplace(part);
    if (!inserted) {
      // already in master counters
      std::get<0>(it->second) += std::get<0>(part.second);
      std::get<1>(it->second) += std::get<1>(part.second);
      std::get<2>(it->second) = std::min(std::get<2>(it->second), std::get<2>(part.second));
      std::get<3>(it->second) = std::max(std::get<3>(it->second), std::get<3>(part.second));
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Destructor
GB03Run::~GB03Run()
{
  //    clear all data members.
  fPartCounter.Reset();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB03Run::CountParticle(const G4String& n, G4double w, G4double e)
{
  fPartCounter.Count(n, w, e);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB03Run::Merge(const G4Run* aRun)
{
  const GB03Run* localRun = static_cast<const GB03Run*>(aRun);
  //=======================================================
  // Merge Partice Counters
  //=======================================================
  fPartCounter.Merge(localRun->fPartCounter);

  G4Run::Merge(aRun);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
