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
#include "SiPMHit.hh"

#include <cstdlib>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SiPMHit::SiPMHit(G4String aVolumeName, G4int aCopyNumSensor,
                 G4int aCopyNumCell)
    : fVolumeName(aVolumeName), fCopyNumCell(aCopyNumCell),
      fCopyNumSensor(aCopyNumSensor) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SiPMHit::Digitise(const G4double aTimeWindow,
                       const G4double aToaThreshold) {

  // process energy deposits
  if (fEdep.empty()) {
    fIsValidHit = false;
    return;
  }

  std::sort(fEdep.begin(), fEdep.end(),
            [](const std::pair<G4double, G4double> &left,
               const std::pair<G4double, G4double> &right) {
              return left.second < right.second; // second = time
            });

  G4double firstHitTime = fEdep[0].second;
  fEdepDigi = 0;
  for (size_t i = 0; i < fEdep.size(); ++i) {
    if (aTimeWindow < 0 || fEdep[i].second < firstHitTime + aTimeWindow)
      fEdepDigi += fEdep[i].first;
#ifdef DEBUG
    else
      std::cout << "Rejecting hit (" << fEdep[i].first << " keV)  at time "
                << fEdep[i].second << "(start: " << firstHitTime << ")"
                << std::endl;
#endif
    if ((fTimeOfArrival == -1) && (fEdepDigi > aToaThreshold))
      fTimeOfArrival = fEdep[i].second;
  }
  fIsValidHit = (fEdepDigi > 0);

  // non ionizing part, does not contribute to TOAs
  if (fEdepNonIonizing.empty())
    return;

  std::sort(fEdepNonIonizing.begin(), fEdepNonIonizing.end(),
            [](const std::pair<G4double, G4double> &left,
               const std::pair<G4double, G4double> &right) {
              return left.second < right.second; // second = time
            });

  fEdepNonIonizingDigi = 0;
  for (size_t i = 0; i < fEdepNonIonizing.size(); ++i) {
    if (aTimeWindow == -1 ||
        fEdepNonIonizing[i].second < firstHitTime + aTimeWindow)
      fEdepNonIonizingDigi += fEdepNonIonizing[i].first;
  }
}
