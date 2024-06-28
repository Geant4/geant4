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

#include "BoundedBrownianAction.hh"

#include "G4DNABoundingBox.hh"
#include "G4Molecule.hh"
#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BoundedBrownianAction::BoundedBrownianAction() : G4VUserBrownianAction() {}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double BoundedBrownianAction::GetDistanceToBoundary(const G4Track& track)
{
  if (!fpBoundingBox->contains(track.GetPosition())) {
    G4ExceptionDescription errMsg;
    errMsg << "Point is out of box : " << *fpBoundingBox
           << " of particle : " << GetIT(track)->GetName() << "(" << track.GetTrackID()
           << ") : " << track.GetPosition();
    G4Exception(
      "BoundedBrownianAction::GetDistanceToBoundary"
      "BoundedBrownianAction",
      "BoundedBrownianAction", FatalErrorInArgument, errMsg);
  }

  auto dx = std::min(track.GetPosition().getX() - fpBoundingBox->Getxlo(),
                     fpBoundingBox->Getxhi() - track.GetPosition().getX());
  auto dy = std::min(track.GetPosition().getY() - fpBoundingBox->Getylo(),
                     fpBoundingBox->Getyhi() - track.GetPosition().getY());
  auto dz = std::min(track.GetPosition().getZ() - fpBoundingBox->Getzlo(),
                     fpBoundingBox->Getzhi() - track.GetPosition().getZ());
  return std::min({dx, dy, dz});
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BoundedBrownianAction::Transport(G4ThreeVector& nextposition, G4Track*)
{
  BouncingAction(nextposition);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BoundedBrownianAction::BouncingAction(G4ThreeVector& nextposition) const
{
  // This algorithm is implemented based on Karamitros, Mathieu et al.2020,arXiv:2006.14225 (2020)
  //  https://doi.org/10.48550/arXiv.2006.14225
  G4double RxM = fpBoundingBox->Getxhi();
  G4double RyM = fpBoundingBox->Getyhi();
  G4double RzM = fpBoundingBox->Getzhi();

  G4double Rxm = fpBoundingBox->Getxlo();
  G4double Rym = fpBoundingBox->Getylo();
  G4double Rzm = fpBoundingBox->Getzlo();
  G4double Lx = RxM - Rxm;
  G4double Ly = RyM - Rym;
  G4double Lz = RzM - Rzm;
  G4double trunx = trunc(std::abs(nextposition.getX() - Rxm) / Lx);
  G4double truny = trunc(std::abs(nextposition.getY() - Rym) / Ly);
  G4double trunz = trunc(std::abs(nextposition.getZ() - Rzm) / Lz);
  G4double hx = std::fmod(trunx, 2);
  G4double hy = std::fmod(truny, 2);
  G4double hz = std::fmod(trunz, 2);
  G4double x = Rxm + hx * Lx + (1 - 2 * hx) * std::abs(std::fmod(nextposition.getX() - Rxm, Lx));
  G4double y = Rym + hy * Ly + (1 - 2 * hy) * std::abs(std::fmod(nextposition.getY() - Rym, Ly));
  G4double z = Rzm + hz * Lz + (1 - 2 * hz) * std::abs(std::fmod(nextposition.getZ() - Rzm, Lz));
  nextposition.set(x, y, z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
