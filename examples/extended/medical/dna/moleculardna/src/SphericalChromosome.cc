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
/// file: SphericalChromosome.cc
/// brief: Implementation of virt chromosome class for Spherical chromosomes

#include <utility>

#include "SphericalChromosome.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4String SphericalChromosome::fShape = "sphere";

SphericalChromosome::SphericalChromosome(const G4String& name,
                                         G4ThreeVector pos,
                                         const G4double& radius)
  : VirtualChromosome(name)
  , fCenter(std::move(pos))
  , fRadius(radius)
  , fRotation(G4RotationMatrix())
{
  fInverseRotation = fRotation.inverse();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SphericalChromosome::SphericalChromosome(const G4String& name,
                                         G4ThreeVector pos,
                                         const G4double& radius,
                                         G4RotationMatrix rot)
  : VirtualChromosome(name)
  , fCenter(std::move(pos))
  , fRadius(radius)
  , fRotation(std::move(rot))
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SphericalChromosome::~SphericalChromosome() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool SphericalChromosome::PointInChromosome(G4ThreeVector const& pos)
{
  G4ThreeVector rpos = pos - fCenter;
  rpos               = fInverseRotation(rpos);
  G4bool radius_ok   = (rpos.mag2() < fRadius * fRadius);
  return radius_ok;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector SphericalChromosome::RandomPointInChromosome()
{
  G4ThreeVector point = fRadius * G4UniformRand() * G4RandomDirection();
  return fRotation(point) + fCenter;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SphericalChromosome::Print()
{
  G4cout << "type: " << fShape << G4endl;
  G4cout << "radius: " << fRadius << G4endl;
  G4cout << "center: " << fCenter << G4endl;
  G4cout << "rotation: " << fRotation.getPhi() << " " << fRotation.getTheta()
  << " " << fRotation.getPhi() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
