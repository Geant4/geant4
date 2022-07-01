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
/// file: EllipticalChromosome.cc
/// brief: Implementation of virt chromosome class for rod chromosomes

#include "EllipticalChromosome.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4PhysicalConstants.hh"

const G4String EllipticalChromosome::fShape = "ellipse";

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EllipticalChromosome::EllipticalChromosome(const G4String& name,
                                           G4ThreeVector pos,
                                           const G4double& sx,
                                           const G4double& sy,
                                           const G4double& sz)
  : VirtualChromosome(name)
  , fCenter(std::move(pos))
  , fSemiX(std::abs(sx))
  , fSemiY(std::abs(sy))
  , fSemiZ(std::abs(sz))
  , fRotation(G4RotationMatrix())
{
  fInverseRotation = fRotation.inverse();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EllipticalChromosome::EllipticalChromosome(
  const G4String& name, G4ThreeVector pos, const G4double& sx,
  const G4double& sy, const G4double& sz, G4RotationMatrix rot)
  : VirtualChromosome(name)
  , fCenter(std::move(pos))
  , fSemiX(std::abs(sx))
  , fSemiY(std::abs(sy))
  , fSemiZ(std::abs(sz))
  , fRotation(std::move(rot))
{
  fInverseRotation = fRotation.inverse();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EllipticalChromosome::~EllipticalChromosome() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool EllipticalChromosome::PointInChromosome(G4ThreeVector const& pos)
{
  G4ThreeVector rpos = pos - fCenter;
  rpos               = fInverseRotation(rpos);

  //====================dousatsu=========================
  G4bool ok = ((std::pow(rpos.getX(), 2) / std::pow(fSemiX, 2) +
                std::pow(rpos.getY(), 2) / std::pow(fSemiY, 2) +
                std::pow(rpos.getZ(), 2) / std::pow(fSemiZ, 2)) <= 1.);
  //====================dousatsu=========================
  //// TODO: This is a box, not an ellipsoid
  //// Check that (x/a)**2 + (y/b)**2 + (z/c)**2 <= 1 instead
  // G4bool ok = (std::abs(rpos.getX()) <= fSemiX);
  // ok = ok && (std::abs(rpos.getY()) <= fSemiY);
  // ok = ok && (std::abs(rpos.getZ()) <= fSemiZ);

  return ok;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// NOTE: Not uniformly random
G4ThreeVector EllipticalChromosome::RandomPointInChromosome()
{
  G4ThreeVector point = G4UniformRand() * G4RandomDirection();
  point.setX(point.getX() * fSemiX);
  point.setY(point.getY() * fSemiY);
  point.setZ(point.getZ() * fSemiZ);
  return fCenter + fRotation(point);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
