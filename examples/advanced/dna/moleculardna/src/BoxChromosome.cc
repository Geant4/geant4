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
/// file: BoxChromosome.cc
/// brief: Implementation of virt chromosome class for Box chromosomes

#include "BoxChromosome.hh"

#include "G4RandomDirection.hh"
#include "Randomize.hh"

#include <utility>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4String BoxChromosome::fShape = "box";

BoxChromosome::BoxChromosome(const G4String& name, G4ThreeVector pos, const G4double& xdim,
                             const G4double& ydim, const G4double& zdim)
  : VirtualChromosome(name),
    fCenter(std::move(pos)),
    fXdim(xdim),
    fYdim(ydim),
    fZdim(zdim),
    fRotation(G4RotationMatrix())
{
  fInverseRotation = fRotation.inverse();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BoxChromosome::BoxChromosome(const G4String& name, G4ThreeVector pos, const G4double& xdim,
                             const G4double& ydim, const G4double& zdim, G4RotationMatrix rot)
  : VirtualChromosome(name),
    fCenter(std::move(pos)),
    fXdim(xdim),
    fYdim(ydim),
    fZdim(zdim),
    fRotation(std::move(rot))
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BoxChromosome::~BoxChromosome() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool BoxChromosome::PointInChromosome(G4ThreeVector const& pos)
{
  G4ThreeVector rpos = pos - fCenter;
  rpos = fInverseRotation(rpos);
  G4bool ok = (std::abs(rpos.getX()) <= fXdim);
  ok = ok && (std::abs(rpos.getY()) <= fYdim);
  ok = ok && (std::abs(rpos.getZ()) <= fZdim);
  return ok;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector BoxChromosome::RandomPointInChromosome()
{
  G4ThreeVector point = G4UniformRand() * G4RandomDirection();
  point.setX(point.getX() * fXdim);
  point.setY(point.getY() * fYdim);
  point.setZ(point.getZ() * fZdim);
  return fRotation(point) + fCenter;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
