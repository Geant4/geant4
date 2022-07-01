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
/// file: RodChromosome.cc
/// brief: Implementation of virt chromosome class for rod chromosomes

#include <utility>

#include "RodChromosome.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4String RodChromosome::fShape = "rod";

RodChromosome::RodChromosome(const G4String& name, G4ThreeVector pos,
                             const G4double& radius, const G4double& height)
  : VirtualChromosome(name)
  , fCenter(std::move(pos))
  , fRadius(radius)
  , fHeight(height)
  , fRotation(G4RotationMatrix())
{
  fInverseRotation = fRotation.inverse();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RodChromosome::RodChromosome(const G4String& name, G4ThreeVector pos,
                             const G4double& radius, const G4double& height,
                             G4RotationMatrix rot)
  : VirtualChromosome(name)
  , fCenter(std::move(pos))
  , fRadius(radius)
  , fHeight(height)
  , fRotation(std::move(rot))
{
  fInverseRotation = fRotation.inverse();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RodChromosome::~RodChromosome() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool RodChromosome::PointInChromosome(G4ThreeVector const& pos)
{
  G4ThreeVector rpos = pos - fCenter;
  rpos               = fInverseRotation(rpos);

  // in the cylinder?
  bool height_ok;
  bool radius_ok;
  G4double height;
  G4double rad2;

  height = std::abs(rpos.getZ());
  rad2   = rpos.getX() * rpos.getX() + rpos.getY() * rpos.getY();

  height_ok = (height < fHeight);
  radius_ok = (rad2 < (fRadius * fRadius));

  G4bool in_cylinder = height_ok && radius_ok;
  // in the circles at each end?
  G4ThreeVector pos1 = rpos - G4ThreeVector(0, 0, height);
  G4bool in_sphere1  = (pos1.mag2() < (fRadius * fRadius));
  G4ThreeVector pos2 = rpos + G4ThreeVector(0, 0, height);
  G4bool in_sphere2  = (pos2.mag2() < (fRadius * fRadius));

  return (in_cylinder || in_sphere1 || in_sphere2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// NOTE: Not uniformly random
G4ThreeVector RodChromosome::RandomPointInChromosome()
{
  G4ThreeVector point;
  if(G4UniformRand() < 0.5)
  {
    // in cylinder
    G4double z     = 2 * (G4UniformRand() - 0.5) * fHeight;
    G4double theta = twopi * G4UniformRand();
    G4double r     = fRadius * std::pow(G4UniformRand(), 0.5);
    G4double x     = r * std::cos(theta);
    G4double y     = r * std::sin(theta);
    point          = G4ThreeVector(x, y, z);
  }
  else
  {
    // in spheres
    point = fRadius * std::pow(G4UniformRand(), 0.5) * G4RandomDirection();
    if(point.getZ() < 0)
    {
      // bottom half
      point.setZ(point.getZ() - fHeight);
    }
    else
    {
      // top half
      point.setZ(point.getZ() + fHeight);
    }
  }
  return fRotation(point) + fCenter;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
