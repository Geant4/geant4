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
/// file: CylindricalChromosome.cc
/// brief: Implementation of virt chromosome class for cylindrical chromosomes

#include <utility>

#include "CylindricalChromosome.hh"

#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4String CylindricalChromosome::fShape = "cyl";

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CylindricalChromosome::CylindricalChromosome(const G4String& name,
                                             G4ThreeVector pos,
                                             const G4double& radius,
                                             const G4double& height)
  : VirtualChromosome(name)
  , fCenter(std::move(pos))
  , fRadius(radius)
  , fHeight(height)
  , fRotation(G4RotationMatrix())
{
  fInverseRotation = fRotation.inverse();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CylindricalChromosome::CylindricalChromosome(const G4String& name,
                                             G4ThreeVector pos,
                                             const G4double& radius,
                                             const G4double& height,
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

bool CylindricalChromosome::PointInChromosome(G4ThreeVector const& pos)
{
  G4ThreeVector rpos = pos - fCenter;
  rpos               = fInverseRotation(rpos);

  G4bool height_ok;
  G4bool radius_ok;
  G4double height;
  G4double rad2;

  height = std::abs(rpos.getZ());
  rad2   = rpos.getX() * rpos.getX() + rpos.getY() * rpos.getY();

  height_ok = (height < fHeight);
  radius_ok = (rad2 < (fRadius * fRadius));

  return height_ok && radius_ok;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector CylindricalChromosome::RandomPointInChromosome()
{
  G4ThreeVector point;
  G4double z     = 2 * (G4UniformRand() - 0.5) * fHeight;
  G4double theta = twopi * G4UniformRand();
  G4double r     = fRadius * std::pow(G4UniformRand(), 0.5);
  G4double x     = r * std::cos(theta);
  G4double y     = r * std::sin(theta);
  point          = G4ThreeVector(x, y, z);
  point          = fRotation(point) + fCenter;

  return point;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
