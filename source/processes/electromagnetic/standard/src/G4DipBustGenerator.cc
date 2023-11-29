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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4DipBustGenerator
//
// Author:  Vladimir Grichine
//
// Creation date: 17 May 2011
//
// Modifications:
//
// 17.07.2018  optimisations, using G4Pow in SampleCosTheta() method  M.Novak
//
// Class Description:
//
// Bremsstrahlung Angular Distribution Generation
// suggested  the dipole approximation in the rest frame of electron
// busted in the laboratory frame.
//
// Class Description: End
//
// -------------------------------------------------------------------
//

#include "G4DipBustGenerator.hh"
#include "Randomize.hh"
// #include "G4Log.hh"
// #include "G4Exp.hh"
#include "G4Pow.hh"
#include <CLHEP/Units/PhysicalConstants.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DipBustGenerator::G4DipBustGenerator(const G4String&)
  : G4VEmAngularDistribution("DipBustGen")
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DipBustGenerator::~G4DipBustGenerator() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DipBustGenerator::SampleCosTheta(const G4double kinEnergy)
{
  const G4double c        = 4. - 8.*G4UniformRand();
  const G4double delta    = 0.5*(std::sqrt(c*c+4.) + std::abs(c));
  const G4double signc    = (c < 0.) ? -1.0 : 1.0;
  // const G4double cofA = -signc*G4Exp(G4Log(delta)/3.0);
  const G4double cofA     = -signc*G4Pow::GetInstance()->A13(delta);
  const G4double cosTheta = std::min(1.,std::max(-1.,cofA - 1./cofA));
  const G4double tau      = kinEnergy/CLHEP::electron_mass_c2;
  const G4double beta     = std::sqrt(tau*(tau + 2.))/(tau + 1.);

  return (cosTheta + beta)/(1. + cosTheta*beta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector&
G4DipBustGenerator::SampleDirection(const G4DynamicParticle* dp, G4double,
                                     G4int, const G4Material*)
{
  const G4double cosTheta = SampleCosTheta(dp->GetKineticEnergy());

  const G4double sinTheta = std::sqrt((1. - cosTheta)*(1. + cosTheta));
  const G4double phi      = CLHEP::twopi*G4UniformRand();

  fLocalDirection.set(sinTheta*std::cos(phi), sinTheta*std::sin(phi),cosTheta);
  fLocalDirection.rotateUz(dp->GetMomentumDirection());

  return fLocalDirection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DipBustGenerator::PolarAngle(G4double eTkin,
					                               G4double, // final_energy
					                               G4int ) // Z
{
  const G4double cosTheta = SampleCosTheta(eTkin);
  G4double theta = std::acos(cosTheta);
  theta = std::min(std::max(theta, 0.), CLHEP::pi);
  return theta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DipBustGenerator::SamplePairDirections(const G4DynamicParticle* dp,
					                                     G4double elecKinEnergy,
					                                     G4double posiKinEnergy,
					                                     G4ThreeVector& dirElectron,
					                                     G4ThreeVector& dirPositron,
					                                     G4int, const G4Material*)
{
  const G4double phi  = CLHEP::twopi * G4UniformRand();
  const G4double sinp = std::sin(phi);
  const G4double cosp = std::cos(phi);

  G4double cost = SampleCosTheta(elecKinEnergy);
  G4double sint = std::sqrt((1. - cost)*(1. + cost));

  dirElectron.set(sint*cosp, sint*sinp, cost);
  dirElectron.rotateUz(dp->GetMomentumDirection());

  cost = SampleCosTheta(posiKinEnergy);
  sint = std::sqrt((1. - cost)*(1. + cost));

  dirPositron.set(-sint*cosp, -sint*sinp, cost);
  dirPositron.rotateUz(dp->GetMomentumDirection());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DipBustGenerator::PrintGeneratorInformation() const
{
  G4cout << "\n" << G4endl;
  G4cout << "Angular Generator based on classical formula from" << G4endl;
  G4cout << "J.D. Jackson, Classical Electrodynamics, Wiley, New York 1975"
	 << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
