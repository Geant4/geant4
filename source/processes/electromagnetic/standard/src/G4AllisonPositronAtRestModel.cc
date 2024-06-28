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
// GEANT4 Class file
//
//
// File name:     G4AllisonPositronAtRestModel
//
// Author:        Vladimir Ivanchenko
// 
// Creation date: 14 May 2024
//
// -------------------------------------------------------------------
//

#include "G4AllisonPositronAtRestModel.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "Randomize.hh"
#include "G4Gamma.hh"
#include "G4RandomDirection.hh"
#include "G4ThreeVector.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4AllisonPositronAtRestModel::G4AllisonPositronAtRestModel()
  : G4VPositronAtRestModel("Allison")
{}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4AllisonPositronAtRestModel::SampleSecondaries(
             std::vector<G4DynamicParticle*>& secParticles,
             G4double&, const G4Material* material) const
{
  const G4double eGamma = CLHEP::electron_mass_c2;

  // In rest frame of positronium gammas are back to back
  const G4ThreeVector& dir1 = G4RandomDirection();
  const G4ThreeVector& dir2 = -dir1;
  auto aGamma1 = new G4DynamicParticle(G4Gamma::Gamma(),dir1,eGamma);
  auto aGamma2 = new G4DynamicParticle(G4Gamma::Gamma(),dir2,eGamma);

  // In rest frame the gammas are polarised perpendicular to each other - see
  // Pryce and Ward, Nature No 4065 (1947) p.435.
  // Snyder et al, Physical Review 73 (1948) p.440.
  G4ThreeVector pol1 = (G4RandomDirection().cross(dir1)).unit();
  G4ThreeVector pol2 = (pol1.cross(dir2)).unit();

  // A positron in matter slows down and combines with an atomic electron to
  // make a neutral atom called positronium, about half the size of a normal
  // atom. I expect that when the energy of the positron is small enough,
  // less than the binding energy of positronium (6.8 eV), it is
  // energetically favourable for an electron from the outer orbitals of a
  // nearby atom or molecule to transfer and bind to the positron, as in an
  // ionic bond, leaving behind a mildly ionised nearby atom/molecule. I
  // would expect the positronium to come away with a kinetic energy of a
  // few eV on average. In its para (spin 0) state it annihilates into two
  // photons, which in the rest frame of the positronium are collinear
  // (back-to-back) due to momentum conservation. Because of the motion of the
  // positronium, photons will be not quite back-to-back in the laboratory.

  // The positroniuim acquires an energy of order its binding energy and
  // doesn't have time to thermalise. Nevertheless, here we approximate its
  // energy distribution by a Maxwell-Boltzman with mean energy <KE>. In terms
  // of a more familiar concept of temperature, and the law of equipartition
  // of energy of translational motion, <KE>=3kT/2. Each component of velocity
  // has a distribution exp(-mv^2/2kT), which is a Gaussian of mean zero
  // and variance kT/m=2<KE>/3m, where m is the positronium mass.

  const G4double meanEnergyPerIonPair = material->GetIonisation()->GetMeanEnergyPerIonPair();
  const G4double& meanKE = meanEnergyPerIonPair;  // Just an alias

  if (meanKE > 0.) {  // Positronium has motion
    // Mass of positronium
    const G4double mass = 2.*CLHEP::electron_mass_c2;
    // Mean <KE>=3kT/2, as described above
    // const G4double T = 2.*meanKE/(3.*k_Boltzmann);
    // Component velocities: Gaussian, variance kT/m=2<KE>/3m.
    const G4double sigmav = std::sqrt(2.*meanKE/(3.*mass));
    // This is in units where c=1
    const G4double vx = G4RandGauss::shoot(0.,sigmav);
    const G4double vy = G4RandGauss::shoot(0.,sigmav);
    const G4double vz = G4RandGauss::shoot(0.,sigmav);
    const G4ThreeVector v(vx,vy,vz);  // In unit where c=1
    const G4ThreeVector& beta = v;    // so beta=v/c=v
    aGamma1->Set4Momentum(aGamma1->Get4Momentum().boost(beta));
    aGamma2->Set4Momentum(aGamma2->Get4Momentum().boost(beta));

    // Rotate polarisation vectors
    const G4ThreeVector& newDir1 = aGamma1->GetMomentumDirection();
    const G4ThreeVector& newDir2 = aGamma2->GetMomentumDirection();
    const G4ThreeVector& axis1 = dir1.cross(newDir1);  // No need to be unit
    const G4ThreeVector& axis2 = dir2.cross(newDir2);  // No need to be unit
    const G4double& angle1 = std::acos(dir1*newDir1);
    const G4double& angle2 = std::acos(dir2*newDir2);
    pol1.rotate(axis1, angle1);
    pol2.rotate(axis2, angle2);
  }

  // use constructors optimal for massless particle
  aGamma1->SetPolarization(pol1);
  aGamma2->SetPolarization(pol2);

  secParticles.push_back(aGamma1);
  secParticles.push_back(aGamma2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4AllisonPositronAtRestModel::PrintGeneratorInformation() const
{
  G4cout << "\n" << G4endl;
  G4cout << "Allison AtRest positron 2-gamma annihilation model." << G4endl;
  G4cout << "Takes into account positronium motion in the media." << G4endl;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
