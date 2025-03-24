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
// File name:     G4OrePowellAtRestModel
//
// Author:        I.Semeniouk & D.Bernard
//
// Creation date: 04 Juin 2024
//
// -------------------------------------------------------------------
//

#include "G4OrePowellAtRestModel.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "Randomize.hh"
#include "G4Gamma.hh"
#include "G4RandomDirection.hh"
#include "G4ThreeVector.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4OrePowellAtRestModel::G4OrePowellAtRestModel() : G4VPositronAtRestModel("OrePowell") {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4OrePowellAtRestModel::SampleSecondaries(
             std::vector<G4DynamicParticle*>& secParticles,
             G4double&, const G4Material*) const
{
  static const G4double PositronMass = CLHEP::electron_mass_c2;
  const G4double ymax = 8.1;

  CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();

  G4double cos12;
  G4double cos13;
  G4double sin12;
  G4double sin13;
  G4double r1;
  G4double r2;
  G4double r3;
  G4double pdf;

  G4double rndmv2[2];
  G4double rndmv1;

  do {
    rndmv1 = rndmEngine->flat();

    do {
      rndmEngine->flatArray(2, rndmv2);
      // energies of photon1 and photon2 normalized to electron rest mass
      r1 = rndmv2[0];
      r2 = rndmv2[1];
      // energy conservation, with positronium assumed = 2 * electron rest mass
      r3 = 2.0 - (r1+r2);
      // cosine of angles between photons, from momentum conservation
      cos12=(r3*r3 - r1*r1 -r2*r2)/(2*r1*r2);
      cos13=(r2*r2 - r1*r1 -r3*r3)/(2*r1*r3);
      // request both cosines < 1.
    } while ( std::abs(cos12) > 1 || std::abs(cos13) > 1 );

    sin12 =  std::sqrt((1 + cos12)*(1 - cos12));
    sin13 = -std::sqrt((1 + cos13)*(1 - cos13));

    G4double cos23=cos12*cos13+sin12*sin13;

    pdf = (1 - cos12)*(1 - cos12) + (1 - cos13)*(1 - cos13) + (1 - cos23)*(1 - cos23);
  } while ( pdf < ymax * rndmv1 );
  // END of Sampling

  // photon directions in the decay plane, photon 1 along z, x perp to the plane.
  G4ThreeVector PhotonMomentum1(0., 0., 1.);
  G4ThreeVector PhotonMomentum2(0.,sin12,cos12);
  G4ThreeVector PhotonMomentum3(0.,sin13,cos13);

  // Random x direction ( rotate decay plane along Z axis)

  G4double phi = CLHEP::twopi * G4UniformRand();
  PhotonMomentum2.rotateZ(phi);
  PhotonMomentum3.rotateZ(phi);


  // First Gamma direction
  G4ThreeVector dir1 = G4RandomDirection();

  PhotonMomentum1.rotateUz(dir1);
  PhotonMomentum2.rotateUz(dir1);
  PhotonMomentum3.rotateUz(dir1);

  auto aGamma1 = new G4DynamicParticle(G4Gamma::Gamma(), PhotonMomentum1,
				       r1 * PositronMass);

  //Random polarization
  G4double phi1 = CLHEP::twopi * G4UniformRand();
  G4ThreeVector pol1(std::cos(phi1),std::sin(phi1),0.0);
  pol1.rotateUz(PhotonMomentum1);
  aGamma1->SetPolarization(pol1);
  secParticles.push_back(aGamma1);

  auto aGamma2 = new G4DynamicParticle(G4Gamma::Gamma(), PhotonMomentum2,
				       r2 * PositronMass);

  G4double phi2 = CLHEP::twopi * G4UniformRand();
  G4ThreeVector pol2(std::cos(phi2),std::sin(phi2),0.0);
  pol2.rotateUz(PhotonMomentum2);
  aGamma2->SetPolarization(pol2);
  secParticles.push_back(aGamma2);

  auto aGamma3 = new G4DynamicParticle(G4Gamma::Gamma(), PhotonMomentum3,
				       r3 * PositronMass);

  G4double phi3 = CLHEP::twopi * G4UniformRand();
  G4ThreeVector pol3(std::cos(phi3),std::sin(phi3),0.0);
  pol3.rotateUz(PhotonMomentum3);
  aGamma3->SetPolarization(pol3);
  secParticles.push_back(aGamma3);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4OrePowellAtRestModel::PrintGeneratorInformation() const
{
  G4cout << "Orel Powell AtRest positron 3-gamma annihilation model" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
