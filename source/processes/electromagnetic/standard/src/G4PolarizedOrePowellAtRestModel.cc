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
// File name:     G4PolarizedOrePowellAtRestModel
//
// Author:        I.Semeniouk & D.Bernard
//
// Creation date: 26 July 2024
//
// -------------------------------------------------------------------
//

#include "G4PolarizedOrePowellAtRestModel.hh"

#include "G4DynamicParticle.hh"
#include "G4Gamma.hh"
#include "G4Material.hh"
#include "G4PhysicalConstants.hh"
#include "G4RandomDirection.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PolarizedOrePowellAtRestModel::G4PolarizedOrePowellAtRestModel()
  : G4VPositronAtRestModel("OrePowellPolarized")
{
  // DEBUG
  //G4cout << "G4PolarizedOrePowellAtRestModel in contractor" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PolarizedOrePowellAtRestModel::SampleSecondaries(
  std::vector<G4DynamicParticle*>& secParticles, G4double&, const G4Material*) const
{
  const G4double PositronMass = CLHEP::electron_mass_c2;
  const G4double ymax = 22.0;
  const G4double sq2 = std::sqrt(2.);

  CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();

  // OrthoPositronium polarization

  // Random mState in 1:1:1 proportion, must be a user parameter
  G4int  mPos = static_cast<G4int>(std::floor(3.0 * rndmEngine->flat()) - 1.0);

  const G4complex iComplex(0., 1.);

  // The u vector defined in global coordinate system
  //quantization along z - axis, must be a user parameter
  G4complex ux, uy, uz;

  if (mPos == 0) {
    uz = 1.;
    uy = 0.;
    ux = 0.;
  }
  else if (mPos == 1) {
    uz = 0.;
    uy = iComplex / sq2;
    ux = 1. / sq2;
  }
  else if (mPos == -1) {
    uz = 0.;
    uy = -iComplex / sq2;
    ux = 1. / sq2;
  }

  G4double r1, r2, r3;
  G4double cos12, cos13;
  G4double pdf;

  G4double rndmv2[2];
  G4double rndmv1;

  G4ThreeVector PhotonDir1, PhotonDir2, PhotonDir3;
  G4ThreeVector PhotonPol1, PhotonPol2, PhotonPol3;

  // rotate from local to world coordinate;
  G4RotationMatrix LtoW;

  do {
    rndmv1 = rndmEngine->flat();

    do {
      rndmEngine->flatArray(2, rndmv2);
      // energies of photon1 and photon2 normalized to electron rest mass
      r1 = rndmv2[0];
      r2 = rndmv2[1];

      // energy conservation, with positronium assumed = 2 * electron rest mass
      r3 = 2.0 - (r1 + r2);
      // cosine of angles between photons, from momentum conservation
      cos12 = (r3 * r3 - r1 * r1 - r2 * r2) / (2 * r1 * r2);
      cos13 = (r2 * r2 - r1 * r1 - r3 * r3) / (2 * r1 * r3);
      // request both cosines < 1.
    } while (std::abs(cos12) > 1 || std::abs(cos13) > 1);

    G4double sin12 =  std::sqrt((1 + cos12)*(1 - cos12));
    G4double sin13 = -std::sqrt((1 + cos13)*(1 - cos13));

    // Photon directions in the decay plane (y-z), photon 1 along z
    PhotonDir1 = {0., 0., 1.};
    PhotonDir2 = {0., sin12, cos12};
    PhotonDir3 = {0., sin13, cos13};

    // Polarization

    // direction perpendicular to the photon direction, still in the decay plane
    G4ThreeVector PhotonPerp1(0., 1., 0.);
    G4ThreeVector PhotonPerp2(0., cos12, -sin12);
    G4ThreeVector PhotonPerp3(0., cos13, -sin13);

    G4ThreeVector xDir(1., 0., 0.);

    // photon polarization vector (a_i in Ore & Powell)
    G4double eta1 = CLHEP::pi * G4UniformRand();
    PhotonPol1 = std::cos(eta1) * PhotonPerp1 + std::sin(eta1) * xDir;

    G4double eta2 = CLHEP::pi * G4UniformRand();
    PhotonPol2 = std::cos(eta2) * PhotonPerp2 + std::sin(eta2) * xDir;

    G4double eta3 = CLHEP::pi * G4UniformRand();
    PhotonPol3 = std::cos(eta3) * PhotonPerp3 + std::sin(eta3) * xDir;

    // direction perpendicular to the photon direction and to  polarization  (a'_i in Ore & Powell)
    G4ThreeVector PhotonPrime1 = PhotonPol1.cross(PhotonDir1);
    G4ThreeVector PhotonPrime2 = PhotonPol2.cross(PhotonDir2);
    G4ThreeVector PhotonPrime3 = PhotonPol3.cross(PhotonDir3);

    // transition matrix elements (t_i in Ore & Powell)
    G4ThreeVector t1(-PhotonPol3 * (PhotonPol1.dot(PhotonPol2))
                     + PhotonPol1 * (PhotonPrime2.dot(PhotonPrime3))
                     - PhotonPrime2 * (PhotonPrime3.dot(PhotonPol1))
                     - PhotonPrime3 * (PhotonPol1.dot(PhotonPrime2)));
    G4ThreeVector t2(-PhotonPol1 * (PhotonPol2.dot(PhotonPol3))
                     + PhotonPol2 * (PhotonPrime3.dot(PhotonPrime1))
                     - PhotonPrime3 * (PhotonPrime1.dot(PhotonPol2))
                     - PhotonPrime1 * (PhotonPol2.dot(PhotonPrime3)));
    G4ThreeVector t3(-PhotonPol2 * (PhotonPol3.dot(PhotonPol1))
                     + PhotonPol3 * (PhotonPrime1.dot(PhotonPrime2))
                     - PhotonPrime1 * (PhotonPrime2.dot(PhotonPol3))
                     - PhotonPrime2 * (PhotonPol3.dot(PhotonPrime1)));

    // Transformation matrix

    G4ThreeVector z = G4RandomDirection();
    auto tmp = G4RandomDirection();
    G4ThreeVector x = z.cross(tmp).unit();
    G4ThreeVector y = z.cross(x);

    LtoW = G4RotationMatrix(x, y, z);

    // G4cout << x << y << z << G4endl;

    G4ThreeVector t = t1 + t2 + t3;
    t.transform(LtoW); // t in World coordinate system 

    G4complex H = t(0) * ux + t(1) * uy + t(2) * uz;
    pdf = std::abs(conj(H) * H);
  } while (pdf < ymax * rndmv1);
  // END of Sampling

  PhotonDir1.transform(LtoW);
  PhotonDir2.transform(LtoW);
  PhotonDir3.transform(LtoW);

  PhotonPol1.transform(LtoW);
  PhotonPol2.transform(LtoW);
  PhotonPol3.transform(LtoW);

  auto aGamma1 = new G4DynamicParticle(G4Gamma::Gamma(), PhotonDir1, r1 * PositronMass);
  aGamma1->SetPolarization(PhotonPol1);
  secParticles.push_back(aGamma1);

  auto aGamma2 = new G4DynamicParticle(G4Gamma::Gamma(), PhotonDir2, r2 * PositronMass);
  aGamma2->SetPolarization(PhotonPol2);
  secParticles.push_back(aGamma2);

  auto aGamma3 = new G4DynamicParticle(G4Gamma::Gamma(), PhotonDir3, r3 * PositronMass);
  aGamma3->SetPolarization(PhotonPol3);
  secParticles.push_back(aGamma3);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PolarizedOrePowellAtRestModel::PrintGeneratorInformation() const
{
  G4cout << "Polarized Ore Powell AtRest positron 3-gamma annihilation model" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
