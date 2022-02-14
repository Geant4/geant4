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
// Geant4 Class file
//
// File name:     G4PolarizedPhotoElectricXS
//
// Author:        Karim Laihem

#include "G4PolarizedPhotoElectricXS.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4PolarizedPhotoElectricXS::G4PolarizedPhotoElectricXS()
{
  fFinalElectronPolarization = G4StokesVector::ZERO;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4PolarizedPhotoElectricXS::~G4PolarizedPhotoElectricXS() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4PolarizedPhotoElectricXS::Initialize(G4double aGammaE, G4double aLept0E,
                                            G4double sinTheta,
                                            const G4StokesVector& beamPol,
                                            const G4StokesVector& /*p1*/,
                                            G4int /*flag*/)
{
  // Polarization transfer to e- in PhotoelectricEffect.
  G4double Gfactor   = aLept0E / CLHEP::electron_mass_c2 + 1.;
  G4double Gfactor_2 = Gfactor * Gfactor;

  G4double beta = std::sqrt(1. - 1. / (Gfactor_2));

  G4double Stokes_P3 = beamPol.z();

  G4double Lept0E  = aLept0E / CLHEP::electron_mass_c2 + 1.;
  G4double Lept0E2 = Lept0E * Lept0E;
  G4double GammaE  = aGammaE / CLHEP::electron_mass_c2;

  G4double cosTheta = std::sqrt(1. - sinTheta * sinTheta);

  G4double I_Lepton0 =
    1.0 +
    (1. / GammaE) * ((2. / (GammaE * Lept0E * (1 - beta * cosTheta))) - 1.);

  G4double A_Lepton0 =
    (Lept0E / (Lept0E + 1)) *
    (2.0 / (GammaE * Lept0E) + beta * cosTheta +
     (2.0 / ((GammaE * Lept0E2) * (1.0 - beta * cosTheta)))) /
    I_Lepton0;

  G4double B_Lepton0 = (Lept0E / (Lept0E + 1.0)) * beta * sinTheta *
                       (2.0 / (GammaE * Lept0E * (1 - beta * cosTheta)) - 1.0) /
                       I_Lepton0;

  fFinalElectronPolarization.setX(Stokes_P3 * B_Lepton0);
  fFinalElectronPolarization.setY(0.);
  fFinalElectronPolarization.setZ(Stokes_P3 * A_Lepton0);

  if((fFinalElectronPolarization.x() * fFinalElectronPolarization.x() +
      fFinalElectronPolarization.y() * fFinalElectronPolarization.y() +
      fFinalElectronPolarization.z() * fFinalElectronPolarization.z()) > 1.)

  {
    G4ExceptionDescription ed;
    ed << "Warning: PhotoelectricEffect Problem in pol-transfer photon to "
          "lepton:Px2 + Py2 + Pz2 > 1\n";
    ed << "Polarization transfer forced to be total and similar as incoming "
          "Photo\n";
    G4Exception("G4PolarizedPhotoElectricXS::Initialize", "pol023", JustWarning,
                ed);
    fFinalElectronPolarization = beamPol;  // to be safe
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4PolarizedPhotoElectricXS::XSection(const G4StokesVector& /*pol2*/,
                                              const G4StokesVector& /*pol3*/)
{
  G4ExceptionDescription ed;
  ed << "ERROR dummy routine G4PolarizedPhotoElectricXS::XSection() "
        "called\n";
  G4Exception("G4PolarizedPhotoElectricXS::XSection", "pol024", FatalException,
              ed);
  return 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4StokesVector G4PolarizedPhotoElectricXS::GetPol2()
{
  return fFinalElectronPolarization;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4StokesVector G4PolarizedPhotoElectricXS::GetPol3()
{
  return G4StokesVector();
}
