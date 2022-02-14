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
// -------------------------------------------------------------------
//
// Geant4 Class file
//
// File name:     G4PolarizedGammaConversionXS
//
// Author:        Andreas Schaelicke on the base of Karim Laihems code

#include "G4PolarizedGammaConversionXS.hh"

G4double G4PolarizedGammaConversionXS::SCRN[2][19] = {
  { 0.5, 1.0, 2.0, 4.0, 8.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0,
    60.0, 70.0, 80.0, 90.0, 100.0, 120.0 },
  { 0.0145, 0.0490, 0.1400, 0.3312, 0.6758, 1.126, 1.367, 1.564, 1.731, 1.875,
    2.001, 2.114, 2.216, 2.393, 2.545, 2.676, 2.793, 2.897, 3.078 }
};

G4PolarizedGammaConversionXS::G4PolarizedGammaConversionXS()
{
  fFinalElectronPolarization = G4StokesVector::ZERO;
  fFinalPositronPolarization = G4StokesVector::ZERO;
}

G4PolarizedGammaConversionXS::~G4PolarizedGammaConversionXS() {}

void G4PolarizedGammaConversionXS::Initialize(
  G4double aGammaE, G4double aLept0E, G4double sintheta,
  const G4StokesVector& beamPol, const G4StokesVector& /*p1*/, G4int /*flag*/)
{
  G4double aLept1E = aGammaE - aLept0E;

  G4double Stokes_P3 = beamPol.z();

  G4double Lept0E  = aLept0E / CLHEP::electron_mass_c2 + 1.;
  G4double Lept0E2 = Lept0E * Lept0E;
  G4double GammaE  = aGammaE / CLHEP::electron_mass_c2;
  G4double Lept1E  = aLept1E / CLHEP::electron_mass_c2 - 1.;
  G4double Lept1E2 = Lept1E * Lept1E;

  // *******  Gamma Transvers Momentum
  G4double TMom = std::sqrt(Lept0E2 - 1.) * sintheta;
  G4double u    = TMom;
  G4double u2   = u * u;
  G4double Xsi  = 1. / (1. + u2);
  G4double Xsi2 = Xsi * Xsi;

  G4double delta =
    12. * std::pow(fZ, 1. / 3.) * Lept0E * Lept1E * Xsi / (121. * GammaE);
  G4double GG = 0.;

  if(delta < 0.5)
  {
    GG = std::log(2. * Lept0E * Lept1E / GammaE) - 2. - fCoul;
  }
  else if(delta < 120.)
  {
    for(G4int j = 1; j < 19; ++j)
    {
      if(SCRN[0][j] >= delta)
      {
        GG = std::log(2. * Lept0E * Lept1E / GammaE) - 2. - fCoul -
             (SCRN[1][j - 1] + (delta - SCRN[0][j - 1]) *
                                 (SCRN[1][j] - SCRN[1][j - 1]) /
                                 (SCRN[0][j] - SCRN[0][j - 1]));
        break;
      }
    }
  }
  else
  {
    G4double alpha_sc = (111. * std::pow(fZ, -1. / 3.)) / Xsi;
    GG                = std::log(alpha_sc) - 2. - fCoul;
  }

  if(GG < -1.)
    GG = -1.;

  G4double I_Lepton = (Lept0E2 + Lept1E2) * (3 + 2 * GG) +
                      2. * Lept0E * Lept1E * (1. + 4. * u2 * Xsi2 * GG);

  G4double L_Lepton1 = GammaE *
                       ((Lept0E - Lept1E) * (3. + 2. * GG) +
                        2 * Lept1E * (1. + 4. * u2 * Xsi2 * GG)) /
                       I_Lepton;

  G4double T_Lepton1 =
    4. * GammaE * Lept1E * Xsi * u * (1. - 2. * Xsi) * GG / I_Lepton;

  G4double Stokes_S1 = (Stokes_P3 * T_Lepton1);
  G4double Stokes_S2 = 0.;
  G4double Stokes_S3 = (Stokes_P3 * L_Lepton1);

  fFinalElectronPolarization.setX(Stokes_S1);
  fFinalElectronPolarization.setY(Stokes_S2);
  fFinalElectronPolarization.setZ(Stokes_S3);

  if(fFinalElectronPolarization.mag2() > 1.)
  {
    G4ExceptionDescription ed;
    ed << "\t" << fFinalElectronPolarization << "\t GG\t" << GG << "\t delta\t"
       << delta << "\n";
    G4Exception("G4PolarizedGammaConversionXS::Initialize", "pol022",
                JustWarning, ed);
    fFinalElectronPolarization.setX(0.);
    fFinalElectronPolarization.setY(0.);
    fFinalElectronPolarization.setZ(Stokes_S3);
    if(Stokes_S3 > 1.)
      fFinalElectronPolarization.setZ(1.);
  }

  G4double L_Lepton2 = GammaE *
                       ((Lept1E - Lept0E) * (3. + 2. * GG) +
                        2 * Lept0E * (1. + 4. * u2 * Xsi2 * GG)) /
                       I_Lepton;

  G4double T_Lepton2 =
    4. * GammaE * Lept0E * Xsi * u * (1. - 2. * Xsi) * GG / I_Lepton;

  G4double Stokes_SS1 = (Stokes_P3 * T_Lepton2);
  G4double Stokes_SS2 = 0.;
  G4double Stokes_SS3 = (Stokes_P3 * L_Lepton2);

  fFinalPositronPolarization.SetPhoton();

  fFinalPositronPolarization.setX(Stokes_SS1);
  fFinalPositronPolarization.setY(Stokes_SS2);
  fFinalPositronPolarization.setZ(Stokes_SS3);

  if(fFinalPositronPolarization.mag2() > 1.)
  {
    G4ExceptionDescription ed;
    ed << "\t" << fFinalPositronPolarization << "\t GG\t" << GG << "\t delta\t"
       << delta << "\n";
    G4Exception("G4PolarizedGammaConversionXS::Initialize", "pol023",
                JustWarning, ed);
  }
}

G4double G4PolarizedGammaConversionXS::XSection(const G4StokesVector& /*pol2*/,
                                                const G4StokesVector& /*pol3*/)
{
  G4ExceptionDescription ed;
  ed << "ERROR dummy routine G4PolarizedGammaConversionXS::XSection "
        "called \n";
  G4Exception("G4PolarizedGammaConversionXS::Initialize", "pol024",
              FatalException, ed);
  return 0.;
}

// return expected mean polarisation
G4StokesVector G4PolarizedGammaConversionXS::GetPol2()
{
  // electron/positron
  return fFinalElectronPolarization;
}

G4StokesVector G4PolarizedGammaConversionXS::GetPol3()
{
  // photon
  return fFinalPositronPolarization;
}
