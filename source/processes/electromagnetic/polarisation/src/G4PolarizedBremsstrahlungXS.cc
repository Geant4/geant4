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
// Geant4 class file
//
// File name:     G4PolarizedBremsstrahlungXS
//
// Author:        Andreas Schaelicke on the base of Karim Laihems code

#include "G4PolarizedBremsstrahlungXS.hh"

#include "G4PhysicalConstants.hh"

G4double G4PolarizedBremsstrahlungXS::SCRN[2][19] = {
  { 0.5, 1.0, 2.0, 4.0, 8.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0,
    60.0, 70.0, 80.0, 90.0, 100.0, 120.0 },
  { 0.0145, 0.0490, 0.1400, 0.3312, 0.6758, 1.126, 1.367, 1.564, 1.731, 1.875,
    2.001, 2.114, 2.216, 2.393, 2.545, 2.676, 2.793, 2.897, 3.078 }
};

G4PolarizedBremsstrahlungXS::G4PolarizedBremsstrahlungXS()
{
  fFinalLeptonPolarization = G4StokesVector::ZERO;
  fFinalGammaPolarization  = G4StokesVector::ZERO;
}

G4PolarizedBremsstrahlungXS::~G4PolarizedBremsstrahlungXS() {}

void G4PolarizedBremsstrahlungXS::Initialize(G4double aLept0E, G4double aGammaE,
                                             G4double sintheta,
                                             const G4StokesVector& beamPol,
                                             const G4StokesVector& /*p1*/,
                                             G4int /*flag*/)
{
  G4double aLept1E = aLept0E - aGammaE;

  G4double Stokes_S1 = beamPol.x();
  G4double Stokes_S2 = beamPol.y();
  G4double Stokes_S3 = beamPol.z();

  G4double Lept0E  = aLept0E / electron_mass_c2 + 1.;
  G4double Lept0E2 = Lept0E * Lept0E;
  G4double GammaE  = aGammaE / electron_mass_c2;
  G4double GammaE2 = GammaE * GammaE;
  G4double Lept1E  = aLept1E / electron_mass_c2 + 1.;
  G4double Lept1E2 = Lept1E * Lept1E;

  // *******  Gamma Transverse Momentum
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
  else if(delta < 120)
  {
    for(G4int j = 1; j < 19; ++j)
    {
      if(SCRN[0][j] >= delta)
      {
        GG = std::log(2 * Lept0E * Lept1E / GammaE) - 2 - fCoul -
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

  G4double I_Lept = (Lept0E2 + Lept1E2) * (3. + 2. * GG) -
                    2 * Lept0E * Lept1E * (1. + 4. * u2 * Xsi2 * GG);
  G4double F_Lept =
    Lept1E * 4. * GammaE * u * Xsi * (1. - 2 * Xsi) * GG / I_Lept;
  G4double E_Lept =
    Lept0E * 4. * GammaE * u * Xsi * (2. * Xsi - 1.) * GG / I_Lept;
  G4double M_Lept =
    4. * Lept0E * Lept1E * (1. + GG - 2. * Xsi2 * u2 * GG) / I_Lept;
  G4double P_Lept =
    GammaE2 * (1. + 8. * GG * (Xsi - 0.5) * (Xsi - 0.5)) / I_Lept;

  G4double Stokes_SS1 = M_Lept * Stokes_S1 + E_Lept * Stokes_S3;
  G4double Stokes_SS2 = M_Lept * Stokes_S2;
  G4double Stokes_SS3 = (M_Lept + P_Lept) * Stokes_S3 + F_Lept * Stokes_S1;

  fFinalLeptonPolarization.setX(Stokes_SS1);
  fFinalLeptonPolarization.setY(Stokes_SS2);
  fFinalLeptonPolarization.setZ(Stokes_SS3);

  if(fFinalLeptonPolarization.mag2() > 1.)
  {
    G4ExceptionDescription ed;
    ed << " WARNING in pol-brem fFinalLeptonPolarization \n";
    ed << "\t" << fFinalLeptonPolarization << "\t GG\t" << GG << "\t delta\t"
       << delta;
    G4Exception("G4PolarizedBremsstrahlungXS::Initialize", "pol014",
                JustWarning, ed);
    fFinalLeptonPolarization.setX(0);
    fFinalLeptonPolarization.setY(0);
    fFinalLeptonPolarization.setZ(Stokes_SS3);
    if(Stokes_SS3 > 1)
      fFinalLeptonPolarization.setZ(1);
  }

  G4double I_Gamma = (Lept0E2 + Lept1E2) * (3. + 2. * GG) -
                     2. * Lept0E * Lept1E * (1. + 4. * u2 * Xsi2 * GG);
  G4double D_Gamma = 8. * Lept0E * Lept1E * u2 * Xsi2 * GG / I_Gamma;
  G4double L_Gamma = GammaE *
                     ((Lept0E + Lept1E) * (3. + 2. * GG) -
                      2. * Lept1E * (1. + 4. * u2 * Xsi2 * GG)) /
                     I_Gamma;
  G4double T_Gamma =
    4. * GammaE * Lept1E * Xsi * u * (2. * Xsi - 1.) * GG / I_Gamma;

  G4double Stokes_P1 = D_Gamma;
  G4double Stokes_P2 = 0.;
  G4double Stokes_P3 = (Stokes_S3 * L_Gamma + Stokes_S1 * T_Gamma);

  fFinalGammaPolarization.SetPhoton();

  fFinalGammaPolarization.setX(Stokes_P1);
  fFinalGammaPolarization.setY(Stokes_P2);
  fFinalGammaPolarization.setZ(Stokes_P3);

  if(fFinalGammaPolarization.mag2() > 1.)
  {
    G4ExceptionDescription ed;
    ed << " WARNING in pol-brem fFinalGammaPolarization \n";
    ed << "\t" << fFinalGammaPolarization << "\t GG\t" << GG << "\t delta\t"
       << delta;
    G4Exception("G4PolarizedBremsstrahlungXS::Initialize", "pol015",
                JustWarning, ed);
  }
}

G4double G4PolarizedBremsstrahlungXS::XSection(const G4StokesVector& /*pol2*/,
                                               const G4StokesVector& /*pol3*/)
{
  G4ExceptionDescription ed;
  ed << "ERROR dummy routine G4PolarizedBremsstrahlungXS::XSection "
        "called.\n";
  G4Exception("G4PolarizedBremsstrahlungXS::XSection", "pol016", FatalException,
              ed);

  return 0.;
}

// return expected mean polarisation
G4StokesVector G4PolarizedBremsstrahlungXS::GetPol2()
{
  // electron/positron
  return fFinalLeptonPolarization;
}

G4StokesVector G4PolarizedBremsstrahlungXS::GetPol3()
{
  // photon
  return fFinalGammaPolarization;
}
