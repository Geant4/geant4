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
// File name:     G4PolarizedAnnihilationXS
//
// Author:        Andreas Schaelicke
//
// Class Description:
//   * calculates the differential cross section in e+e- -> gamma gamma

#include "G4PolarizedAnnihilationXS.hh"

#include "G4PhysicalConstants.hh"

G4PolarizedAnnihilationXS::G4PolarizedAnnihilationXS()
  : polxx(0.)
  , polyy(0.)
  , polzz(0.)
  , polxz(0.)
  , polzx(0.)
  , polxy(0.)
  , polyx(0.)
  , polyz(0.)
  , polzy(0.)
  , fPhi0(0.)
{
  fPhi2  = G4ThreeVector(0., 0., 0.);
  fPhi3  = G4ThreeVector(0., 0., 0.);
  fDice  = 0.;
  fPolXS = 0.;
  fUnpXS = 0.;
  ISPxx = ISPyy = ISPzz = ISPnd = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4PolarizedAnnihilationXS::~G4PolarizedAnnihilationXS() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4PolarizedAnnihilationXS::TotalXS()
{
  // total cross section is sum of
  //  + unpol xsec  "sigma0"
  //  + longitudinal polarised cross section "sigma_zz" times
  //  pol_3(e-)*pol_3(e+)
  //  + transverse contribution "(sigma_xx+sigma_yy)/2" times
  //  pol_T(e-)*pol_T(e+)
  //     (Note: if both beams are transversely polarised, i.e. pol_T(e-)!=0 and
  //      pol_T(e+)!=0, and sigma_xx!=sigma_yy, then the diff. cross section
  //      will exhibit a azimuthal asymmetry even if pol_T(e-)*pol_T(e+)==0)
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4PolarizedAnnihilationXS::Initialize(
  G4double eps, G4double gam,
  G4double,                    // phi
  const G4StokesVector& pol0,  // positron polarization
  const G4StokesVector& pol1,  // electron polarization
  G4int flag)
{
  G4double diffXSFactor = re2 / (gam - 1.);
  DefineCoefficients(pol0, pol1);

  // prepare dicing
  fDice = 0.;
  G4double symmXS =
    0.125 * ((-1. / sqr(gam + 1.)) / sqr(eps) +
             ((sqr(gam) + 4. * gam - 1.) / sqr(gam + 1.)) / eps - 1.);

  G4ThreeVector epsVector(1. / sqr(eps), 1. / eps, 1.);
  G4ThreeVector oneEpsVector(1. / sqr(1. - eps), 1. / (1. - eps), 1.);
  G4ThreeVector sumEpsVector(epsVector + oneEpsVector);
  G4ThreeVector difEpsVector(epsVector - oneEpsVector);
  G4ThreeVector calcVector(0., 0., 0.);

  // temporary variables
  G4double helpVar2 = 0., helpVar1 = 0.;

  // unpolarised contribution
  helpVar1   = (gam * gam + 4. * gam + 1.) / sqr(gam + 1.);
  helpVar2   = -1. / sqr(gam + 1.);
  calcVector = G4ThreeVector(helpVar2, helpVar1, -1.);
  fUnpXS     = 0.125 * calcVector * sumEpsVector;

  // initial particles polarised contribution
  helpVar2   = 1. / sqr(gam + 1.);
  helpVar1   = -(gam * gam + 4. * gam + 1.) / sqr(gam + 1.);
  calcVector = G4ThreeVector(helpVar2, helpVar1, 0.5 * (gam + 3.));
  ISPxx      = 0.25 * (calcVector * sumEpsVector) / (gam - 1.);

  helpVar1   = 1. / sqr(gam + 1.);
  calcVector = G4ThreeVector(-helpVar1, 2. * gam * helpVar1, -1.);
  ISPyy      = 0.125 * calcVector * sumEpsVector;

  helpVar1   = 1. / (gam - 1.);
  helpVar2   = 1. / sqr(gam + 1.);
  calcVector = G4ThreeVector(
    -(gam * gam + 1.) * helpVar2,
    (gam * gam * (gam + 1.) + 7. * gam + 3.) * helpVar2, -(gam + 3.));
  ISPzz = 0.125 * helpVar1 * (calcVector * sumEpsVector);

  helpVar1   = std::sqrt(std::fabs(eps * (1. - eps) * 2. * (gam + 1.) - 1.));
  calcVector = G4ThreeVector(-1. / (gam * gam - 1.), 2. / (gam - 1.), 0.);
  ISPnd      = 0.125 * (calcVector * difEpsVector) * helpVar1;

  fPolXS = ISPxx * polxx + ISPyy * polyy + ISPzz * polzz;
  fPolXS += ISPnd * (polzx + polxz);
  fPhi0 = fUnpXS + fPolXS;
  fDice = symmXS;

  if(polzz != 0.)
  {
    fDice *= (1. + (polzz * ISPzz / fUnpXS));
    if(fDice < 0.)
      fDice = 0.;
  }
  // prepare final state coefficients
  if(flag == 2)
  {
    // circular polarisation
    G4double circ1 = 0., circ2 = 0., circ3 = 0.;
    helpVar1 = 8. * sqr(1. - eps) * sqr(eps) * (gam - 1.) * sqr(gam + 1.) /
               std::sqrt(gam * gam - 1.);
    helpVar2 = sqr(gam + 1.) * sqr(eps) * (-2. * eps + 3.) -
               (gam * gam + 3. * gam + 2.) * eps;
    circ1 = helpVar2 + gam;
    circ1 /= helpVar1;
    circ2 = helpVar2 + 1.;
    circ2 /= helpVar1;
    helpVar1 = std::sqrt(std::fabs(eps * (1. - eps) * 2. * (gam + 1.) - 1.));
    helpVar1 /= std::sqrt(gam * gam - 1.);
    calcVector = G4ThreeVector(1., -2. * gam, 0.);
    circ3      = 0.125 * (calcVector * sumEpsVector) / (gam + 1.);
    circ3 *= helpVar1;

    fPhi2.setZ(circ2 * pol1.z() + circ1 * pol0.z() +
               circ3 * (pol1.x() + pol0.x()));
    fPhi3.setZ(-circ1 * pol1.z() - circ2 * pol0.z() -
               circ3 * (pol1.x() + pol0.x()));

    // common to both linear polarisation
    calcVector          = G4ThreeVector(-1., 2. * gam, 0.);
    G4double linearZero = 0.125 * (calcVector * sumEpsVector) / sqr(gam + 1.);

    //        Linear Polarisation #1
    helpVar1 = std::sqrt(std::fabs(2. * (gam + 1.) * (1. - eps) * eps - 1.)) /
               ((gam + 1.) * eps * (1. - eps));
    helpVar2 = helpVar1 * helpVar1;

    // photon 1
    G4double diagContrib = 0.125 * helpVar2 * (polxx + polyy - polzz);
    G4double nonDiagContrib =
      0.125 * helpVar1 * (-polxz / (1. - eps) + polzx / eps);

    fPhi2.setX(linearZero + diagContrib + nonDiagContrib);

    // photon 2
    nonDiagContrib = 0.125 * helpVar1 * (polxz / eps - polzx / (1. - eps));

    fPhi3.setX(linearZero + diagContrib + nonDiagContrib);

    //        Linear Polarisation #2
    helpVar1 =
      std::sqrt(gam * gam - 1.) * (2. * (gam + 1.) * eps * (1. - eps) - 1.);
    helpVar1 /= 8. * sqr(1. - eps) * sqr(eps) * sqr(gam + 1.) * (gam - 1.);
    helpVar2 = std::sqrt((gam * gam - 1.) *
                         std::fabs(2. * (gam + 1.) * eps * (1. - eps) - 1.));
    helpVar2 /= 8. * sqr(1. - eps) * sqr(eps) * sqr(gam + 1.) * (gam - 1.);

    G4double contrib21 = (-polxy + polyx) * helpVar1;
    G4double contrib32 =
      -(eps * (gam + 1.) - 1.) * polyz + (eps * (gam + 1.) - gam) * polzy;

    contrib32 *= helpVar2;
    fPhi2.setY(contrib21 + contrib32);

    contrib32 =
      -(eps * (gam + 1.) - gam) * polyz + (eps * (gam + 1.) - 1.) * polzy;
    contrib32 *= helpVar2;
    fPhi3.setY(contrib21 + contrib32);
  }
  fPhi0 *= diffXSFactor;
  fPhi2 *= diffXSFactor;
  fPhi3 *= diffXSFactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4PolarizedAnnihilationXS::XSection(const G4StokesVector& pol2,
                                             const G4StokesVector& pol3)
{
  G4double xs = fPhi0 + pol2 * fPhi2 + pol3 * fPhi3;
  return xs;
}

// calculate total cross section
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4PolarizedAnnihilationXS::TotalXSection(G4double, G4double,
                                                  G4double gam,
                                                  const G4StokesVector& pol0,
                                                  const G4StokesVector& pol1)
{
  G4double totalXSFactor = pi * re2 / (gam + 1.);  // atomic number ignored
  DefineCoefficients(pol0, pol1);

  G4double xs = 0.;

  G4double gam2     = gam * gam;
  G4double sqrtgam1 = std::sqrt(gam2 - 1.);
  G4double logMEM   = std::log(gam + sqrtgam1);
  G4double unpME    = (gam * (gam + 4.) + 1.) * logMEM;
  unpME += -(gam + 3.) * sqrtgam1;
  unpME /= 4. * (gam2 - 1.);
  G4double longPart = (3 + gam * (gam * (gam + 1.) + 7.)) * logMEM;
  longPart += -(5. + gam * (3 * gam + 4.)) * sqrtgam1;
  longPart /= 4. * sqr(gam - 1.) * (gam + 1.);
  G4double tranPart = -(5 * gam + 1.) * logMEM;
  tranPart += (gam + 5.) * sqrtgam1;
  tranPart /= 4. * sqr(gam - 1.) * (gam + 1.);

  xs += unpME;
  xs += polzz * longPart;
  xs += (polxx + polyy) * tranPart;

  return xs * totalXSFactor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4StokesVector G4PolarizedAnnihilationXS::GetPol2()
{
  // Note, mean polarization can not contain correlation effects.
  return G4StokesVector(1. / fPhi0 * fPhi2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4StokesVector G4PolarizedAnnihilationXS::GetPol3()
{
  // Note, mean polarization can not contain correlation effects.
  return G4StokesVector(1. / fPhi0 * fPhi3);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4PolarizedAnnihilationXS::DefineCoefficients(const G4StokesVector& pol0,
                                                   const G4StokesVector& pol1)
{
  polxx = pol0.x() * pol1.x();
  polyy = pol0.y() * pol1.y();
  polzz = pol0.z() * pol1.z();

  polxz = pol0.x() * pol1.z();
  polzx = pol0.z() * pol1.x();

  polyz = pol0.y() * pol1.z();
  polzy = pol0.z() * pol1.y();

  polxy = pol0.x() * pol1.y();
  polyx = pol0.y() * pol1.x();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4PolarizedAnnihilationXS::GetXmin(G4double y)
{
  return 0.5 * (1. - std::sqrt((y - 1.) / (y + 1.)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4PolarizedAnnihilationXS::GetXmax(G4double y)
{
  return 0.5 * (1. + std::sqrt((y - 1.) / (y + 1.)));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4PolarizedAnnihilationXS::DiceEpsilon() { return fDice; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4PolarizedAnnihilationXS::getVar(G4int choice)
{
  if(choice == -1)
    return fPolXS / fUnpXS;
  if(choice == 0)
    return fUnpXS;
  if(choice == 1)
    return ISPxx;
  if(choice == 2)
    return ISPyy;
  if(choice == 3)
    return ISPzz;
  if(choice == 4)
    return ISPnd;
  return 0;
}
