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
// File name:     G4PolarizedComptonXS
//
// Author:        Andreas Schaelicke
//
// Class Description:
//   determine the  polarization of the final state in a Compton scattering
//   process employing the differential cross section by F.W.Lipps & H.A.Tolhoek
//   ( Physica 20 (1954) 395 )
//   recalculated by P.Starovoitov

#include "G4PolarizedComptonXS.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4PolarizedComptonXS::G4PolarizedComptonXS()
{
  SetYmin(0.);

  fPhi0  = 0.;
  fPolXS = 0.;
  fUnpXS = 0.;
  fPhi2  = G4ThreeVector(0., 0., 0.);
  fPhi3  = G4ThreeVector(0., 0., 0.);
  polxx = polyy = polzz = polxz = polzx = polyz = polzy = polxy = polyx = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4PolarizedComptonXS::~G4PolarizedComptonXS() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4PolarizedComptonXS::Initialize(G4double eps, G4double X,
                                      G4double,  // phi
                                      const G4StokesVector& pol0,
                                      const G4StokesVector& pol1, G4int flag)
{
  G4double cosT = 1. - (1. / eps - 1.) / X;
  if(cosT > 1. + 1.e-8)
    cosT = 1.;
  else if(cosT < -1. - 1.e-8)
    cosT = -1.;
  G4double cosT2 = cosT * cosT;
  G4double cosT3 = cosT2 * cosT;
  G4double sinT2 = 1. - cosT2;
  if(sinT2 > 1. + 1.e-8)
    sinT2 = 1.;
  else if(sinT2 < 0.)
    sinT2 = 0.;
  G4double sinT  = std::sqrt(sinT2);
  G4double cos2T = 2. * cosT2 - 1.;
  G4double sin2T = 2. * sinT * cosT;
  G4double eps2  = sqr(eps);
  DefineCoefficients(pol0, pol1);
  G4double diffXSFactor = re2 / (4. * X);

  // unpolarized Cross Section
  fUnpXS = (eps2 + 1. - eps * sinT2) / (2. * eps);
  // initial polarization dependence
  fPolXS = -sinT2 * pol0.x() + (1. - eps) * sinT * polzx +
           ((eps2 - 1.) / eps) * cosT * polzz;
  fPolXS *= 0.5;

  fPhi0 = fUnpXS + fPolXS;

  if(flag == 2)
  {
    // polarization of outgoing photon
    G4double phi21 = -sinT2 + 0.5 * (cos2T + 3.) * pol0.x() -
                     ((1. - eps) / eps) * sinT * polzx;
    phi21 *= 0.5;
    G4double phi22 = cosT * pol0.y() + ((1. - eps) / (2. * eps)) * sinT * polzy;
    G4double phi23 = ((eps2 + 1.) / eps) * cosT * pol0.z() -
                     ((1. - eps) / eps) * (eps * cosT2 + 1.) * pol1.z();
    phi23 += 0.5 * (1. - eps) * sin2T * pol1.x();
    phi23 += (eps - 1.) * (-sinT2 * polxz + sinT * polyy - 0.5 * sin2T * polxx);
    phi23 *= 0.5;

    fPhi2 = G4ThreeVector(phi21, phi22, phi23);

    // polarization of outgoing electron
    G4double phi32 = -sinT2 * polxy + ((1. - eps) / eps) * sinT * polyz +
                     0.5 * (cos2T + 3.) * pol1.y();
    phi32 *= 0.5;

    G4double phi31    = 0.;
    G4double phi31add = 0.;
    G4double phi33    = 0.;
    G4double phi33add = 0.;

    if((1. - eps) > 1.e-12)
    {
      G4double helpVar = std::sqrt(eps2 - 2. * cosT * eps + 1.);

      phi31 = (1. - eps) * (1. + cosT) * sinT * pol0.z();
      phi31 +=
        (-eps * cosT3 + eps * cosT2 + (eps - 2.) * cosT + eps) * pol1.x();
      phi31 += -(eps * cosT2 - eps * cosT + cosT + 1.) * sinT * pol1.z();
      phi31 /= 2. * helpVar;

      phi31add = -eps * sqr(1. - cosT) * (1. + cosT) * polxx;
      phi31add += (1. - eps) * sinT2 * polyy;
      phi31add += -(-eps2 + cosT * (cosT * eps - eps + 1.) * eps + eps - 1.) *
                  sinT * polxz / eps;
      phi31add /= 2. * helpVar;

      phi33 = ((1. - eps) / eps) *
              (-eps * cosT2 + eps * (eps + 1.) * cosT - 1.) * pol0.z();
      phi33 += -(eps * cosT2 + (1. - eps) * eps * cosT + 1.) * sinT * pol1.x();
      phi33 +=
        -(-eps2 * cosT3 + eps * (eps2 - eps + 1.) * cosT2 - cosT + eps2) *
        pol1.z() / eps;
      phi33 /= -2. * helpVar;

      phi33add = (eps * (eps - cosT - 1.) * cosT + 1.) * sinT * polxx;
      phi33add += -(-eps2 + cosT * eps + eps - 1.) * sinT2 * polxz;
      phi33add += (eps - 1.) * (cosT - eps) * sinT * polyy;
      phi33add /= -2. * helpVar;
    }
    else
    {
      phi31 = -pol1.z() -
              (X - 1.) * std::sqrt(1. - eps) * pol1.x() / std::sqrt(2. * X);
      phi31add = -(-X * X * pol1.z() - 2. * X * (2. * pol0.z() - pol1.z()) -
                   (4. * pol0.x() + 5.) * pol1.z()) *
                 (1. - eps) / (4. * X);

      phi33 = pol1.x() -
              (X - 1.) * std::sqrt(1. - eps) * pol1.z() / std::sqrt(2. * X);
      phi33add = -(X * X - 2. * X + 4. * pol0.x() + 5.) * (1. - eps) *
                 pol1.x() / (4. * X);
    }
    fPhi3 = G4ThreeVector(phi31 + phi31add, phi32, phi33 + phi33add);
  }
  fUnpXS *= diffXSFactor;
  fPolXS *= diffXSFactor;
  fPhi0 *= diffXSFactor;
  fPhi2 *= diffXSFactor;
  fPhi3 *= diffXSFactor;
}

G4double G4PolarizedComptonXS::XSection(const G4StokesVector& pol2,
                                        const G4StokesVector& pol3)
{
  G4bool gammaPol2    = !(pol2 == G4StokesVector::ZERO);
  G4bool electronPol3 = !(pol3 == G4StokesVector::ZERO);

  G4double phi = 0.;
  // polarization independent part
  phi += fPhi0;

  if(gammaPol2)
  {
    // part depending on the polarization of the final photon
    phi += fPhi2 * pol2;
  }

  if(electronPol3)
  {
    // part depending on the polarization of the final electron
    phi += fPhi3 * pol3;
  }

  // return cross section.
  return phi;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4PolarizedComptonXS::TotalXSection(G4double /*xmin*/,
                                             G4double /*xmax*/, G4double k0,
                                             const G4StokesVector& pol0,
                                             const G4StokesVector& pol1)
{
  G4double k1 = 1. + 2. * k0;

  G4double unit = fZ * CLHEP::pi * CLHEP::classic_electr_radius *
                  CLHEP::classic_electr_radius;

  G4double pre = unit / (sqr(k0) * sqr(1. + 2. * k0));

  G4double xs_0 = ((k0 - 2.) * k0 - 2.) * sqr(k1) * std::log(k1) +
                  2. * k0 * (k0 * (k0 + 1.) * (k0 + 8.) + 2.);
  G4double xs_pol = (k0 + 1.) * sqr(k1) * std::log(k1) -
                    2. * k0 * (5. * sqr(k0) + 4. * k0 + 1.);

  return pre * (xs_0 / k0 + pol0.p3() * pol1.z() * xs_pol);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4StokesVector G4PolarizedComptonXS::GetPol2()
{
  // Note, mean polarization can not contain correlation effects.
  return G4StokesVector(1. / fPhi0 * fPhi2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4StokesVector G4PolarizedComptonXS::GetPol3()
{
  // Note, mean polarization can not contain correlation effects.
  return G4StokesVector(1. / fPhi0 * fPhi3);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4PolarizedComptonXS::DefineCoefficients(const G4StokesVector& pol0,
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
