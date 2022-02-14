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
// File name:     G4StokesVector
//
// Author:        Andreas Schaelicke
//
// Class Description:
//   Provides Stokes vector representation employed in polarized processes.

#include "G4StokesVector.hh"

#include "G4PolarizationHelper.hh"
#include "Randomize.hh"

const G4StokesVector G4StokesVector::ZERO =
  G4StokesVector(G4ThreeVector(0., 0., 0.));
const G4StokesVector G4StokesVector::P1 =
  G4StokesVector(G4ThreeVector(1., 0., 0.));
const G4StokesVector G4StokesVector::P2 =
  G4StokesVector(G4ThreeVector(0., 1., 0.));
const G4StokesVector G4StokesVector::P3 =
  G4StokesVector(G4ThreeVector(0., 0., 1.));
const G4StokesVector G4StokesVector::M1 =
  G4StokesVector(G4ThreeVector(-1., 0., 0.));
const G4StokesVector G4StokesVector::M2 =
  G4StokesVector(G4ThreeVector(0., -1., 0.));
const G4StokesVector G4StokesVector::M3 =
  G4StokesVector(G4ThreeVector(0., 0., -1.));

G4StokesVector::G4StokesVector()
  : G4ThreeVector()
  , fIsPhoton(false)
{}

G4StokesVector::G4StokesVector(const G4ThreeVector& v)
  : G4ThreeVector(v)
  , fIsPhoton(false)
{}

G4bool G4StokesVector::IsZero() const { return *this == ZERO; }

void G4StokesVector::RotateAz(G4ThreeVector nInteractionFrame,
                              G4ThreeVector particleDirection)
{
  G4ThreeVector yParticleFrame =
    G4PolarizationHelper::GetParticleFrameY(particleDirection);

  G4double cosphi = yParticleFrame * nInteractionFrame;
  if(cosphi > (1. + 1.e-8) || cosphi < (-1. - 1.e-8))
  {
    G4ExceptionDescription ed;
    ed << " warning G4StokesVector::RotateAz  cosphi>1 or cosphi<-1\n"
       << " cosphi=" << cosphi << "\n"
       << " zAxis=" << particleDirection << " (" << particleDirection.mag()
       << ")\n"
       << " yAxis=" << yParticleFrame << " (" << yParticleFrame.mag() << ")\n"
       << " nAxis=" << nInteractionFrame << " (" << nInteractionFrame.mag()
       << ")\n";
    G4Exception("G4StokesVector::RotateAz", "pol030", JustWarning, ed);
  }
  if(cosphi > 1.)
    cosphi = 1.;
  else if(cosphi < -1.)
    cosphi = -1.;

  G4double hel =
    (yParticleFrame.cross(nInteractionFrame) * particleDirection) > 0. ? 1.
                                                                       : -1.;

  G4double sinphi = hel * std::sqrt(1. - cosphi * cosphi);

  RotateAz(cosphi, sinphi);
}

void G4StokesVector::InvRotateAz(G4ThreeVector nInteractionFrame,
                                 G4ThreeVector particleDirection)
{
  // note if incoming particle is on z-axis,
  // we might encounter some nummerical problems, since
  // nInteratonFrame and yParticleFrame are actually (almost) the same momentum
  // and the normalization is only good to 10^-12 !

  G4ThreeVector yParticleFrame =
    G4PolarizationHelper::GetParticleFrameY(particleDirection);
  G4double cosphi = yParticleFrame * nInteractionFrame;

  if(cosphi > 1. + 1.e-8 || cosphi < -1. - 1.e-8)
  {
    G4ExceptionDescription ed;
    ed << " warning G4StokesVector::RotateAz  cosphi>1 or cosphi<-1\n";
    G4Exception("G4StokesVector::InvRotateAz", "pol030", JustWarning, ed);
  }
  if(cosphi > 1.)
    cosphi = 1.;
  else if(cosphi < -1.)
    cosphi = -1.;

  // check sign once more!
  G4double hel =
    (yParticleFrame.cross(nInteractionFrame) * particleDirection) > 0. ? 1.
                                                                       : -1.;
  G4double sinphi = hel * std::sqrt(std::fabs(1. - cosphi * cosphi));
  RotateAz(cosphi, -sinphi);
}

void G4StokesVector::RotateAz(G4double cosphi, G4double sinphi)
{
  if(!fIsPhoton)
  {
    G4double xsi1 = cosphi * p1() + sinphi * p2();
    G4double xsi2 = -sinphi * p1() + cosphi * p2();
    setX(xsi1);
    setY(xsi2);
    return;
  }

  G4double sin2phi = 2. * cosphi * sinphi;
  G4double cos2phi = cosphi * cosphi - sinphi * sinphi;

  G4double xsi1 = cos2phi * p1() + sin2phi * p2();
  G4double xsi2 = -sin2phi * p1() + cos2phi * p2();
  setX(xsi1);
  setY(xsi2);
}

G4double G4StokesVector::GetBeta()
{
  G4double bet = getPhi();
  if(fIsPhoton)
  {
    bet *= 0.5;
  }
  return bet;
}

void G4StokesVector::DiceUniform()
{
  G4double costheta = 2. * G4UniformRand() - 1.;
  G4double sintheta = std::sqrt(1. - costheta * costheta);
  G4double aphi     = 2. * CLHEP::pi * G4UniformRand();
  setX(std::sin(aphi) * sintheta);
  setY(std::cos(aphi) * sintheta);
  setZ(costheta);
}

void G4StokesVector::DiceP1()
{
  if(G4UniformRand() > 0.5)
    setX(1.);
  else
    setX(-1.);
  setY(0.);
  setZ(0.);
}

void G4StokesVector::DiceP2()
{
  setX(0.);
  if(G4UniformRand() > 0.5)
    setY(1.);
  else
    setY(-1.);
  setZ(0.);
}

void G4StokesVector::DiceP3()
{
  setX(0.);
  setY(0.);
  if(G4UniformRand() > 0.5)
    setZ(1.);
  else
    setZ(-1.);
}

void G4StokesVector::FlipP3() { setZ(-z()); }

G4ThreeVector G4StokesVector::PolError(const G4StokesVector& sum2, long n)
{
  // delta x = sqrt[ ( <x^2> - <x>^2 )/(n-1) ]
  G4ThreeVector mean   = (1. / n) * G4ThreeVector(*this);
  G4ThreeVector polsqr = G4StokesVector(mean).PolSqr();
  G4ThreeVector result =
    G4StokesVector((1. / (n - 1.) * ((1. / n) * sum2 - polsqr))).PolSqrt();
  return result;
}

G4ThreeVector G4StokesVector::PolDiv(const G4StokesVector& b)
{
  return G4ThreeVector(b.x() != 0. ? x() / b.x() : 11111.,
                       b.y() != 0. ? y() / b.y() : 11111.,
                       b.z() != 0. ? z() / b.z() : 11111.);
}
