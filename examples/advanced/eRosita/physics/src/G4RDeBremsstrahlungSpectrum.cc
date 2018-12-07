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
// File name:     G4RDeBremsstrahlungSpectrum
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
//
// Creation date: 29 September 2001
//
// Modifications:
// 10.10.01  MGP  Revision to improve code quality and consistency with design
// 15.11.01  VI   Update spectrum model Bethe-Haitler spectrum at high energy
// 30.05.02  VI   Update interpolation between 2 last energy points in the
//                parametrisation
// 21.02.03  V.Ivanchenko    Energy bins are defined in the constructor
// 28.02.03  V.Ivanchenko    Filename is defined in the constructor
//
// -------------------------------------------------------------------

#include "G4RDeBremsstrahlungSpectrum.hh"
#include "G4RDBremsstrahlungParameters.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"


G4RDeBremsstrahlungSpectrum::G4RDeBremsstrahlungSpectrum(const G4DataVector& bins,
  const G4String& name):G4RDVEnergySpectrum(),
  lowestE(0.1*eV),
  xp(bins)
{
  length = xp.size();
  theBRparam = new G4RDBremsstrahlungParameters(name,length+1);

  verbose = 0;
}


G4RDeBremsstrahlungSpectrum::~G4RDeBremsstrahlungSpectrum()
{
  delete theBRparam;
}


G4double G4RDeBremsstrahlungSpectrum::Probability(G4int Z,
                                                G4double tmin,
                                                G4double tmax,
                                                G4double e,
                                                G4int,
                                       const G4ParticleDefinition*) const
{
  G4double tm = std::min(tmax, e);
  G4double t0 = std::max(tmin, lowestE);
  if(t0 >= tm) return 0.0;

  t0 /= e;
  tm /= e;

  G4double z0 = lowestE/e;
  G4DataVector p;

  // Access parameters
  for (size_t i=0; i<=length; i++) {
    p.push_back(theBRparam->Parameter(i, Z, e));
  }

  G4double x  = IntSpectrum(t0, tm, p);
  G4double y  = IntSpectrum(z0, 1.0, p);


  if(1 < verbose) {
    G4cout << "tcut(MeV)= " << tmin/MeV
           << "; tMax(MeV)= " << tmax/MeV
           << "; t0= " << t0
           << "; tm= " << tm
           << "; xp[0]= " << xp[0]
           << "; z= " << z0
           << "; val= " << x
           << "; nor= " << y
           << G4endl;
  }
  p.clear();

  if(y > 0.0) x /= y;
  else        x  = 0.0;
  //  if(x < 0.0) x  = 0.0;

  return x;
}


G4double G4RDeBremsstrahlungSpectrum::AverageEnergy(G4int Z,
                                                  G4double tmin,
                                                  G4double tmax,
                                                  G4double e,
                                                  G4int,
						  const G4ParticleDefinition*) const
{
  G4double tm = std::min(tmax, e);
  G4double t0 = std::max(tmin, lowestE);
  if(t0 >= tm) return 0.0;

  t0 /= e;
  tm /= e;

  G4double z0 = lowestE/e;

  G4DataVector p;

  // Access parameters
  for (size_t i=0; i<=length; i++) {
      p.push_back(theBRparam->Parameter(i, Z, e));
  }

  G4double x  = AverageValue(t0, tm, p);
  G4double y  = IntSpectrum(z0, 1.0, p);

  // Add integrant over lowest energies
  G4double zmin = tmin/e;
  if(zmin < t0) {
    G4double c  = std::sqrt(theBRparam->ParameterC(Z));
    x += p[0]*(t0 - zmin - c*(std::atan(t0/c) - std::atan(zmin/c)));
  }
  x *= e;

  if(1 < verbose) {
    G4cout << "tcut(MeV)= " << tmin/MeV
           << "; tMax(MeV)= " << tmax/MeV
           << "; e(MeV)= " << e/MeV
           << "; t0= " << t0
           << "; tm= " << tm
           << "; y= " << y
           << "; x= " << x
           << G4endl;
  }
  p.clear();

  if(y > 0.0) x /= y;
  else        x  = 0.0;
  //  if(x < 0.0) x  = 0.0;

  return x;
}


G4double G4RDeBremsstrahlungSpectrum::SampleEnergy(G4int Z,
                                                 G4double tmin,
                                                 G4double tmax,
                                                 G4double e,
                                                 G4int,
						 const G4ParticleDefinition*) const
{
  G4double tm = std::min(tmax, e);
  G4double t0 = std::max(tmin, lowestE);
  if(t0 >= tm) return 0.0;

  t0 /= e;
  tm /= e;

  G4DataVector p;

  for (size_t i=0; i<=length; i++) {
    p.push_back(theBRparam->Parameter(i, Z, e));
  }
  G4double amaj = std::max(p[length], 1. - (p[1] - p[0])*xp[0]/(xp[1] - xp[0]) );

  G4double amax = std::log(tm);
  G4double amin = std::log(t0);
  G4double tgam, q, fun;

  do {
    G4double x = amin + G4UniformRand()*(amax - amin);
    tgam = std::exp(x);
    fun = Function(tgam, p);

    if(fun > amaj) {
          G4cout << "WARNING in G4RDeBremsstrahlungSpectrum::SampleEnergy:"
                 << " Majoranta " << amaj
                 << " < " << fun
                 << G4endl;
    }

    q = amaj * G4UniformRand();
  } while (q > fun);

  tgam *= e;

  p.clear();

  return tgam;
}

G4double G4RDeBremsstrahlungSpectrum::IntSpectrum(G4double xMin,
                                                G4double xMax,
						const G4DataVector& p) const
{
  G4double x1 = std::min(xMin, xp[0]);
  G4double x2 = std::min(xMax, xp[0]);
  G4double sum = 0.0;

  if(x1 < x2) {
    G4double k = (p[1] - p[0])/(xp[1] - xp[0]);
    sum += (1. - k*xp[0])*std::log(x2/x1) + k*(x2 - x1);
  }

  for (size_t i=0; i<length-1; i++) {
    x1 = std::max(xMin, xp[i]);
    x2 = std::min(xMax, xp[i+1]);
    if(x1 < x2) {
      G4double z1 = p[i];
      G4double z2 = p[i+1];
      sum += z2 - z1 + std::log(x2/x1)*(z1*x2 - z2*x1)/(x2 - x1);
    }
  }
  if(sum < 0.0) sum = 0.0;
  return sum;
}

G4double G4RDeBremsstrahlungSpectrum::AverageValue(G4double xMin,
                                                 G4double xMax,
						 const G4DataVector& p) const
{
  G4double x1 = std::min(xMin, xp[0]);
  G4double x2 = std::min(xMax, xp[0]);
  G4double z1 = x1;
  G4double z2 = x2;
  G4double sum = 0.0;

  if(x1 < x2) {
    G4double k = (p[1] - p[0])/(xp[1] - xp[0]);
    sum += (z2 - z1)*(1. - k*xp[0]);
    z1 *= x1;
    z2 *= x2;
    sum += 0.5*k*(z2 - z1);
  }

  for (size_t i=0; i<length-1; i++) {
    x1 = std::max(xMin, xp[i]);
    x2 = std::min(xMax, xp[i+1]);
    if(x1 < x2) {
      z1 = p[i];
      z2 = p[i+1];
      sum += 0.5*(z2 - z1)*(x2 + x1) + z1*x2 - z2*x1;
    }
  }
  if(sum < 0.0) sum = 0.0;
  return sum;
}

G4double G4RDeBremsstrahlungSpectrum::Function(G4double x,
					     const G4DataVector& p) const
{
  G4double f = 0.0;

  if(x <= xp[0]) {
    f = p[0] + (p[1] - p[0])*(x - xp[0])/(xp[1] - xp[0]);

  } else {

    for (size_t i=0; i<length-1; i++) {

      if(x <= xp[i+1]) {
        f = p[i] + (p[i+1] - p[i])*(x - xp[i])/(xp[i+1] - xp[i]);
        break;
      }
    }
  }
  return f;
}

void G4RDeBremsstrahlungSpectrum::PrintData() const
{ theBRparam->PrintData(); }

G4double G4RDeBremsstrahlungSpectrum::Excitation(G4int , G4double ) const
{
  return 0.0;
}

G4double G4RDeBremsstrahlungSpectrum::MaxEnergyOfSecondaries(G4double kineticEnergy,
							   G4int, // Z,
							   const G4ParticleDefinition*) const
{
  return kineticEnergy;
}
