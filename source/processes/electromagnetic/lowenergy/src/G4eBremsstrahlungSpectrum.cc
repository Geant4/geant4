//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4eBremsstrahlungSpectrum.cc,v 1.7 2002-05-30 17:53:09 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4eBremsstrahlungSpectrum
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
//
// -------------------------------------------------------------------

#include "G4eBremsstrahlungSpectrum.hh"
#include "G4BremsstrahlungParameters.hh"
#include "Randomize.hh"


G4eBremsstrahlungSpectrum::G4eBremsstrahlungSpectrum():
  G4VEnergySpectrum(),
  lowestE(0.1*eV),
  length(15)
{
  theBRparam = new G4BremsstrahlungParameters();
  xp.clear();
  for(size_t i=0; i<length; i++) {
    G4double x = 0.1*((G4double)i);
    if(i == 0)  x = 0.01;
    if(i == 10) x = 0.95;
    if(i == 11) x = 0.97;
    if(i == 12) x = 0.99;
    if(i == 13) x = 0.995;
    if(i == 14) x = 1.0;
    xp.push_back(x);
  }
  verbose = 0;
}


G4eBremsstrahlungSpectrum::~G4eBremsstrahlungSpectrum() 
{
  delete theBRparam;
}


G4double G4eBremsstrahlungSpectrum::Probability(G4int Z, 
                                                G4double tmin, 
                                                G4double tmax, 
                                                G4double e,
                                                G4int,
                                       const G4ParticleDefinition*) const
{
  G4double tm = G4std::min(tmax, e);
  G4double t0 = G4std::max(tmin, lowestE);
  if(t0 >= tm) return 0.0;

  t0 /= e;
  tm /= e;

  G4double z = lowestE/e;
  G4DataVector p;

  // Access parameters
  for (size_t i=0; i<=length; i++) {
    p.push_back(theBRparam->Parameter(i, Z, e));
  }
 
  G4double x  = IntSpectrum(t0, tm, p);
  G4double y  = IntSpectrum(z, 1.0, p); 


  if(1 < verbose) {
    G4cout << "tcut(MeV)= " << tmin/MeV 
           << "; tMax(MeV)= " << tmax/MeV 
           << "; t0= " << t0 
           << "; tm= " << tm 
           << "; xp[0]= " << xp[0]
           << "; z= " << z 
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


G4double G4eBremsstrahlungSpectrum::AverageEnergy(G4int Z,
                                                  G4double tmin, 
                                                  G4double tmax, 
                                                  G4double e,
                                                  G4int,
                                         const G4ParticleDefinition*) const
{
  G4double tm = G4std::min(tmax, e);
  G4double t0 = G4std::max(tmin, lowestE);
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

  G4double c  = sqrt(theBRparam->ParameterC(Z));
  G4double f  = Function(z0, p);
  x += 0.5*f*z0*(z0 - c*atan(z0/c));

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


G4double G4eBremsstrahlungSpectrum::SampleEnergy(G4int Z,
                                                 G4double tmin, 
                                                 G4double tmax, 
                                                 G4double e,
                                                 G4int,
                                        const G4ParticleDefinition*) const
{
  G4double tm = G4std::min(tmax, e);
  G4double t0 = G4std::max(tmin, lowestE);
  if(t0 >= tm) return 0.0;

  t0 /= e;
  tm /= e;

  G4DataVector p;

  for (size_t i=0; i<=length; i++) {
    p.push_back(theBRparam->Parameter(i, Z, e));
  }
  G4double amaj = G4std::max(p[15], 1. - (p[1] - p[0])/9.);

  G4double amax = log(tm);
  G4double amin = log(t0);
  G4double tgam, q, fun;

  do {
    G4double x = amin + G4UniformRand()*(amax - amin);
    tgam = exp(x);
    fun = Function(tgam, p);

    if(fun > amaj) {
          G4cout << "WARNING in G4eBremsstrahlungSpectrum::SampleEnergy:" 
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

G4double G4eBremsstrahlungSpectrum::IntSpectrum(G4double xMin, 
                                                G4double xMax,
                                          const G4DataVector& p) const
{
  G4double x1 = G4std::min(xMin, xp[0]);
  G4double x2 = G4std::min(xMax, xp[0]);
  G4double sum = 0.0;
  
  if(x1 < x2) {
    G4double k = (p[1] - p[0])/0.09;
    sum += (1. - k*xp[0])*log(x2/x1) + k*(x2 - x1);
  }

  for (size_t i=0; i<length-1; i++) {
    x1 = G4std::max(xMin, xp[i]);
    x2 = G4std::min(xMax, xp[i+1]);
    if(x1 < x2) {
      G4double z1 = p[i];
      G4double z2 = p[i+1];
      sum += z2 - z1 + log(x2/x1)*(z1*x2 - z2*x1)/(x2 - x1);
    }
  }
  if(sum < 0.0) sum = 0.0;
  return sum;
}

G4double G4eBremsstrahlungSpectrum::AverageValue(G4double xMin, 
                                                 G4double xMax,
			                   const G4DataVector& p) const
{
  G4double x1 = G4std::min(xMin, xp[0]);
  G4double x2 = G4std::min(xMax, xp[0]);
  G4double z1 = x1;
  G4double z2 = x2;
  G4double sum = 0.0;
  
  if(x1 < x2) {
    G4double k = (p[1] - p[0])/0.09;
    sum += (z2 - z1)*(1. - k*xp[0]);
    z1 *= x1;
    z2 *= x2;
    sum += 0.5*k*(z2 - z1);
  }

  for (size_t i=0; i<length-1; i++) {
    x1 = G4std::max(xMin, xp[i]);
    x2 = G4std::min(xMax, xp[i+1]);
    if(x1 < x2) {
      z1 = p[i];
      z2 = p[i+1];
      sum += 0.5*(z2 - z1)*(x2 + x1) + z1*x2 - z2*x1;
    }
  }
  if(sum < 0.0) sum = 0.0;
  return sum;
} 
  
G4double G4eBremsstrahlungSpectrum::Function(G4double x, 
                                       const G4DataVector& p) const
{
  G4double f = 0.0;
  
  if(x <= xp[0]) {
    f = 1. + (p[1] - p[0])*(x - xp[0])/0.09;

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

void G4eBremsstrahlungSpectrum::PrintData() const
{ theBRparam->PrintData(); }

G4double G4eBremsstrahlungSpectrum::Excitation(G4int Z, G4double kineticEnergy) const
{
  return 0.0;
}

G4double G4eBremsstrahlungSpectrum::MaxEnergyOfSecondaries(G4double kineticEnergy,
							   G4int Z,
							   const G4ParticleDefinition*) const
{ 
  return kineticEnergy; 
}
