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
// $Id: G4eBremsstrahlungSpectrum.cc,v 1.2 2001-11-15 11:33:33 vnivanch Exp $
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
//
// -------------------------------------------------------------------

#include "G4eBremsstrahlungSpectrum.hh"
#include "G4BremsstrahlungParameters.hh"
#include "Randomize.hh"


G4eBremsstrahlungSpectrum::G4eBremsstrahlungSpectrum():
  G4VEnergySpectrum(),
  lowestE(0.1*eV)
{
  theBRparam = new G4BremsstrahlungParameters();
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

  G4double a  = theBRparam->ParameterA(Z, e);
  G4double b  = 0.5*theBRparam->ParameterB(Z, e);

  t0 /= e;
  tm /= e;
 
  G4double z = lowestE/e;

  G4double x  = log(tm/t0) + (tm - t0)*(a + b*(tm + t0));
  G4double y  = log(1./z)  + (1. - z) *(a + b*(1. + z));     

  if(y > 0.0) x /= y;
  else        x  = 0.0;
  if(x < 0.0) x  = 0.0;

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

  G4double a  = theBRparam->ParameterA(Z, e);
  G4double b  = theBRparam->ParameterB(Z, e);
  G4double c  = theBRparam->ParameterC(Z);

  t0 /= e;
  tm /= e;

  G4double z = lowestE/e;

  G4double x  = (tm - t0)*(1.0 + 0.5*a*(tm + t0) 
                                 + b*(tm*tm +t0*tm + t0*t0)/3.0 );
  G4double z2 = z*z/c;
  x += (1. + z2)*(1. + a*z + b*z*z)*0.5*c*(z2 - log(1. + z2))/z2;

  x *= e; 

  G4double y  = log(1./z)  + (1. - z) *(a + b*(1. + z));     

  if(y > 0.0) x /= y;
  else        x  = 0.0;
  if(x < 0.0) x  = 0.0;

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

  G4double a  = theBRparam->ParameterA(Z, e);
  G4double b  = theBRparam->ParameterB(Z, e);

  t0 /= e;
  tm /= e;

  G4double amaj = G4std::max(1. + a*t0 + b*t0*t0, 1. + a*tm + b*tm*tm);
  G4double t    = -a*0.5/b;
  if(t > t0 && t < tm) amaj = G4std::max(amaj, 1.0 + 0.5*t*a);
  G4double amax = log(tm);
  G4double amin = log(t0);
  G4double tgam, q;
  do {
    G4double x = amin + G4UniformRand()*(amax - amin);
    tgam = exp(x);
    q = amaj * G4UniformRand();
  } while (q > 1. + a*tgam + b*tgam*tgam);

  tgam *= e;

  return tgam; 
}

void G4eBremsstrahlungSpectrum::PrintData() const
{ theBRparam->PrintData(); }
