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
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4BremsstrahlungElectronSpectrum
//
// Author:        V.Ivanchenko (Vladimir.Ivanchenko@cern.ch)
// 
// Creation date: 29 September 2001
//
// Modifications: 
// 10 Oct 2001  MGP  Revision to improve code quality and consistency with design
//
// -------------------------------------------------------------------

#include "G4BremsstrahlungElectronSpectrum.hh"
#include "G4BremsstrahlungParameters.hh"
#include "Randomize.hh"


G4BremsstrahlungElectronSpectrum::G4BremsstrahlungElectronSpectrum():
  G4VEnergySpectrum(),
  lowestE(0.1*eV)
{
  theBRparam = new G4BremsstrahlungParameters();
}


G4BremsstrahlungElectronSpectrum::~G4BremsstrahlungElectronSpectrum() 
{
  delete theBRparam;
}


G4double G4BremsstrahlungElectronSpectrum::Probability(G4int Z, 
                                                       G4double tmin, 
                                                       G4double tmax, 
                                                       G4double e,
                                                       G4int,
                                       const G4ParticleDefinition*) const
{
  G4double tm = G4std::min(tmax, e);
  G4double t0 = G4std::max(tmin, lowestE);
  if(tmin >= tm) return 0.0;

  G4double a  = theBRparam->ParameterA(Z, e);
  G4double b  = theBRparam->ParameterB(Z, e);

  G4double x  = a*log(tm/t0) + b*(tm - t0);
  G4double y  = a*log(e/lowestE)   + b*(e - lowestE);

  if(y > 0.0) x /= y;
  else        x  = 0.0;
  if(x < 0.0) x  = 0.0;

  return x; 
}


G4double G4BremsstrahlungElectronSpectrum::AverageEnergy(G4int Z,
                                                         G4double tmin, 
                                                         G4double tmax, 
                                                         G4double e,
                                                         G4int,
                                         const G4ParticleDefinition*) const
{
  G4double tm = G4std::min(tmax, e);
  G4double t0 = G4std::max(tmin, lowestE);
  if(tmin >= tm) return 0.0;

  G4double a  = theBRparam->ParameterA(Z, e);
  G4double b  = theBRparam->ParameterB(Z, e);

  G4double x  = (tm - t0)*(a + 0.5*b*(tm + t0));
  G4double y  = a*log(e/lowestE)   + b*(e - lowestE);

  if(y > 0.0) x /= y;
  else        x  = 0.0;
  if(x < 0.0) x  = 0.0;

  return x; 
}


G4double G4BremsstrahlungElectronSpectrum::SampleEnergy(G4int Z,
                                                        G4double tmin, 
                                                        G4double tmax, 
                                                        G4double e,
                                                        G4int,
                                        const G4ParticleDefinition*) const
{
  G4double tm = G4std::min(tmax, e);
  G4double t0 = G4std::max(tmin, lowestE);
  if(tmin >= tm) return 0.0;

  G4double a  = theBRparam->ParameterA(Z, e);
  G4double b  = theBRparam->ParameterB(Z, e);

  // Sample gamma energy
  // The emitted gamma energy is from EEDL data fitted with A/E+B function.
  // Original formula A/E+B+C*E and sampling methods are reported by 
  // J.Stepanek formula has been modified by A. Forti and S. Giani. 
  G4double amaj = a + b*tm;
  G4double amax = log(tm);
  G4double amin = log(t0);
  G4double tgam, q;
  do {
    G4double x = amin + G4UniformRand()*(amax - amin);
    tgam = exp(x);
    q = amaj * G4UniformRand();
  } while (q >= a + b*tgam);

  return tgam; 
}





