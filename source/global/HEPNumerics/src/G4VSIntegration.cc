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
// G4VSIntegration
//
// Created 03.03.2025 V.Ivanchenko
//
// --------------------------------------------------------------------

#include "G4VSIntegration.hh"
#include "Randomize.hh"
#include "G4Log.hh"

void
G4VSIntegration::InitialiseIntegrator(G4double acc, G4double f1, G4double f2,
				      G4double de, G4double dmin, G4double dmax)
{
  if (acc > 0.0) { fAcc = acc; }
  if (f1 > 0.0 && f1 < 1.0) { fFactor1 = f1; }
  if (f2 > 1.0 && f2 < 5.0) { fFactor2 = f2; }
  if (de > 0.0) { fDelta = de; }
  fMinDelta = (dmin <= de && dmin > 0.0) ? dmin : de;
  fMaxDelta = (dmax > de) ? dmax : de;
  if (fVerbose > 2) {
    G4cout << "### G4VSIntegration::InitialiseIntegrator: "
           << "fAcc=" << fAcc << " fFact1=" << fFactor1
	   << " fFact2=" << fFactor2 << " dE=" << fDelta
	   << " dEmin=" << fMinDelta << " dEmax=" << fMaxDelta << G4endl;
  }
}

G4double
G4VSIntegration::ComputeIntegral(const G4double emin, const G4double  emax)
{
  G4double res = 0.0;
  if (emin >= emax) { return res; }
  fEmin = emin;
  fEmax = emax;

  // preparing smart binning
  G4int nbin = G4lrint((emax - emin)/fDelta) + 1;
  nbin = std::max(nbin, 4);
  G4double edelta = (emax - emin)/static_cast<G4double>(nbin);
  nbin += nbin;

  // prepare integration
  G4double x(emin), y(0.0);
  fPmax = ProbabilityDensityFunction(x);
  G4double problast = fPmax;
#ifdef G4VERBOSE
  if (fVerbose > 1) {
    G4cout << "### G4VSIntegration::ComputeIntegral: "
           << "Pmax=" << fPmax << " Emin=" << emin
           << " Emax=" << emax << " dE=" << edelta << " nbin=" << nbin 
           << G4endl;
  }
#endif

  fE1 = fE2 = emax;
  fP1 = fP2 = 0.0;
  G4bool endpoint = false;
  x += edelta;
  for (G4int i=0; i<=nbin; ++i) {
    // the last point - may be earlier due to dynamic interval
    if (x >= emax) {
      edelta += emax - x;
      x = emax;
      endpoint = true;
    }
    y = ProbabilityDensityFunction(x);
#ifdef G4VERBOSE
    if (fVerbose > 2) { 
      G4cout << "    " << i << ".  E=" << x << "  prob=" << y
             << " Edel=" << edelta << G4endl;
    } 
#endif
    if (y >= fPmax) {
      fPmax = y;
      // for the case of 2nd maximum 
      fP1 = fP2 = 0.0;
      fE1 = fE2 = emax;
    } else if (!endpoint) {
      if (0.0 == fP1) {
	// definition of the 2nd area
        if (y < fFactor1*fPmax) {
          fE1 = x;
          fP1 = y;
        }
      // for the case of 2nd maximum in the 2nd area
      } else if (y > fP1) {
	fP2 = 0.0;
	fP1 = y;
	fE2 = emax;
	// definition of the 3d area
      } else if (0.0 == fP2 && y < fFactor1*fP1) {
        fE2 = x;
        fP2 = y;
	// extra maximum inside the 3d area
      } else if (0.0 < fP2 && y > fP2) {
        fP2 = y;
      }
    }
    
    G4double del = (y + problast)*edelta*0.5;
    res += del;
    // end of the loop 
    if (del < fAcc*res || endpoint) { break; }
    problast = y;

    // smart next step definition
    if (del != res && del > 0.8*res && 0.7*edelta > fMinDelta) { 
      edelta *= 0.7;
    } else if (del < 0.1*res && 1.5*edelta < fMaxDelta) { 
      edelta *= 1.5;
    }
    x += edelta;
  }
#ifdef G4VERBOSE
  if (fVerbose > 1) {
    G4cout << "### G4VSIntegration::ComputeIntegral: "
           << "I=" << res << " E1=" << fE1 << " E2=" << fE2
           << " Pmax=" << fPmax << " P1=" << fP1 << " P2=" << fP2 
           << G4endl;
  }
#endif
  return res;
}

G4double G4VSIntegration::SampleValue()
{
  // should never happen
  if (fEmin >= fEmax) { return fEmin; }

  fPmax *= fFactor2;

  // two regions with flat and one with exponential majorant 
  G4double b = 1.0;
  G4double p3 = 0.0;
  G4double Q3 = 0.0;

  // 2d and 3d areas may be considered
  if (fP2 > 0.0 && fE2 < fEmax) {
    Q3 = 2*ProbabilityDensityFunction(fEmax - 0.5*(fEmax - fE2));
    // exclude 3d area from sampling
    if (4*Q3 > fP2 || Q3 <= 0.0) {
      fE2 = fEmax;
      // 3d area is considered
    } else {
      b = 2*G4Log(fP2/Q3)/(fEmax - fE2);
      p3 = (fP2 - Q3)/b;
    }
  }
  G4double p1 = (fE1 - fEmin)*fPmax;
  G4double p2 = (fE2 - fE1)*fP1;
  if (p2 < 1.e-8*p1) {
    p2 = 0.0;
    p1 = (fE2 - fEmin)*fPmax;
    fE1 = fE2;
  }
  G4double sum = p1 + p2 + p3;
  G4double del1 = (fE1 - fEmin)/p1;
  G4double del2 = (p2 > 0.0) ? (fE2 - fE1)/p2 : 0.0;

  CLHEP::HepRandomEngine* rndm = G4Random::getTheEngine();
  const G4int nmax = 1000;
  G4double e, gmax, gg;
  G4int n = 0;
  do {
    ++n;
    G4double q = rndm->flat();
    G4double p = sum*q;
    G4int idx = 0;
    if (p <= p1) {
      gmax = fPmax;
      e = del1*p + fEmin;
    } else if (p <= p1 + p2) {
      gmax = fP1;
      e = del2*(p - p1) + fE1;
      idx = 1;
    } else {
      G4double x = 1.0 - rndm->flat()*(1.0 - Q3/fP2);
      e = fE2 - G4Log(x)/b;
      gmax = fP2*x;
      idx = 2;
    }
    gg = ProbabilityDensityFunction(e);
    if ((gg > gmax || n >= nmax) && fVerbose > 0) {
      ++fnWarn;
      if (fnWarn < fWarnLimit) {
	G4cout << "### G4VSIntegration::SampleValue() for " << ModelName()
	       << " in area=" << idx << " n=" << n << " gg/gmax=" << gg/gmax
	       << " prob=" << gg << " gmax=" << gmax << G4endl; 
	G4cout << "    E=" << e << " Emin=" << fEmin << " Emax=" << fEmax
	       << " E1=" << fE1 << " E2=" << fE2 << " Fmax=" << fPmax
	       << " F1=" << fP1 << G4endl;
      }
    }
  } while(gmax*rndm->flat() > gg && n < nmax);
#ifdef G4VERBOSE
  if (fVerbose > 1) {
    G4cout << "### G4VSIntegration::SampleValue for " << ModelName()
           << " E=" << e << " Ntry=" << n
           << " Emin=" << fEmin << " Emax=" << fEmax << G4endl;
  }
#endif
  return e;
}

const G4String& G4VSIntegration::ModelName() const
{
  return dummy;
}
