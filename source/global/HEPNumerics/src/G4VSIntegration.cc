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
  nbin = std::max(nbin, 6);
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
      // shifted energy limit between the 1st and 2nd areas
      } else if (y > fP1) {
	fP2 = 0.0;
	fE1 = x;
	fP1 = y;
	fE2 = emax;
	// definition of the 3d area
      } else if (0.0 == fP2 && y < fFactor1*fP1) {
        fE2 = x;
        fP2 = y;
	// extra maximum inside the 3d area
        // shifted energy limit between the 2nd and 3d areas
      } else if (0.0 < fP2 && y > fP2) {
	if (y > fP1) { fP1 = y; }
	fE2 = x;
	fP2 = y;
      }
    }
    
    G4double del = (y + problast)*edelta*0.5;
    res += del;
    
    // end of the loop condition
    if ((del < fAcc*res && 0 < fP2) || endpoint) { break; }
    problast = y;

    // smart next step definition
    if (del != res) {
      if (del > 0.8*res && 0.7*edelta > fMinDelta) { 
	edelta *= 0.7;
      } else if (del < 0.1*res && 1.5*edelta < fMaxDelta) { 
	edelta *= 1.5;
      }
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

  // if 3d region is considered it is subdivided on 2 parts
  // so sampling may be performed in 4 energy intervals
  G4double p3 = 0.0;
  G4double p4 = 0.0;
  G4double E3 = fEmax;
  G4double Q0 = fPmax*fFactor2;
  G4double Q1 = fP1*fFactor2;
  G4double Q2 = fP2*fFactor2;
  G4double Q3 = 0.0;
  
  // for some distributions it may happens that there is a local maximum
  // closed to maximal energy
  G4double Q4 = ProbabilityDensityFunction(fEmax - 0.02*(fEmax - fEmin))*fFactor2;
  Q0 = std::max(Q0, Q4);
  if (Q1 > 0.0) { Q1 = std::max(Q1, Q4); }

  // integral under the 1st and the 2nd areas 
  G4double p1 = (fE1 - fEmin)*Q0;
  G4double p2 = (fE2 - fE1)*Q1;
  
  // if p2 is very small the 2nd area should not be considered
  if (p2 < 1.e-8*p1) {
    p2 = 0.0;
    p1 = (fE2 - fEmin)*Q0;
    fE1 = fE2;
  }

  // 3d area may be considered
  if (Q2 > 0.0 && fE2 < fEmax) {
    E3 = fEmax - 0.5*(fEmax - fE2);
    Q3 = ProbabilityDensityFunction(E3)*fFactor2;
    Q2 = std::max(Q2, Q3);
    Q3 = std::max(Q3, Q4);
    p3 = (E3 - fE2)*Q2;
    p4 = (fEmax - E3)*Q3;
  }
  // sampling in 4 areas, probabilities may be zero except p1
  G4double sum = p1 + p2 + p3 + p4;
  G4double del1 = (fE1 - fEmin)/p1;
  G4double del2 = (p2 > 0.0) ? (fE2 - fE1)/p2 : 0.0;
  G4double del3 = (p3 > 0.0) ? (E3 - fE2)/p3 : 0.0;
  G4double del4 = (p4 > 0.0) ? (fEmax - E3)/p4 : 0.0;

  CLHEP::HepRandomEngine* rndm = G4Random::getTheEngine();
  const G4int nmax = 100000;
  G4double e, gmax, gg;
  for (G4int n=0; n < nmax; ++n) {
    G4double q = rndm->flat();
    G4double p = sum*q;
    G4int idx = 0;
    if (p <= p1) {
      gmax = Q0;
      e = del1*p + fEmin;
    } else if (p <= p1 + p2) {
      gmax = Q1;
      e = del2*(p - p1) + fE1;
      idx = 1;
    } else if (p <= p1 + p2 + p3) {
      gmax = Q2;
      e = del3*(p - p1 - p2) + fE2;
      idx = 2;
    } else {
      gmax = Q3;
      e = del4*(p - p1 - p2 - p3) + E3;
      idx = 3;
    }
    gg = ProbabilityDensityFunction(e);
    if ((gg > gmax || n >= nmax) && fVerbose > 0) {
      ++fnWarn;
      if (fnWarn < fWarnLimit) {
	G4cout << "### G4VSIntegration::SampleValue() for " << ModelName()
	       << " in area=" << idx << " n=" << n << " gg/gmax=" << gg/gmax
	       << " prob=" << gg << " gmax=" << gmax << G4endl; 
	G4cout << "    E=" << e << " Emin=" << fEmin << " Emax=" << fEmax
	       << " E1=" << fE1 << " E2=" << fE2 << " E3=" << E3
	       << " F0=" << Q0 << " F1=" << Q1 << " F2=" << Q2 << " F3=" << Q3
	       << " F4=" << Q4 << G4endl;
      }
    }
    if (gmax*rndm->flat() <= gg) {
#ifdef G4VERBOSE
      if (fVerbose > 1) {
	G4cout << "### G4VSIntegration::SampleValue for " << ModelName()
	       << " E=" << e << " Ntry=" << n
	       << " Emin=" << fEmin << " Emax=" << fEmax << G4endl;
      }
#endif
      return e;
    }
  }
  // if sampling not converged, then sample uniformaly in the 1st energy region
  e = fEmin + rndm->flat()*(fE1 - fEmin);
#ifdef G4VERBOSE
  if (fVerbose > 1) {
    G4cout << "### G4VSIntegration::SampleValue for " << ModelName()
           << " E=" << e << " Ntry=" << nmax
           << " Emin=" << fEmin << " Emax=" << fEmax << G4endl;
  }
#endif
  return e;
}

const G4String& G4VSIntegration::ModelName() const
{
  return dummy;
}
