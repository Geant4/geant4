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
// $Id: G4VEmissionProbability.cc 66241 2012-12-13 18:34:42Z gunter $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// Modifications:
// 28.10.2010 V.Ivanchenko defined members in constructor and cleaned up

#include "G4VEmissionProbability.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"
#include "Randomize.hh"

G4VEmissionProbability::G4VEmissionProbability(G4int Z, G4int A)
  :OPTxs(3),fVerbose(0),theZ(Z),theA(A), 
   LevelDensity(0.1),elimit(CLHEP::MeV),accuracy(0.02) 
{
  fG4pow = G4Pow::GetInstance();
  fPairCorr = G4PairingCorrection::GetInstance();
  length = nfilled = 0;
  emin = emax = eCoulomb = probmax = eprobmax = totProbability = 0.0;
}

G4VEmissionProbability::~G4VEmissionProbability() 
{}

void G4VEmissionProbability::Initialise()
{
  G4DeexPrecoParameters* param = G4NuclearLevelData::GetInstance()->GetParameters();
  OPTxs = param->GetDeexModelType();
  LevelDensity = param->GetLevelDensity();
}

void G4VEmissionProbability::ResetIntegrator(size_t nbin, G4double de, G4double eps)
{
  if(nbin > 0) {
    fProb.clear();
    fEner.clear();
    fEner.resize(nbin+1, 0.0);
    fProb.resize(nbin+1, 0.0);
    length = nbin;
  }
  if(de > 0.0) { elimit = de; }
  if(accuracy > 0.0) { accuracy = eps; }
}

G4double G4VEmissionProbability::ComputeProbability(G4double, G4double)
{
  return 0.0;
}

G4double 
G4VEmissionProbability::IntegrateProbability(G4double elow, G4double ehigh, G4double cb)
{
  totProbability = 0.0;
  if(elow >= ehigh) { return totProbability; }

  emin = elow;
  emax = ehigh;
  eCoulomb = cb;

  size_t nbin = (size_t)((emax - emin)/elimit) + 1;
  G4double edelta = elimit;
  if(nbin < 3) { 
    nbin = 3; 
    edelta = (emax - emin)/(G4double)nbin;
  }
  if(nbin > length) { 
    fEner.resize(nbin + 1); 
    fProb.resize(nbin + 1); 
    length = nbin;
  }

  G4double x(emin), del, y; 
  G4double problast = ComputeProbability(x, eCoulomb);
  eprobmax= emin;
  probmax = problast;
  fEner[0] = emin;
  fProb[0] = problast;
  size_t i(0);
  if(fVerbose > 1) {
    G4cout << "### G4VEmissionProbability::IntegrateProbability: " 
	   << " Emax= " << emax << " QB= " << cb << " nbin= " << nbin 
	   << G4endl;
    G4cout << "    0.  E= " << emin << "  prob= " << problast << G4endl;
  }
  static const G4double edeltamin = 0.2*CLHEP::MeV;
  static const G4double edeltamax = 2*CLHEP::MeV;
  G4bool peak = true;
  do {
    ++i;
    x += edelta;
    if(x > emax) { 
      x = emax; 
      edelta = emax - fEner[i-1];
    }
    y = ComputeProbability(x, eCoulomb);
    if(fVerbose > 1) { 
      G4cout << "    " << i << ".  E= " << x << "  prob= " << y 
	     << " Edel= " << edelta << G4endl;
    } 
    fEner[i] = x;
    fProb[i] = y;
    del = (y + problast)*edelta*0.5;
    totProbability += del;

    // smart step definition
    if(del < accuracy*totProbability) { break; }
    else if(del != totProbability && del > 0.8*totProbability 
	    && edelta > edeltamin) 
      { edelta *= 0.5; } 
    else if(del < 0.1*totProbability && edelta < edeltamax) 
      { edelta *= 2.0; } 
    if(y > probmax) {
      probmax = y;
      eprobmax = x;
    } else if(peak && y < probmax) { 
      edelta *= 2.0;
      peak = false;
    } 
    problast = y;
    // Loop checking, 10-Mar-2017, Vladimir Ivanchenko
  } while(i < nbin && x < emax);

  nfilled = i;
  if(fVerbose > 1) { G4cout << " Probability= " << totProbability << G4endl; }
  return totProbability;
}

G4double G4VEmissionProbability::SampleEnergy()
{
  static const G4double fact = 1.1;
  probmax *= fact;
  if(0.0 == fProb[nfilled] && nfilled > 2) { --nfilled; }
  for(size_t i=0; i<=nfilled; ++i) {
    fProb[i] *= fact;
  }
  G4double ekin(0.0), s1(1.0), s2(0.0), ksi(0.0), psi(1.0), z(0.0);

  if(fVerbose > 1) {
    G4cout << "### G4VEmissionProbability::SampleEnergy: " 
	   << " Emin= " << emin << " Emax= " << emax 
	   << " Nf= " << nfilled << G4endl;
  }
  G4double x0 = emax;
  G4double x1 = fEner[nfilled-1];
  G4double x2 = fEner[nfilled];
  G4double y0 = probmax;
  G4double y1 = fProb[nfilled-1];
  G4double y2 = fProb[nfilled];
  G4bool islog(false);
  // a condition if a special treatment of falling down 
  // distribution is needed 
  // 
  G4double ee = 0.5*(emax - emin);
  if(nfilled > 5 && y2 > 0.0 && y2 < 0.1*y0 && y1 > 0.0 && x1 - emin < ee) { 
    ksi = G4Log(y1/y2)/G4Log(x2/x1);
    x0 = x2*G4Exp(-G4Log(y0/y2)/ksi);
    // general condition to have two sampling area
    if(x0 < x1 && x0 > eprobmax && x0 - emin < ee) {
      // condition when the first area does not exist
      if(x0 <= emin) {
	s1  = 0.0;
	s2  = 1.0;
	y0 *= G4Exp(G4Log(x0/emin)*ksi);
	x0  = emin;
      } else { 
	s1  = (x0 - emin)*y0;
      }
      // parameters of the second area
      if(std::abs(1.0 - ksi) < 0.1) {
        z  = G4Log(emax/x0);
        if(s1 > 0.0) { s2 = x0*y0*z; }
        islog = true;
      } else {
	psi = 1.0/(1.0 - ksi);
	z   = G4Exp(G4Log(emax/x0)*(1. - ksi)) - 1.;
	if(s1 > 0.0) { s2  = std::max(y0*x0*z*psi, 0.0); }
      }
    }
  }
  G4double sum = s1 + s2;
  if(fVerbose > 1) {
    G4cout << " Epeak= " << eprobmax << " e0= " << x0 << " e1= " << x1 << " e2= " << x2  
	   <<  " psi= " << psi << G4endl; 
    G4cout << "   s1= " << s1 << " s2= " << s2 
	   << " y0= " << y0 << " y1= " << y1 << " y2= " << y2 << " z= " << z 
	   << G4endl; 
    for(size_t i=0; i<=nfilled; ++i) {
      G4cout << i << ".   E= " << fEner[i] << "  P= " << fProb[i] << G4endl;
    }
  }

  CLHEP::HepRandomEngine* rndm = G4Random::getTheEngine();
  static const G4int nmax = 100;
  G4double gr, g;
  for(size_t i=0; i<nmax; ++i) {
    G4double q = sum*rndm->flat();
    if(s2 <= 0.0 || q <= s1) {
      gr = y0;
      ekin = emin + (x0 - emin)*q/s1;
    } else if(islog) {
      ekin = x0*G4Exp((q - s1)*z/s2);
      gr   = y0*x0/ekin; 
    } else {
      ekin = x0*G4Exp(G4Log(1.0 + (q - s1)*z/s2)*psi);
      gr   = y0*G4Exp(G4Log(x0/ekin)*ksi);
    }
    g = ComputeProbability(ekin, eCoulomb);
    if(fVerbose > 1) {
      G4cout << "    " << i
	     << ". prob= " << g << " probmax= " << gr
	     << " Ekin= " << ekin << G4endl;
    }
    if(g > gr && fVerbose > 0) {
      G4cout << "### G4VEmissionProbability::SampleEnergy Warning i= " << i
	     << " prob/probmax= " << g/gr << "  rndm= " << q
	     << "\n    prob= " << g << " probmax= " << gr
	     << " z= " << z << " ksi= " << ksi
	     << "\n    Ekin= " << ekin << " Emin= " << emin
	     << " Emax= " << emax << " E0= " << x0 << G4endl;
      G4cout << "    s1= " << s1 << " s2= " << s2 
	     << " y0= " << y0 << " y1= " << y1 << " y2= " << y2 << G4endl; 
      for(size_t j=0; j<=nfilled; ++j) {
	G4cout << j << ".   E= " << fEner[j] << "  P= " << fProb[j] << G4endl;
      }
    }
    if(gr*rndm->flat() <= g) { break; }
  }
  return ekin;
}

