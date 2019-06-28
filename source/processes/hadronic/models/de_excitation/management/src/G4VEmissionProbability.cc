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
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// Modifications:
// 28.10.2010 V.Ivanchenko defined members in constructor and cleaned up

#include "G4VEmissionProbability.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"
#include "Randomize.hh"
#include "G4Exp.hh"
#include "G4Log.hh"

G4VEmissionProbability::G4VEmissionProbability(G4int Z, G4int A)
  :OPTxs(3),fVerbose(1),theZ(Z),theA(A),elimit(CLHEP::MeV),accuracy(0.02) 
{
  fG4pow = G4Pow::GetInstance();
  fPairCorr = G4NuclearLevelData::GetInstance()->GetPairingCorrection();
  length = nbin = 0;
  emin = emax = eCoulomb = totProbability = probmax = 0.0;
}

G4VEmissionProbability::~G4VEmissionProbability() 
{}

void G4VEmissionProbability::Initialise()
{
  G4DeexPrecoParameters* param = G4NuclearLevelData::GetInstance()->GetParameters();
  OPTxs = param->GetDeexModelType();
}

void G4VEmissionProbability::ResetIntegrator(size_t nbins, G4double de, G4double eps)
{
  if(nbins > 0) { length = nbins; }
  if(de > 0.0)  { elimit = de; }
  if(eps > 0.0) { accuracy = eps; }
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

  G4double edelta = elimit;
  nbin = (size_t)((emax - emin)/edelta) + 1;
  static const G4double edeltamin = 0.2*CLHEP::MeV;
  static const G4double edeltamax = 2*CLHEP::MeV;
  if(nbin < 4) { 
    nbin = 4;
    edelta = (emax - emin)/(G4double)nbin;
  } else if(nbin > length) {
    nbin = length;
  }

  G4double x(emin), del, y; 
  G4double edelmicro= edelta*0.02;
  probmax = ComputeProbability(x + edelmicro, eCoulomb);
  G4double problast = probmax;
  if(fVerbose > 2) {
    G4cout << "### G4VEmissionProbability::IntegrateProbability: " 
	   << " Emax= " << emax << " QB= " << cb << " nbin= " << nbin 
	   << G4endl;
    G4cout << "    0.  E= " << emin << "  prob= " << probmax << G4endl;
  }
  for(size_t i=1; i<=nbin; ++i) {
    x += edelta;
    if(x > emax) { 
      edelta += (emax - x);
      x = emax; 
    }
    G4bool endpoint = (std::abs(x - emax) < edelmicro) ? true : false;
    G4double xx = endpoint ? x - edelmicro : x;
    y = ComputeProbability(xx, eCoulomb);
    if(fVerbose > 2) { 
      G4cout << "    " << i << ".  E= " << x << "  prob= " << y 
	     << " Edel= " << edelta << G4endl;
    } 
    probmax = std::max(probmax, y);
    del = (y + problast)*edelta*0.5;
    totProbability += del;
    // end of the loop
    if(del < accuracy*totProbability || endpoint) { break; }
    problast = y;

    // smart step definition
    if(del != totProbability && del > 0.8*totProbability && 0.7*edelta > edeltamin) { 
      edelta *= 0.7;
    } else if(del < 0.1*totProbability && 1.5*edelta < edeltamax) { 
      edelta *= 1.5;
    }
  }

  if(fVerbose > 1) { G4cout << " Probability= " << totProbability << G4endl; }
  return totProbability;
}

G4double G4VEmissionProbability::SampleEnergy()
{
  static const G4double fact = 1.05;
  probmax *= fact;

  if(fVerbose > 1) {
    G4cout << "### G4VEmissionProbability::SampleEnergy: " 
	   << " Emin= " << emin << " Emax= " << emax 
	   << " probmax= " << probmax << G4endl;
  }

  CLHEP::HepRandomEngine* rndm = G4Random::getTheEngine();
  static const G4int nmax = 100;
  G4double del = emax - emin;
  G4double ekin, g;
  G4int n = 0;
  do {
    ekin = del*rndm->flat() + emin; 
    ++n;
    g = ComputeProbability(ekin, eCoulomb);
    if(fVerbose > 2) {
      G4cout << "    " << n
	     << ". prob= " << g << " probmax= " << probmax
	     << " Ekin= " << ekin << G4endl;
    }
    if((g > probmax || n > nmax) && fVerbose > 1) {
      G4cout << "### G4VEmissionProbability::SampleEnergy for Z= " << theZ 
             << " A= " << theA 
             << "\n    Warning n= " << n
	     << " prob/probmax= " << g/probmax 
	     << " prob= " << g << " probmax= " << probmax 
	     << "\n    Ekin= " << ekin << " Emin= " << emin
	     << " Emax= " << emax << G4endl;
    }
  } while(probmax*rndm->flat() > g && n < nmax);
  return ekin;
}

