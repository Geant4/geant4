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
// $Id: G4VFermiBreakUp.cc,v 1.5 2006-06-29 20:13:13 gunter Exp $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)
//
// Modifications:
// 01.04.2011 General cleanup by V.Ivanchenko - use only one object
//         theConfigurationList and do not instantiate it at each decay

#include "G4FermiBreakUp.hh"
#include "G4FermiFragmentsPool.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4Pow.hh"

G4FermiBreakUp::G4FermiBreakUp() : thePhaseSpace(nullptr)
{
  g4calc = G4Pow::GetInstance();
  Coef = ConstCoeff = 0.0;

  nmax = 16;
  NormalizedWeights.resize(nmax,0.0);
  massRes.reserve(4);
  frag.resize(4, 0);
  thePool = G4FermiFragmentsPool::Instance(); 
  Initialise();
}

G4FermiBreakUp::~G4FermiBreakUp()
{}

void G4FermiBreakUp::Initialise()
{
  if(thePhaseSpace != nullptr) { return; }
  // Kappa = V/V_0 it is used in calculation of Coulomb energy
  // Nuclear radius r0 (is a model parameter)
  G4double Kappa = 6.0;
  G4double r0 = 1.3*CLHEP::fermi;

  Coef = 0.6*(CLHEP::elm_coupling/r0)/g4calc->Z13(1+G4int(Kappa));
  ConstCoeff = g4calc->powN(r0/hbarc,3)*Kappa/(6.0*pi*pi);

  thePhaseSpace = thePool->GetFermiPhaseSpaceDecay();
}

G4bool G4FermiBreakUp::IsApplicable(G4int Z, G4int A, G4double mass) const
{
  return thePool->IsApplicable(Z, A, mass);
}

G4double G4FermiBreakUp::CoulombBarrier(
  const std::vector<const G4VFermiFragment*>* conf)
{
  //  Calculates Coulomb Barrier (MeV) for given channel with K fragments.
  G4int SumA = 0;
  G4int SumZ = 0;
  G4double CoulombEnergy = 0.;
  size_t nn = conf->size();
  for (size_t i=0; i<nn; ++i) {
    G4int z = (*conf)[i]->GetZ();
    G4int a = (*conf)[i]->GetA();
    CoulombEnergy += G4double(z*z)/g4calc->Z13(a);
    SumA += a;
    SumZ += z;
  }
  CoulombEnergy -= SumZ*SumZ/g4calc->Z13(SumA);
  return -Coef * CoulombEnergy;
}

void G4FermiBreakUp::BreakFragment(G4FragmentVector* theResult, 
				   G4Fragment* theNucleus)
{
  if(thePool == nullptr) { Initialise(); }
  // Calculate Momenta of K fragments
  G4double M = theNucleus->GetMomentum().m();
  const std::vector<const G4VFermiFragment*>* conf = 
    SelectConfiguration(theNucleus->GetZ_asInt(), 
			theNucleus->GetA_asInt(), M);

  // should never happen
  if(!conf) { 
    theResult->push_back(theNucleus);
    return; 
  }

  size_t nn = conf->size();

  // should never happen
  if(0 == nn) {
    theResult->push_back(theNucleus);
    return; 
  }

  G4LorentzVector fourMomentum = theNucleus->GetMomentum();

  // one unstable fragment
  if(1 == nn) {
    (*conf)[0]->FillFragment(theResult, fourMomentum);

    // normal case
  } else {
    G4ThreeVector boostVector = fourMomentum.boostVector();  
    massRes.clear();
    for(size_t i=0; i<nn; ++i) {
      massRes.push_back( (*conf)[i]->GetTotalEnergy() );
    }
    std::vector<G4LorentzVector*>* mom = thePhaseSpace->Decay(M, massRes);

    //  size_t nmom = mom->size();  
    // G4cout << "nmom= " << nmom << G4endl;

    // Go back to the Lab Frame
    for (size_t j=0; j<nn; ++j) {    
      (*mom)[j]->boost(boostVector); 
      (*conf)[j]->FillFragment(theResult, *((*mom)[j]));
      delete (*mom)[j];
    }
    delete mom;
  }
  delete theNucleus;
}

const std::vector<const G4VFermiFragment*>* 
G4FermiBreakUp::SelectConfiguration(G4int Z, G4int A, G4double mass)
{
  const std::vector<const G4VFermiFragment*>* res = 0;
  //   new std::vector<const G4VFermiFragment*>;  
  const std::vector<const G4FermiConfiguration*>* conflist = 
    thePool->GetConfigurationList(Z, A, mass);
  if(!conflist) { return res; }
  size_t nn = conflist->size();
  if(0 < nn) { 
    size_t idx = 0;
    if(1 < nn) {
      if(nn > nmax) {
        nmax = nn;
	NormalizedWeights.resize(nmax,0.0);
      }
      G4double prob = 0.0;
      for(size_t i=0; i<nn; ++i) { 
	prob += DecayProbability(A, mass, (*conflist)[i]);
	NormalizedWeights[i] = prob;
      }
      prob *= G4UniformRand();
      for(idx=0; idx<nn; ++idx) { 
	if(NormalizedWeights[idx] >= prob) { break; }
      }
    }
    //const std::vector<const G4VFermiFragment*> flist = 
    res = (*conflist)[idx]->GetFragmentList();
    //size_t nf = flist.size();
    //for(size_t i=0; i<nf; ++i) { res->push_back(flist[i]); }
    //G4cout << "FermiBreakUp: " << nn << " conf; selected idx= " 
    //	 << idx << "  Nprod= " << nf << G4endl; 
  }
  delete conflist;
  return res;
}

G4double G4FermiBreakUp::DecayProbability(G4int A, G4double TotalE, 
					  const G4FermiConfiguration* conf)
  // Decay probability  for a given channel with K fragments
{
  // A: Atomic Weight
  // TotalE: Total energy of nucleus

  G4double KineticEnergy = TotalE; 
  const std::vector<const G4VFermiFragment*>* flist = 
    conf->GetFragmentList();

  // Number of fragments
  size_t K = flist->size();
  if(K > frag.size()) { frag.resize(K, 0); }

  for (size_t i=0; i<K; ++i) {
    frag[i] = (*flist)[i];
    KineticEnergy -= 
      ((*flist)[i]->GetFragmentMass() + (*flist)[i]->GetExcitationEnergy());
  }

  // Check that there is enough energy to produce K fragments
  KineticEnergy -= CoulombBarrier(flist);
  if (KineticEnergy <= 0.0) { return 0.0; }

  // Spin factor S_n
  G4double S_n = 1.0;

  // mass factors
  G4double ProdMass = 1.0;
  G4double SumMass = 0.0;

  for (size_t i=0; i<K; ++i) {
    G4double mass = (*flist)[i]->GetFragmentMass();
    ProdMass *= mass;
    SumMass += mass;
    S_n *= (*flist)[i]->GetPolarization();
  }

  G4double MassFactor = ProdMass/SumMass;
  MassFactor *= std::sqrt(MassFactor);  
  
  G4double Coeff = g4calc->powN(ConstCoeff*A, K-1);

  //JMQ 111009 Bug fixed: gamma function for odd K was wrong  by a factor 2
  //VI  251014 Use G4Pow
  G4double Gamma = 1.0;
  G4double Energ = twopi*KineticEnergy;

  // integer argument of Gamma function
  if ((K-1)%2 != 1) {
    G4int N = 3*(K-1)/2;
    Gamma = g4calc->factorial(N - 1);
    Energ = g4calc->powN(Energ, N);
    
    // n + 1/2 argument of Gamma function
    // http://en.wikipedia.org/wiki/Gamma_function
  } else   { 
    G4int n2 = 3*K - 4;
    G4int n1 = n2/2;

    static const G4double sqrtpi = std::sqrt(CLHEP::pi);
    Gamma = sqrtpi*g4calc->factorial(n2)/
      (g4calc->powN(4.0,n1)*g4calc->factorial(n1));
    Energ = g4calc->powN(Energ, n1)*std::sqrt(Energ);
  }
  
  // Permutation Factor G_n
  // search for identical fragments
  G4double G_n = 1.0;
  for(size_t i=0; i<K-1; ++i) {
    if(frag[i]) {
      G4int nf = 1;
      for(size_t j=i+1; j<K; ++j) {
        if(frag[i] == frag[j]) {
	  ++nf;
          frag[j] = 0;
	}
      }
      if(1 < nf) { G_n *= g4calc->factorial(nf); }
    }
  }

  G4double Weight = Coeff*MassFactor*S_n*Energ/(G_n*Gamma*KineticEnergy);

  return Weight; 
}

