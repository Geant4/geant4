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
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)
//
// Modifications:
// 01.04.2011 General cleanup by V.Ivanchenko - more clean usage of static
// 23.04.2011 V.Ivanchenko: make this class to be responsible for
//            selection of decay channel and decay

#include <set>

#include "G4FermiConfigurationList.hh"
#include "G4FermiFragmentsPool.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4Pow.hh"

const G4double G4FermiConfigurationList::Kappa = 6.0;
const G4double G4FermiConfigurationList::r0 = 1.3*CLHEP::fermi;

G4FermiConfigurationList::G4FermiConfigurationList()
{
  thePool = G4FermiFragmentsPool::Instance(); 
  g4pow = G4Pow::GetInstance();
  Coef = 0.6*(CLHEP::elm_coupling/r0)/g4pow->Z13(1+G4int(Kappa));
  ConstCoeff = g4pow->powN(r0/hbarc,3)*Kappa*std::sqrt(2.0/pi)/3.0;

  // 16 is the max number
  nmax = 50;
  NormalizedWeights.resize(nmax,0.0);
}

G4FermiConfigurationList::~G4FermiConfigurationList()
{}

G4double 
G4FermiConfigurationList::CoulombBarrier(
  const std::vector<const G4VFermiFragment*>& conf)
{
  //  Calculates Coulomb Barrier (MeV) for given channel with K fragments.
  G4int SumA = 0;
  G4int SumZ = 0;
  G4double CoulombEnergy = 0.;
  size_t nn = conf.size();
  for (size_t i=0; i<nn; ++i) {
    G4int z = conf[i]->GetZ();
    G4int a = conf[i]->GetA();
    CoulombEnergy += G4double(z*z)/g4pow->Z13(a);
    SumA += a;
    SumZ += z;
  }
  CoulombEnergy -= SumZ*SumZ/g4pow->Z13(SumA);
  return -Coef * CoulombEnergy;
}

G4double 
G4FermiConfigurationList::DecayProbability(G4int A, G4double TotalE, 
					   G4FermiConfiguration* conf)
  // Decay probability  for a given channel with K fragments
{
  // A: Atomic Weight
  // TotalE: Total energy of nucleus

  G4double KineticEnergy = TotalE; // MeV
  G4double ProdMass = 1.0;
  G4double SumMass = 0.0;
  G4double S_n = 1.0;
  std::set<G4int>      combSet;
  std::multiset<G4int> combmSet;

  const std::vector<const G4VFermiFragment*> flist = 
    conf->GetFragmentList();

  // Number of fragments
  G4int K = flist.size();

  for (G4int i=0; i<K; ++i) {
    G4int a = flist[i]->GetA();
    combSet.insert(a);
    combmSet.insert(a);
    G4double mass = flist[i]->GetFragmentMass();
    ProdMass *= mass;
    SumMass += mass;
    // Spin factor S_n
    S_n *= flist[i]->GetPolarization();
    KineticEnergy -= mass + flist[i]->GetExcitationEnergy();
  }

  // Check that there is enough energy to produce K fragments
  KineticEnergy -= CoulombBarrier(flist);
  if (KineticEnergy <= 0.0) { return 0.0; }

  G4double MassFactor = ProdMass/SumMass;
  MassFactor *= std::sqrt(MassFactor);  
  
  // This is the constant (doesn't depend on nucleus) part
  G4double Coeff = g4pow->powN(ConstCoeff*A, K-1);

  //JMQ 111009 Bug fixed: gamma function for odd K was wrong  by a factor 2
  // new implementation explicitely according to standard properties of Gamma function
  // Calculation of 1/Gamma(3(k-1)/2)
  G4double Gamma = 1.0;
  //  G4double arg = 3.0*(K-1)/2.0 - 1.0;
  //  while (arg > 1.1) 
  //    {
  //      Gamma *= arg; 
  //      arg--;
  //    }
  //  if ((K-1)%2 == 1) Gamma *= std::sqrt(pi);

  if ((K-1)%2 != 1)

    {
      G4double arg = 3.0*(K-1)/2.0 - 1.0;
      while (arg > 1.1) 
	{
	  Gamma *= arg; 
	  arg--;
	}
    }
  else   { 
    G4double 	n= 3.0*K/2.0-2.0;
    G4double arg2=2*n-1;
    while (arg2>1.1)
      {
	Gamma *= arg2;
	arg2  -= 2;
      }
    Gamma *= std::sqrt(pi)/g4pow->powZ(2,n);
  }
  // end of new implementation of Gamma function
  
  
  // Permutation Factor G_n
  G4double G_n = 1.0;
  for (std::set<G4int>::iterator itr = combSet.begin(); 
       itr != combSet.end(); ++itr)
    {
      for (G4int ni = combmSet.count(*itr); ni > 1; ni--) { G_n *= ni; }
    }

  G4double Weight = Coeff * MassFactor * (S_n / G_n) / Gamma;
  Weight *=  std::sqrt(g4pow->powN(KineticEnergy,3*(K-1)))/KineticEnergy;

  return Weight; 
}

G4FragmentVector* 
G4FermiConfigurationList::GetFragments(const G4Fragment& theNucleus)
{
  // Calculate Momenta of K fragments
  G4double M = theNucleus.GetMomentum().m();
  const std::vector<const G4VFermiFragment*>* conf = 
    SelectConfiguration(theNucleus.GetZ_asInt(), 
			theNucleus.GetA_asInt(), M);


  G4FragmentVector* theResult = new G4FragmentVector();
  size_t nn = conf->size();
  if(1 >= nn) {
    theResult->push_back(new G4Fragment(theNucleus));
    delete conf;
    return theResult; 
  }

  G4ThreeVector boostVector = theNucleus.GetMomentum().boostVector();  
  std::vector<G4double> mr;
  mr.reserve(nn);
  for(size_t i=0; i<nn; ++i) {
    mr.push_back( (*conf)[i]->GetTotalEnergy() );
  }
  std::vector<G4LorentzVector*>* mom = thePhaseSpace.Decay(M,mr);
  if(!mom) { 
    delete conf;
    return theResult; 
  }

  size_t nmom = mom->size();  

  // Go back to the Lab Frame
  if(0 < nmom) { 
    for (size_t j=0; j<nmom; ++j) {
      G4LorentzVector* FourMomentum = (*mom)[j];
    
      // Lorentz boost
      FourMomentum->boost(boostVector);
      
      G4FragmentVector* fragment = (*conf)[j]->GetFragment(*FourMomentum);
 
      size_t nfrag = fragment->size();
      for (size_t k=0; k<nfrag; ++k) { theResult->push_back((*fragment)[k]); }
      delete fragment;
      delete (*mom)[j];
    }
  }
  
  delete mom;
  delete conf;
  return theResult;
}

const std::vector<const G4VFermiFragment*>* 
G4FermiConfigurationList::SelectConfiguration(G4int Z, G4int A, G4double mass)
{
  std::vector<const G4VFermiFragment*>* res = 
    new std::vector<const G4VFermiFragment*>;  
  const std::vector<G4FermiConfiguration*>* conflist = 
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
    const std::vector<const G4VFermiFragment*> flist = 
      (*conflist)[idx]->GetFragmentList();
    size_t nf = flist.size();
    for(size_t i=0; i<nf; ++i) { res->push_back(flist[i]); }
    //G4cout << "FermiBreakUp: " << nn << " conf; selected idx= " 
    //	 << idx << "  Nprod= " << nf << G4endl; 
  }
  delete conflist;
  return res;
}
