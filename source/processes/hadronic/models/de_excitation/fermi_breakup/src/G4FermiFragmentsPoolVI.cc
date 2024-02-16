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
// FermiBreakUp de-excitation model
// by V. Ivanchenko (July 2016)
//

#include "G4FermiFragmentsPoolVI.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4NuclearLevelData.hh"
#include "G4NucleiProperties.hh"
#include "G4LevelManager.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4FermiBreakUpUtil.hh"
#include <iomanip>

G4FermiFragmentsPoolVI::G4FermiFragmentsPoolVI()
{}

G4FermiFragmentsPoolVI::~G4FermiFragmentsPoolVI()
{
  for (G4int i=0; i<maxA; ++i) {
    for (G4int j=0; j<maxZ; ++j) {
      auto ptr = list_c[j][i];
      if (nullptr != ptr) {
	for ( auto const & p : *ptr) { delete p; }
        delete ptr;
      }
    }
  }
  for (auto const & ptr : fragment_pool) { delete ptr; }
}

G4bool G4FermiFragmentsPoolVI::HasDecay(const G4int Z, const G4int A,
                                        const G4double eexc) const
{
  if (Z < maxZ && A < maxA && nullptr != list_c[Z][A]) {
    for (auto const & ch : *(list_c[Z][A])) {
      if (ch->GetExcitation() <= eexc + fTolerance && ch->NumberPairs() > 0) {
	return true;
      } 
    }
  }
  return false;
}

const G4FermiChannels* 
G4FermiFragmentsPoolVI::ClosestChannels(G4int Z, G4int A, G4double etot) const
{
  const G4FermiChannels* res = nullptr;
  if (Z >= maxZ || A >= maxA) { return res; } 

  auto chan = list_c[Z][A];
  if (nullptr == chan) { return res; }

  G4double demax = 1.e+9;
  for (auto const & ch : *chan) {
    if (ch->NumberPairs() == 0) { continue; }
    G4double de = etot - ch->GetFragment()->GetTotalEnergy();
    // an excitation coincide with a level
    if (std::abs(de) <= fTolerance) {
      return ch;
    }
    if (de >= 0 && de < demax) {
      demax = de;
      res = ch;
    } 
  }
  return res;
}

G4bool G4FermiFragmentsPoolVI::IsInThePool(const G4int Z, const G4int A,  
                                           const G4double exc) const
{
  for (auto const& fr : fragment_pool) {
    if(fr->GetZ() == Z && fr->GetA() == A &&  
       std::abs(exc - fr->GetExcitationEnergy()) < fTolerance) 
      { return true; }
  }
  return false;
}

void G4FermiFragmentsPoolVI::Initialise()
{
  if (isInitialized) { return; }
  isInitialized = true;
  G4DeexPrecoParameters* param = 
    G4NuclearLevelData::GetInstance()->GetParameters();
  fTolerance = 2*CLHEP::eV;
  fElim = param->GetFBUEnergyLimit();

  fragment_pool.reserve(991);

  // stable particles
  fragment_pool.push_back(new G4FermiFragment(1, 0, 1, 0.0, DBL_MAX));
  fragment_pool.push_back(new G4FermiFragment(1, 1, 1, 0.0, DBL_MAX));
  fragment_pool.push_back(new G4FermiFragment(2, 1, 2, 0.0, DBL_MAX));
  fragment_pool.push_back(new G4FermiFragment(3, 1, 1, 0.0, 3.8879e+08));
  fragment_pool.push_back(new G4FermiFragment(3, 2, 1, 0.0, DBL_MAX));
  fragment_pool.push_back(new G4FermiFragment(4, 2, 0, 0.0, DBL_MAX));
  fragment_pool.push_back(new G4FermiFragment(5, 2, 3, 0.0, 7.0325e-22));
  fragment_pool.push_back(new G4FermiFragment(5, 3, 3, 0.0, 3.70493e-22));

  // use level data and construct the pool
  G4NuclearLevelData* ndata = G4NuclearLevelData::GetInstance();
  for (G4int Z=1; Z<maxZ; ++Z) {
    G4int Amin = ndata->GetMinA(Z);
    G4int Amax = std::min(maxA, ndata->GetMaxA(Z)+1);
    for (G4int A=Amin; A<Amax; ++A) {
      const G4LevelManager* man = ndata->GetLevelManager(Z, A);
      if (nullptr != man) {
        std::size_t nn = man->NumberOfTransitions();
        for(std::size_t i=0; i<=nn; ++i) {
          G4double exc = man->LevelEnergy(i);
	  /*          
          G4cout << "++ Z=" << Z << " A=" << A << " Eex=" << exc 
                 << " time(ns)=" << man->LifeTime(i)/ns << " i=" << i
                 << " InPool=" << IsInThePool(Z, A, exc) << G4endl;
	  */
          // only levels below limit are consided 
          if (exc >= fElim) { continue; }
          // only new are considered
          if (IsInThePool(Z, A, exc)) { continue; }
          fragment_pool.push_back(new G4FermiFragment(A, Z, man->TwoSpinParity(i),
                                                      exc, man->LifeTime(i)));
        }
      }
    }
  }
  G4int nfrag = (G4int)fragment_pool.size();
  for (auto const& f : fragment_pool) {
    G4int Z = f->GetZ();
    G4int A = f->GetA();
    if (list_c[Z][A] == nullptr) {
      list_c[Z][A] = new std::vector<G4FermiChannels*>;
    }
    (list_c[Z][A])->push_back(new G4FermiChannels(f));
  }

  // list of fragment pairs ordered by A
  for (G4int i=0; i<nfrag; ++i) {
    const G4FermiFragment* f1 = fragment_pool[i];
    G4int Z1 = f1->GetZ();
    G4int A1 = f1->GetA();
    G4double e1 = f1->GetTotalEnergy();
    for (G4int j=0; j<nfrag; ++j) {
      const G4FermiFragment* f2 = fragment_pool[j];
      G4int Z2 = f2->GetZ();
      G4int A2 = f2->GetA();
      if(A2 < A1 || (A2 == A1 && Z2 < Z1)) { continue; }
      G4int Z = Z1 + Z2;
      G4int A = A1 + A2;

      if(Z >= maxZ || A >= maxA) { continue; } 

      G4double e2 = f2->GetTotalEnergy();
      G4double minE = e1 + e2;
      G4double exc = minE - G4NucleiProperties::GetNuclearMass(A, Z);
      /*
      if(1 == Z) {
        G4cout << "+!+ Z=" << Z << " A=" << A 
             << " Z1=" << Z1 << " A1=" << A1
             << " Z2=" << Z2 << " A2=" << A2 << " Eex=" << exc 
             << "  Qb=" << G4FermiBreakUpUtil::CoulombBarrier(Z1, A1, Z2, A2, exc)
             << " e1=" << e1 << " e2=" << e2
             << " M=" << G4NucleiProperties::GetNuclearMass(A, Z)
             << G4endl; 
      }
      */
      // ignore very excited case
      if (exc > fElim) { continue; }
      auto chan = list_c[Z][A];
      if (nullptr == chan) { continue; }
      std::size_t kmax = chan->size();
      for (std::size_t k=0; k<kmax; ++k) {
	auto ch = (*chan)[k];
        const G4double e0 = ch->GetMass();
        auto f0 = ch->GetFragment();
        if (e0 > minE && G4FermiBreakUpUtil::CheckSpinParity(f1, f2, f0)) { 
          const G4double cb =
            G4FermiBreakUpUtil::CoulombBarrier(Z1, A1, Z2, A2, ch->GetExcitation());
          if (e0 >= minE + cb) {
	    ch->AddChannel(new G4FermiPair(f1, f2));
	  }
        }
      }
    }
  }
  // compute cumulative probabilities
  for (G4int A=1; A<maxA; ++A) {
    for (G4int Z=0; Z<maxZ; ++Z) {
      auto chan = list_c[Z][A];
      if(nullptr == chan) { continue; }
      std::size_t kmax = chan->size();
      for (std::size_t k=0; k<kmax; ++k) {
	auto ch = (*chan)[k];
	auto frag = ch->GetFragment();
	std::size_t nch = ch->NumberPairs();
	if (1 < nch) {
	  const std::vector<G4FermiPair*>& pairs = ch->GetChannels();
	  G4double ptot = 0.0;
	  for (std::size_t i=0; i<nch; ++i) {
	    ptot += G4FermiBreakUpUtil::Probability(frag->GetA(),
                                                    pairs[i]->GetFragment1(),
                                                    pairs[i]->GetFragment2(),
                                                    frag->GetTotalEnergy(),
						    frag->GetExcitationEnergy());
            pairs[i]->SetProbability(ptot);
	  }
	  // normalisation
	  if (0.0 == ptot) {
	    pairs[0]->SetProbability(1.0);
	  } else {
	    ptot = 1./ptot;
	    for (std::size_t i=0; i<nch-1; ++i) {
	      G4double x = ptot*pairs[i]->Probability();
	      pairs[i]->SetProbability(x);
	    }
	    pairs[nch - 1]->SetProbability(1.0);
          }
        }
      }
    }
  }
}

void G4FermiFragmentsPoolVI::DumpFragment(const G4FermiFragment* f) const
{
  if (nullptr != f) {
    G4long prec = G4cout.precision(6);
    G4cout << "   Z=" << f->GetZ() << " A=" << std::setw(2) << f->GetA() 
           << " Mass(GeV)=" << std::setw(8) << f->GetFragmentMass()/GeV
           << " Eexc(MeV)=" << std::setw(7) << f->GetExcitationEnergy()
           << " 2S=" << f->TwoSpinParity() << G4endl;
    G4cout.precision(prec);
  }
}

void G4FermiFragmentsPoolVI::Dump() const
{
  G4cout <<"----------------------------------------------------------------"
         <<G4endl;
  G4cout << "##### List of Fragments in the Fermi Fragment Pool #####" 
         << G4endl;
  std::size_t nfrag = fragment_pool.size();
  G4cout << "      Nfragnents=" << nfrag << " Elim(MeV)=" << fElim/CLHEP::MeV << G4endl; 
  for(std::size_t i=0; i<nfrag; ++i) {
    DumpFragment(fragment_pool[i]);
  }
  G4cout << G4endl;


  G4cout << "----------------------------------------------------------------"
         << G4endl;
  G4cout << "### G4FermiFragmentPoolVI: fragments sorted by A" << G4endl; 

  G4int ama{0};
  G4long prec = G4cout.precision(6);
  for (G4int A=1; A<maxA; ++A) {
    for (G4int Z=0; Z<maxZ; ++Z) {
      auto chan = list_c[Z][A];
      if (nullptr == chan) { continue; }
      std::size_t jmax = chan->size();
      G4cout << " # A=" << A << "  Z=" << Z << "  Nfagments=" << jmax << G4endl;
      for(std::size_t j=0; j<jmax; ++j) {
	auto ch = (*chan)[j];
	if(nullptr == ch) { continue; }
	auto f = ch->GetFragment();
	G4int a1 = f->GetA();
	G4int z1 = f->GetZ();
	std::size_t nch = ch->NumberPairs();
	ama += nch;
	G4cout << "   ("<<a1<<","<<z1<<");  Eex(MeV)= "
	       << f->GetExcitationEnergy() 
	       << " 2S=" << f->TwoSpinParity()
	       << "; Nchannels=" << nch
	       << G4endl; 
	for (std::size_t k=0; k<nch; ++k) {
          auto fpair = ch->GetPair(k);
	  if(nullptr == fpair) { continue; }
	  G4cout << "         (" << fpair->GetFragment1()->GetZ()
		 << ", "  << fpair->GetFragment1()->GetA() 
		 << ", " << fpair->GetFragment1()->GetExcitationEnergy()
		 << ")  ("<< fpair->GetFragment2()->GetZ()
		 << ", "  << fpair->GetFragment2()->GetA() << ", "
		 << fpair->GetFragment2()->GetExcitationEnergy()
		 << ")  prob= " << fpair->Probability()
		 << G4endl;
	}
      }
    }
  }
  G4cout.precision(prec);
  G4cout << " ======== Total number of channels " << ama << "  ======" << G4endl;
}

