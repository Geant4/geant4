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

#include "G4FermiBreakUpVI.hh"
#include "G4FermiBreakUpUtil.hh"
#include "G4FermiFragmentsPoolVI.hh"
#include "G4FermiChannels.hh"
#include "G4FermiPair.hh"
#include "G4PhysicalConstants.hh"
#include "G4NuclearLevelData.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4PhysicsModelCatalog.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4AutoLock.hh"

G4FermiFragmentsPoolVI* G4FermiBreakUpVI::fPool = nullptr;

namespace
{
  G4Mutex theFBUMutex = G4MUTEX_INITIALIZER;
}

G4FermiBreakUpVI::G4FermiBreakUpVI()
{
  frag.reserve(10);
  lvect.reserve(10);
  secID = G4PhysicsModelCatalog::GetModelID("model_G4FermiBreakUpVI");
  prob.resize(12,0.0);
  if (nullptr == fPool) {
    G4AutoLock l(&theFBUMutex);
    if (nullptr == fPool) {
      fPool = new G4FermiFragmentsPoolVI();
      fPool->Initialise();
      isFirst = true;
    }
    l.unlock();
  }
}

G4FermiBreakUpVI::~G4FermiBreakUpVI()
{
  if (isFirst) { 
    delete fPool;
    fPool = nullptr;
  }
}

void G4FermiBreakUpVI::Initialise()
{
  G4DeexPrecoParameters* param = 
    G4NuclearLevelData::GetInstance()->GetParameters();
  fTolerance = param->GetMinExcitation();
  fElim = param->GetFBUEnergyLimit();
  if (verbose > 1) {
    G4cout << "### G4FermiBreakUpVI::Initialise(): the pool is initilized=" 
	   << fPool->IsInitialized() << " fTolerance(eV)=" << fTolerance/CLHEP::eV
           << " Elim(MeV)=" << fElim/CLHEP::MeV << G4endl;
  }
}

G4bool G4FermiBreakUpVI::IsApplicable(G4int Z, G4int A, G4double eexc) const
{
  return (Z < maxZ && A < maxA && eexc <= fElim && fPool->HasDecay(Z, A, eexc));
}

void G4FermiBreakUpVI::BreakFragment(G4FragmentVector* theResult, 
				     G4Fragment* theNucleus)
{
  if (verbose > 1) {
    G4cout << "### G4FermiBreakUpVI::BreakFragment start new fragment " 
           << G4endl;
    G4cout << *theNucleus << G4endl;
  }
  if (!fPool->IsInitialized()) { fPool->Initialise(); } 

  // initial fragment
  G4int Z = theNucleus->GetZ_asInt();
  G4int A = theNucleus->GetA_asInt();
  G4double excitation = theNucleus->GetExcitationEnergy();
  if (!IsApplicable(Z, A, excitation)) { return; }
  G4double mass = theNucleus->GetGroundStateMass() + excitation;
  G4LorentzVector lv0 = theNucleus->GetMomentum();

  // sample first decay of an initial state
  // if not possible to decay - exit
  if (!SampleDecay(Z, A, mass, excitation, lv0)) { return; }

  G4double time = theNucleus->GetCreationTime();
  delete theNucleus;

  static const G4int imax = 100; 

  // loop over vector of Fermi fragments
  // vector may grow at each iteraction
  for (std::size_t i=0; i<frag.size(); ++i) {
    Z = frag[i]->GetZ();
    A = frag[i]->GetA();
    excitation = frag[i]->GetExcitationEnergy();
    lv0 = lvect[i];
    G4bool unstable = IsApplicable(Z, A, excitation);
    if (unstable) { 
      mass = frag[i]->GetTotalEnergy();
      if (verbose > 1) {
	G4cout << "# FermiFrag " << i << ".  Z= " << Z << " A= " << A 
	       << " mass= " << mass << " exc= " 
	       << frag[i]->GetExcitationEnergy() << G4endl;
      }
      unstable = SampleDecay(Z, A, mass, excitation, lv0);
    }
    // stable fragment
    if (!unstable) {
      if(verbose > 1) { G4cout << "   New G4Fragment" << G4endl; }
      G4Fragment* f = new G4Fragment(A, Z, lv0);
      f->SetCreationTime(time);
      f->SetCreatorModelID(secID);
      theResult->push_back(f);
    }
    // limit the loop
    if (i == imax) { break; }
  }
  frag.clear();
  lvect.clear();
}

G4bool G4FermiBreakUpVI::SampleDecay(const G4int Z, const G4int A, const G4double mass,
                                     const G4double exc, G4LorentzVector& lv0)
{
  const G4FermiChannels* chan = fPool->ClosestChannels(Z, A, mass);
  if (nullptr == chan) { return false; }
  std::size_t nn = chan->NumberPairs();
  if (verbose > 1) {
    G4cout << "G4FermiBreakUpVI::SampleDecay " << nn << " channels Eex= " 
	   << chan->GetExcitation() << G4endl;
  }
  if (0 == nn) { return false; }
  if (nn > prob.size()) { prob.resize(nn, 0.0); }

  const G4FermiPair* fpair = nullptr;

  // one unstable fragment
  if (1 == nn) {
    fpair = chan->GetPair(0);

    // more pairs
  } else {
    
    G4double q = G4UniformRand();
    const std::vector<G4FermiPair*>& pvect = chan->GetChannels();
    std::size_t i{0};
    G4bool pre = true; 
    if (std::abs(exc - chan->GetExcitation()) < fTolerance) {
      // static probabilities may be used
      for (; i<nn; ++i) {
	if (q <= pvect[i]->Probability()) {
	  fpair = pvect[i];
	  break;
	}
      }
    } else {
      // recompute probabilities
      pre = false;
      G4double ptot = 0.0;
      for (i=0; i<nn; ++i) {
	ptot += G4FermiBreakUpUtil::Probability(A, pvect[i]->GetFragment1(),
					        pvect[i]->GetFragment2(),
                                                mass, exc);
	prob[i] = ptot;
      }
      ptot *= q;
      for (i=0; i<nn; ++i) {
        if(ptot <= prob[i]) {
          fpair = pvect[i];
	  break;
	}
      }
    }
    if (verbose > 2) {
      G4cout << "Probabilities of 2-body decay: Nchannels=" << nn 
             << " channels; i=" << i << " is selected; predefined=" 
	     << pre << G4endl;
      for (std::size_t j=0; j<nn; ++j) {
	G4cout << j << ". "; 
	if (pre) { G4cout << pvect[j]->Probability(); }
	else { G4cout << prob[j]; }
	G4cout << " Z1= " << pvect[j]->GetFragment1()->GetZ()
	       << " A1= " << pvect[j]->GetFragment1()->GetA()
	       << " Z2= " << pvect[j]->GetFragment2()->GetZ()
	       << " A2= " << pvect[j]->GetFragment2()->GetA()
	       << G4endl;
      }
    }
  }
  if (nullptr == fpair) { return false; }

  auto frag1 = fpair->GetFragment1();
  auto frag2 = fpair->GetFragment2();
  
  G4double mass1 = frag1->GetTotalEnergy();
  G4double mass2 = frag2->GetTotalEnergy();
  if (verbose > 2) {
    G4cout << " M= " << mass << " M1= " << mass1 << "  M2= " << mass2 
	   << " Exc1= " << frag1->GetExcitationEnergy() 
	   << " Exc2= " << frag2->GetExcitationEnergy() << G4endl;
  }
  // sample fragment1
  G4double e1 = 0.5*(mass*mass - mass2*mass2 + mass1*mass1)/mass;
  //G4cout << "ekin1= " << e1 - mass1 << G4endl;
  G4double p1(0.0);
  if (e1 > mass1) {
    p1 = std::sqrt((e1 - mass1)*(e1 + mass1));
  } else {
    e1 = mass1;
  }
  G4LorentzVector lv1 = G4LorentzVector(G4RandomDirection()*p1, e1);

  // compute kinematics
  auto boostVector = lv0.boostVector();  
  lv1.boost(boostVector);
  G4LorentzVector lv2 = lv0 - lv1;

  frag.push_back(frag1);
  frag.push_back(frag2);
  lvect.push_back(lv1);
  lvect.push_back(lv2);

  return true;
}
