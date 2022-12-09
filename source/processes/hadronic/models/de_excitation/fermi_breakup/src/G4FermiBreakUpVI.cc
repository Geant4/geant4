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
#include "G4FermiFragmentsPoolVI.hh"
#include "G4FermiDecayProbability.hh"
#include "G4FermiChannels.hh"
#include "G4FermiPair.hh"
#include "G4RandomDirection.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4PhysicsModelCatalog.hh"

G4FermiFragmentsPoolVI* G4FermiBreakUpVI::thePool = nullptr;

#ifdef G4MULTITHREADED
G4Mutex G4FermiBreakUpVI::FermiBreakUpVIMutex = G4MUTEX_INITIALIZER;
#endif

G4FermiBreakUpVI::G4FermiBreakUpVI() 
  : theDecay(nullptr), rndmEngine(nullptr), maxZ(9), maxA(17), secID(-1)
{
  frag.reserve(10);
  lvect.reserve(10);
  Z = A = spin = 0;
  secID = G4PhysicsModelCatalog::GetModelID("model_G4FermiBreakUpVI");
  mass = elim = excitation = 0.0;
  tolerance = CLHEP::MeV;  
  frag1 = frag2 = nullptr;
  prob.resize(12,0.0);
  Initialise();
}

G4FermiBreakUpVI::~G4FermiBreakUpVI()
{
  if(G4Threading::IsMasterThread()) { 
    delete thePool;
    thePool = nullptr;
  }
}

void G4FermiBreakUpVI::Initialise()
{
  if(verbose > 1) {
    G4cout << "### G4FermiBreakUpVI::Initialise(): " << thePool << G4endl;
  }
  if(thePool == nullptr) { InitialisePool(); }
  theDecay = thePool->FermiDecayProbability();
  elim = thePool->GetEnergyLimit();
}

void G4FermiBreakUpVI::InitialisePool()
{
#ifdef G4MULTITHREADED
  G4MUTEXLOCK(&G4FermiBreakUpVI::FermiBreakUpVIMutex);
#endif
  if(thePool == nullptr) {
    thePool = new G4FermiFragmentsPoolVI();
  }
#ifdef G4MULTITHREADED
  G4MUTEXUNLOCK(&G4FermiBreakUpVI::FermiBreakUpVIMutex);
#endif
}

G4bool G4FermiBreakUpVI::IsApplicable(G4int ZZ, G4int AA, G4double eexc) const
{
  return (ZZ < maxZ && AA < maxA && AA > 0 && eexc <= elim 
  	  && thePool->HasChannels(ZZ, AA, eexc));
}

void G4FermiBreakUpVI::BreakFragment(G4FragmentVector* theResult, 
				     G4Fragment* theNucleus)
{
  if(verbose > 1) {
    G4cout << "### G4FermiBreakUpVI::BreakFragment start new fragment " 
           << G4endl;
    G4cout << *theNucleus << G4endl;
  }

  // initial fragment
  Z = theNucleus->GetZ_asInt();
  A = theNucleus->GetA_asInt();
  excitation = theNucleus->GetExcitationEnergy();
  mass = theNucleus->GetGroundStateMass() + excitation;
  spin = -1;

  lv0 = theNucleus->GetMomentum();
  rndmEngine = G4Random::getTheEngine();

  // sample first decay of an initial state
  // if not possible to decay - exit
  if(!SampleDecay()) {
    return; 
  }

  G4double time = theNucleus->GetCreationTime();
  delete theNucleus;

  static const G4int imax = 100; 

  // loop over vector of Fermi fragments
  // vector may grow at each iteraction
  for(size_t i=0; i<frag.size(); ++i) {
    Z = frag[i]->GetZ();
    A = frag[i]->GetA();
    spin = frag[i]->GetSpin();
    mass = frag[i]->GetTotalEnergy();
    lv0 = lvect[i];
    if(verbose > 1) {
      G4cout << "# FermiFrag " << i << ".  Z= " << Z << " A= " << A 
	     << " mass= " << mass << " exc= " 
	     << frag[i]->GetExcitationEnergy() << G4endl;
    }
    // stable fragment
    if(!SampleDecay()) {
      if(verbose > 1) { G4cout << "   New G4Fragment" << G4endl; }
      G4Fragment* f = new G4Fragment(A, Z, lv0);
      f->SetSpin(0.5*spin);
      f->SetCreationTime(time);
      f->SetCreatorModelID(secID);
      theResult->push_back(f);
    }
    // limit the loop
    if(i == imax) {
      break;
    }
  }
  frag.clear();
  lvect.clear();
}

G4bool G4FermiBreakUpVI::SampleDecay()
{
  const G4FermiChannels* chan = thePool->ClosestChannels(Z, A, mass);
  if(nullptr == chan) { return false; }
  size_t nn = chan->GetNumberOfChannels();
  if(verbose > 1) {
    G4cout << "== SampleDecay " << nn << " channels Eex= " 
	   << chan->GetExcitation() << G4endl;
  }
  if(0 == nn) { return false; }

  const G4FermiPair* fpair = nullptr;

  // one unstable fragment
  if(1 == nn) {
    fpair = chan->GetPair(0);

    // more pairs
  } else {
    
    // in static probabilities may be used
    if(std::abs(excitation - chan->GetExcitation()) < tolerance) {
      fpair = chan->SamplePair(rndmEngine->flat());
      
    } else {

      // recompute probabilities
      const std::vector<const G4FermiPair*>& pvect = chan->GetChannels();
      if(nn > 12) { prob.resize(nn, 0.0); }
      G4double ptot = 0.0;
      if(verbose > 2) { 
	G4cout << "Start recompute probabilities" << G4endl;
      }
      for(size_t i=0; i<nn; ++i) {
	ptot += theDecay->ComputeProbability(Z, A, -1, mass,
					     pvect[i]->GetFragment1(),
					     pvect[i]->GetFragment2());
	prob[i] = ptot;
	if(verbose > 2) { 
	  G4cout << i << ". " << prob[i]
		 << " Z1= " << pvect[i]->GetFragment1()->GetZ()
		 << " A1= " << pvect[i]->GetFragment1()->GetA()
		 << " Z2= " << pvect[i]->GetFragment2()->GetZ()
		 << " A2= " << pvect[i]->GetFragment2()->GetA()
		 << G4endl;
	}
      }
      ptot *= rndmEngine->flat();
      for(size_t i=0; i<nn; ++i) {
        if(ptot <= prob[i] || i+1 == nn) {
          fpair = pvect[i];
	  break;
	}
      }
    }
  }
  if(!fpair) { return false; }

  frag1 = fpair->GetFragment1();
  frag2 = fpair->GetFragment2();
  
  G4double mass1 = frag1->GetTotalEnergy();
  G4double mass2 = frag2->GetTotalEnergy();
  if(verbose > 2) {
    G4cout << " M= " << mass << " M1= " << mass1 << "  M2= " << mass2 
	   << " Exc1= " << frag1->GetExcitationEnergy() 
	   << " Exc2= " << frag2->GetExcitationEnergy() << G4endl;
  }
  // sample fragment1
  G4double e1 = 0.5*(mass*mass - mass2*mass2 + mass1*mass1)/mass;
  //G4cout << "ekin1= " << e1 - mass1 << G4endl;
  G4double p1(0.0);
  if(e1 > mass1) {
    p1 = std::sqrt((e1 - mass1)*(e1 + mass1));
  } else {
    e1 = mass1;
  }
  G4ThreeVector v = G4RandomDirection();
  G4LorentzVector lv1 = G4LorentzVector(v*p1, e1);

  // compute kinematics
  boostVector = lv0.boostVector();  
  lv1.boost(boostVector);

  lv0 -= lv1;

  G4double e2 = lv0.e();
  if(e2 < mass2) {
    lv0.set(0.,0.,0.,mass2);
  }

  frag.push_back(frag1);
  frag.push_back(frag2);
  lvect.push_back(lv1);
  lvect.push_back(lv0);

  return true;
}

