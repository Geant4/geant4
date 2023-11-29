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
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// Alex Howard - added protection for negative probabilities in the sum, 14/2/07
//
// Modif (03 September 2008) by J. M. Quesada for external choice of inverse 
// cross section option
// JMQ (06 September 2008) Also external choices have been added for 
// superimposed Coulomb barrier (if useSICBis set true, by default is false) 
//
// V.Ivanchenko (27 July 2009)  added G4EvaporationDefaultGEMFactory option
// V.Ivanchenko (10 May  2010)  rewrited BreakItUp method: do not make new/delete
//                              photon channel first, fission second,
//                              added G4UnstableFragmentBreakUp to decay 
//                              unphysical fragments (like 2n or 2p),
//                              use Z and A integer
// V.Ivanchenko (22 April 2011) added check if a fragment can be deexcited 
//                              by the FermiBreakUp model
// V.Ivanchenko (23 January 2012) added pointer of G4VPhotonEvaporation 
// V.Ivanchenko (6 May 2013)    added check of existence of residual ion
//                              in the ion table

#include "G4Evaporation.hh"
#include "G4SystemOfUnits.hh"
#include "G4EvaporationFactory.hh"
#include "G4EvaporationGEMFactory.hh"
#include "G4EvaporationGEMFactoryVI.hh"
#include "G4EvaporationDefaultGEMFactory.hh"
#include "G4NistManager.hh"
#include "G4VFermiBreakUp.hh"
#include "G4PhotonEvaporation.hh"
#include "G4VEvaporationChannel.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4NuclearLevelData.hh"
#include "G4LevelManager.hh"
#include "G4UnstableFragmentBreakUp.hh"
#include "Randomize.hh"

G4Evaporation::G4Evaporation(G4VEvaporationChannel* photoEvaporation)  
  : G4VEvaporation(),fVerbose(0),nChannels(0),minExcitation(0.1*keV),
    isInitialised(false)
{
  if(photoEvaporation) { SetPhotonEvaporation(photoEvaporation); }
  else                 { SetPhotonEvaporation(new G4PhotonEvaporation()); }

  channelType = fDummy;
  theChannelFactory = nullptr;

  fLevelData = G4NuclearLevelData::GetInstance();
  theTableOfIons = G4ParticleTable::GetParticleTable()->GetIonTable();
  nist = G4NistManager::Instance();
  unstableBreakUp = new G4UnstableFragmentBreakUp();
  /*
    G4cout << "G4Evaporation() " << this << " thePhotonEvaporation: "
	   << photoEvaporation << "  UnstableFragmentBreakUp: " 
           << unstableBreakUp << G4endl;
  */
}

G4Evaporation::~G4Evaporation()
{
  delete unstableBreakUp;
}

void G4Evaporation::InitialiseChannels()
{
  if(isInitialised) { return; }

  G4DeexPrecoParameters* param = fLevelData->GetParameters(); 
  minExcitation = param->GetMinExcitation();
  fVerbose = param->GetVerbose();
  unstableBreakUp->SetVerbose(fVerbose);

  G4DeexChannelType type = param->GetDeexChannelsType();
  if(type == fCombined) { SetCombinedChannel(); }
  else if(type == fGEM) { SetGEMChannel(); }
  else if(type == fEvaporation) { SetDefaultChannel(); }
  else if(type == fGEMVI) { SetGEMVIChannel(); }

  isInitialised = true;
}

void G4Evaporation::InitialiseChannelFactory()
{
  theChannels = theChannelFactory->GetChannel(); 
  nChannels = theChannels->size();   
  probabilities.resize(nChannels, 0.0);

  if(fVerbose > 1) {
    G4cout << "### G4Evaporation::InitialiseChannelFactory  for " 
	   << nChannels << " channels " << this << G4endl;
  }
  for(size_t i=0; i<nChannels; ++i) {
    (*theChannels)[i]->SetOPTxs(OPTxs);
    (*theChannels)[i]->Initialise();
  }
}

void G4Evaporation::SetDefaultChannel()
{
  if(fEvaporation != channelType) {
    channelType = fEvaporation;
    if(theChannelFactory) {
      CleanChannels();
      delete theChannelFactory;
    }
    theChannelFactory = new G4EvaporationFactory(thePhotonEvaporation);
    InitialiseChannelFactory();
  }
}

void G4Evaporation::SetGEMChannel()
{
  if(fGEM != channelType) {
    channelType = fCombined;
    if(theChannelFactory) {
      CleanChannels();
      delete theChannelFactory;
    }
    theChannelFactory = new G4EvaporationGEMFactory(thePhotonEvaporation);
    InitialiseChannelFactory();
  }
}

void G4Evaporation::SetGEMVIChannel()
{
  if(fGEMVI != channelType) {
    channelType = fGEMVI;
    if(theChannelFactory) {
      CleanChannels();
      delete theChannelFactory;
    }
    theChannelFactory = new G4EvaporationGEMFactoryVI(thePhotonEvaporation);
    InitialiseChannelFactory();
  }
}

void G4Evaporation::SetCombinedChannel()
{
  if(fCombined != channelType) {
    channelType = fCombined;
    if(theChannelFactory) {
      CleanChannels();
      delete theChannelFactory;
    }
    theChannelFactory = 
      new G4EvaporationDefaultGEMFactory(thePhotonEvaporation);
    InitialiseChannelFactory();
  }
}

void G4Evaporation::BreakFragment(G4FragmentVector* theResult, 
				  G4Fragment* theResidualNucleus)
{
  if(!isInitialised) { InitialiseChannels(); }

  G4double totprob, prob, oldprob = 0.0;
  size_t maxchannel, i;

  G4int Amax = theResidualNucleus->GetA_asInt();
  if(fVerbose > 1) {
    G4cout << "### G4Evaporation::BreakItUp loop" << G4endl;
  }
  CLHEP::HepRandomEngine* rndm = G4Random::getTheEngine();

  // Starts loop over evaporated particles, loop is limited by number
  // of nucleons
  for(G4int ia=0; ia<Amax; ++ia) {
 
    // g,n,p and light fragments - evaporation is finished
    G4int Z = theResidualNucleus->GetZ_asInt();
    G4int A = theResidualNucleus->GetA_asInt();
    if(A <= 1) { break; }
    G4double Eex = theResidualNucleus->GetExcitationEnergy();

    // stop deecitation loop if residual can be deexcited by FBU    
    if(theFBU->IsApplicable(Z, A, Eex)) { break; }

    // check if it is stable, then finish evaporation
    G4double abun = nist->GetIsotopeAbundance(Z, A); 
    // stop deecitation loop in the case of a cold stable fragment 
    if(Eex <= minExcitation && 
       (abun > 0.0 || (A == 3 && (Z == 1 || Z == 2)))) { break; }
 
    totprob = 0.0;
    maxchannel = nChannels;
    if(fVerbose > 1) {
      G4cout << "Evaporation# " << ia << " Z= " << Z << " A= " << A 
             << " Eex(MeV)= " << theResidualNucleus->GetExcitationEnergy()
	     << " aban= " << abun << G4endl;
    }
    // loop over evaporation channels
    for(i=0; i<nChannels; ++i) {
      prob = (*theChannels)[i]->GetEmissionProbability(theResidualNucleus);
      if(fVerbose > 1 && prob > 0.0) {
	G4cout << "    Channel# " << i << "  prob= " << prob << G4endl; 
      }
      totprob += prob;
      probabilities[i] = totprob;

      // if two recent probabilities are near zero stop computations
      if(i>=8 && prob > 0.0) {
	if(prob <= totprob*1.e-8 && oldprob <= totprob*1.e-8) {
	  maxchannel = i+1; 
	  break;
	}
      }
      oldprob = prob;
    }

    // photon evaporation in the case of no other channels available
    // do evaporation chain and return back ground state fragment
    if(0.0 < totprob && probabilities[0] == totprob) {
      if(fVerbose > 1) {
	G4cout << "$$$ Start chain of gamma evaporation" << G4endl;
      }
      (*theChannels)[0]->BreakUpChain(theResult, theResidualNucleus);

      // release residual stable fragment 
      if(abun > 0.0) {
	theResidualNucleus->SetLongLived(true);
	break;
      }
      // release residual fragment known to FBU
      Eex = theResidualNucleus->GetExcitationEnergy();
      if(theFBU->IsApplicable(Z, A, Eex)) { break; }

      // release residual fragment with non-zero life time
      if(theResidualNucleus->IsLongLived()) { break; }
      totprob = 0.0;
    }

    if(0.0 == totprob && A < 30) {
      // if residual fragment is exotic, then it forced to decay 
      // if success, then decay product is added to results 
      if(fVerbose > 1) { 
	G4cout << "$$$ Decay exotic fragment" << G4endl; 
      }
      if(unstableBreakUp->BreakUpChain(theResult, theResidualNucleus)) {
        continue;
      }
      // release if it is not possible to decay
      break;
    }

    // select channel
    totprob *= rndm->flat();

    // loop over evaporation channels
    for(i=0; i<maxchannel; ++i) { if(probabilities[i] >= totprob) { break; } }

    if(fVerbose > 1) { G4cout << "$$$ Channel # " << i << G4endl; }
    G4Fragment* frag = (*theChannels)[i]->EmittedFragment(theResidualNucleus);
    if(fVerbose > 2 && frag) { G4cout << "   " << *frag << G4endl; }

    // normaly a fragment should be created
    if(nullptr != frag) { theResult->push_back(frag); }
    else { break; }
  }
}
