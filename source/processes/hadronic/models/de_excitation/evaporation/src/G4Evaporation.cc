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
// $Id$
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

#include "G4Evaporation.hh"
#include "G4SystemOfUnits.hh"
#include "G4EvaporationFactory.hh"
#include "G4EvaporationGEMFactory.hh"
#include "G4EvaporationDefaultGEMFactory.hh"
#include "G4HadronicException.hh"
#include "G4NistManager.hh"
#include "G4FermiFragmentsPool.hh"
#include "G4PhotonEvaporation.hh"

G4Evaporation::G4Evaporation()  
  : theChannels(0),nChannels(0)
{
  SetPhotonEvaporation(new G4PhotonEvaporation());
  theChannelFactory = new G4EvaporationDefaultGEMFactory(thePhotonEvaporation);
  SetParameters();
  InitialiseEvaporation();
}

G4Evaporation::G4Evaporation(G4VEvaporationChannel* photoEvaporation)  
  : theChannels(0),nChannels(0)
{
  if(photoEvaporation) { SetPhotonEvaporation(photoEvaporation); }
  else                 { SetPhotonEvaporation(new G4PhotonEvaporation()); }

  theChannelFactory = new G4EvaporationDefaultGEMFactory(thePhotonEvaporation);
  SetParameters();
  InitialiseEvaporation();
}

/*
G4Evaporation::G4Evaporation(std::vector<G4VEvaporationChannel*>* channels) 
  : theChannels(channels),theChannelFactory(0),nChannels(0)
{
  // are input relible?
  G4bool accepted = true;
  if(!theChannels) { accepted = false; }
  else if(0 == theChannels->size()) { accepted = false; }

  if(accepted) {
    SetPhotonEvaporation((*channels)[0]); 
  } else {
    SetPhotonEvaporation(new G4PhotonEvaporation());
    theChannelFactory = new G4EvaporationDefaultGEMFactory(thePhotonEvaporation); 
  }
  SetParameters();
  InitialiseEvaporation();
}
*/

G4Evaporation::~G4Evaporation()
{
  CleanChannels();
  delete thePhotonEvaporation;
  delete theChannelFactory; 
}

void G4Evaporation::CleanChannels()
{
  for (size_t i=1; i<nChannels; ++i) { delete (*theChannels)[i]; }
  delete theChannels;
}

void G4Evaporation::SetParameters()
{
  nist = G4NistManager::Instance();
  minExcitation = CLHEP::keV;
  maxZforFBU = G4FermiFragmentsPool::Instance()->GetMaxZ();
  maxAforFBU = G4FermiFragmentsPool::Instance()->GetMaxA();
  probabilities.reserve(68);
}

void G4Evaporation::InitialiseEvaporation()
{
  CleanChannels();
  theChannels = theChannelFactory->GetChannel(); 
  nChannels = theChannels->size();   
  probabilities.resize(nChannels, 0.0);
  Initialise();
}

void G4Evaporation::Initialise()
{
  for(size_t i=0; i<nChannels; ++i) {
    (*theChannels)[i]->SetOPTxs(OPTxs);
    (*theChannels)[i]->UseSICB(useSICB);
  }
}

void G4Evaporation::SetDefaultChannel()
{
  delete theChannelFactory;
  theChannelFactory = new G4EvaporationFactory(thePhotonEvaporation);
  InitialiseEvaporation();
}

void G4Evaporation::SetGEMChannel()
{
  delete theChannelFactory;
  theChannelFactory = new G4EvaporationGEMFactory(thePhotonEvaporation);
  InitialiseEvaporation();
}

void G4Evaporation::SetCombinedChannel()
{
  delete theChannelFactory;
  theChannelFactory = new G4EvaporationDefaultGEMFactory(thePhotonEvaporation);
  InitialiseEvaporation();
}

void G4Evaporation::SetPhotonEvaporation(G4VEvaporationChannel* ptr)
{
  if(ptr) { G4VEvaporation::SetPhotonEvaporation(ptr); }
  if(0 < nChannels) { (*theChannels)[0] = ptr; }
}

G4FragmentVector * G4Evaporation::BreakItUp(const G4Fragment &theNucleus)
{
  G4FragmentVector * theResult = new G4FragmentVector;
  G4FragmentVector * theTempResult;
  const G4double Elimit = 3*MeV;

  // The residual nucleus (after evaporation of each fragment)
  G4Fragment* theResidualNucleus = new G4Fragment(theNucleus);

  G4double totprob, prob, oldprob = 0.0;
  size_t maxchannel, i;

  G4int Amax = theResidualNucleus->GetA_asInt();

  // Starts loop over evaporated particles, loop is limited by number
  // of nucleons
  for(G4int ia=0; ia<Amax; ++ia) {
 
    // g,n,p and light fragments - evaporation is finished
    G4int Z = theResidualNucleus->GetZ_asInt();
    G4int A = theResidualNucleus->GetA_asInt();

    // stop deecitation loop if can be deexcited by FBU
    if(maxZforFBU > Z && maxAforFBU >= A) {
      theResult->push_back(theResidualNucleus);
      return theResult;
    }

    // check if it is stable, then finish evaporation
    G4double abun = nist->GetIsotopeAbundance(Z, A); 
    /*
    G4cout << "### G4Evaporation::BreakItUp step " << ia << " Z= " << Z
    	   << " A= " << A << " Eex(MeV)= " 
    	   << theResidualNucleus->GetExcitationEnergy()
    	   << " aban= " << abun << G4endl;
    */
    // stop deecitation loop in the case of the cold stable fragment 
    G4double Eex = theResidualNucleus->GetExcitationEnergy();
    if(Eex <= minExcitation && abun > 0.0) {
      theResult->push_back(theResidualNucleus);
      return theResult;
    }
 
    totprob = 0.0;
    maxchannel = nChannels;
    /*
    G4cout << "### Evaporation loop #" << ia 
    	   << "  Fragment: " << theResidualNucleus << G4endl;
    */
    // loop over evaporation channels
    for(i=0; i<nChannels; ++i) {
      prob = (*theChannels)[i]->GetEmissionProbability(theResidualNucleus);
      //G4cout << "  Channel# " << i << "  prob= " << prob << G4endl; 

      totprob += prob;
      probabilities[i] = totprob;
      // if two recent probabilities are near zero stop computations
      if(i>=8) {
	if(prob <= totprob*1.e-8 && oldprob <= totprob*1.e-8) {
	  maxchannel = i+1; 
	  break;
	}
      }
      oldprob = prob;
      // protection for very excited fragment - avoid GEM
      if(7 == i && Eex > Elimit*A) {
        maxchannel = 8;
        break;
      }
    }

    // photon evaporation in the case of no other channels available
    // do evaporation chain and reset total probability
    if(0.0 < totprob && probabilities[0] == totprob) {
      //G4cout << "Start gamma evaporation" << G4endl;
      theTempResult = (*theChannels)[0]->BreakUpFragment(theResidualNucleus);
      if(theTempResult) {
	size_t nsec = theTempResult->size();
	for(size_t j=0; j<nsec; ++j) {
	  theResult->push_back((*theTempResult)[j]);
	}
        delete theTempResult;
      }
      totprob = 0.0;
    }

    // stable fragnent - evaporation is finished
    if(0.0 == totprob) {

      // if fragment is exotic, then try to decay it
      if(0.0 == abun && Z < 20) {
        //G4cout << "$$$ Decay exotic fragment" << G4endl;
	theTempResult = unstableBreakUp.BreakUpFragment(theResidualNucleus);
	if(theTempResult) {
	  size_t nsec = theTempResult->size();
	  for(size_t j=0; j<nsec; ++j) {
	    theResult->push_back((*theTempResult)[j]);
	  }
	  delete theTempResult;
	}
      }

      // save residual fragment
      theResult->push_back(theResidualNucleus);
      return theResult;
    }


    // select channel
    totprob *= G4UniformRand();
    // loop over evaporation channels
    for(i=0; i<maxchannel; ++i) { if(probabilities[i] >= totprob) { break; } }

    // this should not happen
    if(i >= nChannels) { i = nChannels - 1; }


    // single photon evaporation, primary pointer is kept
    if(0 == i) {
      //G4cout << "Single gamma" << G4endl;
      G4Fragment* gamma = (*theChannels)[0]->EmittedFragment(theResidualNucleus);
      if(gamma) { theResult->push_back(gamma); }

      // fission, return results to the main loop if fission is succesful
    } else if(1 == i) {
      //G4cout << "Fission" << G4endl;
      theTempResult = (*theChannels)[1]->BreakUp(*theResidualNucleus);
      if(theTempResult) {
	size_t nsec = theTempResult->size();
        G4bool deletePrimary = true;
	for(size_t j=0; j<nsec; ++j) {
          if(theResidualNucleus == (*theTempResult)[j]) { deletePrimary = false; }
	  theResult->push_back((*theTempResult)[j]);
	}
        if(deletePrimary) { delete theResidualNucleus; }
        delete theTempResult;
	return theResult;
      }

      // other channels
    } else {
      //G4cout << "Channel # " << i << G4endl;
      theTempResult = (*theChannels)[i]->BreakUp(*theResidualNucleus);
      if(theTempResult) {
	size_t nsec = theTempResult->size();
        if(nsec > 0) {
          --nsec;
	  for(size_t j=0; j<nsec; ++j) {
	    theResult->push_back((*theTempResult)[j]);
	  }
	  // if the residual change its pointer 
	  // then delete previous residual fragment and update to the new
          if(theResidualNucleus != (*theTempResult)[nsec] ) { 
	    delete theResidualNucleus; 
	    theResidualNucleus = (*theTempResult)[nsec];
	  }
	}
        delete theTempResult;
      }
    }
  }

  // loop is stopped, save residual
  theResult->push_back(theResidualNucleus);
  
  return theResult;
}

