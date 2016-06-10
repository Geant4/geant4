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
// $Id: G4PhotonEvaporation.cc 70091 2013-05-23 08:55:18Z gcosmo $
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4PhotonEvaporation
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 23 October 1998
//
// Modifications: 
//      
// 8 March 2002, Fan Lei (flei@space.qinetiq.com)
//   
//        Implementation of Internal Convertion process in discrete deexcitation
//        The following public methods have been added. 
//
//       void SetICM (G4bool);
//       void CallFromRDM(G4bool);
//       void SetMaxHalfLife(G4double) ;
//       void SetEOccupancy( G4ElectronOccupancy  eOccupancy) ;
//       G4ElectronOccupancy GetEOccupancy () ;
//
// 11 May 2010, V.Ivanchenko add implementation of EmittedFragment and 
//                           BreakUpFragment methods; cleanup logic
//
// -------------------------------------------------------------------
//

#include "G4PhotonEvaporation.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Gamma.hh"
#include "G4LorentzVector.hh"
#include "G4VGammaTransition.hh"
#include "G4Fragment.hh"
#include "G4FragmentVector.hh"
#include "G4ContinuumGammaDeexcitation.hh"
#include "G4DiscreteGammaDeexcitation.hh"
#include "G4E1Probability.hh"

G4PhotonEvaporation::G4PhotonEvaporation(const G4String & aName,
					 G4EvaporationChannelType timeType)
  : G4VEvaporationChannel(aName, timeType),
   verbose(0), myOwnProbAlgorithm(true),
   eOccupancy(0), vShellNumber(-1), gammaE(0.)
{ 
  probAlgorithm = new G4E1Probability;
  contDeexcitation = new G4ContinuumGammaDeexcitation;
  G4DiscreteGammaDeexcitation* p = new G4DiscreteGammaDeexcitation();
  discrDeexcitation = p;

  p->SetICM(false);

  // Time limits
  G4double timeLimit = DBL_MAX;
  char* env = getenv("G4AddTimeLimitToPhotonEvaporation"); 
  if(env) { timeLimit = 1.e-16*second; }

  // Time for short-cut simulation of photon evaporation
  p->SetTimeLimit(timeLimit);

  // Time limit for isomere production 
  SetMaxHalfLife(1.e-6*second);

  nucleus = 0;
}

G4PhotonEvaporation::~G4PhotonEvaporation()
{ 
  if(myOwnProbAlgorithm) delete probAlgorithm;
  delete discrDeexcitation;
  delete contDeexcitation;
}

G4Fragment* G4PhotonEvaporation::EmittedFragment(G4Fragment* aNucleus)
{
  //G4cout << "G4PhotonEvaporation::EmittedFragment" << G4endl;
  nucleus = aNucleus;
  
  // Do one photon emission by the continues deexcitation  
  contDeexcitation->SetNucleus(nucleus);
  contDeexcitation->Initialize();

  if(contDeexcitation->CanDoTransition()) {  
    G4Fragment* gamma = contDeexcitation->GenerateGamma();
    if(gamma) { 
      if (verbose > 0) {
	G4cout << "G4PhotonEvaporation::EmittedFragment continium deex: "   
	       << gamma << G4endl;
        G4cout << "   Residual: " << nucleus << G4endl;
      }
      return gamma; 
    }
  }

  // Do one photon emission by the discrete deexcitation 
  discrDeexcitation->SetNucleus(nucleus);
  discrDeexcitation->Initialize();

  if(discrDeexcitation->CanDoTransition()) {  
    G4Fragment* gamma = discrDeexcitation->GenerateGamma();
    if(gamma) { 
      if (verbose > 0) {
	G4cout << "G4PhotonEvaporation::EmittedFragment discrete deex: "   
	       << gamma << G4endl;
        G4cout << "   Residual: " << nucleus << G4endl;
      }
      return gamma; 
    }
  }

  if (verbose > 0) {
    G4cout << "G4PhotonEvaporation unable emit gamma: " 
	   << nucleus << G4endl;
  }
  return 0;
}

G4FragmentVector* G4PhotonEvaporation::BreakUpFragment(G4Fragment* aNucleus)
{
  //G4cout << "G4PhotonEvaporation::BreakUpFragment" << G4endl;
  // The same pointer of primary nucleus
  nucleus = aNucleus;
  contDeexcitation->SetNucleus(nucleus);
  discrDeexcitation->SetNucleus(nucleus);

  // Do the whole gamma chain 
  G4FragmentVector* products = contDeexcitation->DoChain();  
  if( !products ) { products = new G4FragmentVector(); }

  if (verbose > 0) {
    G4cout << "G4PhotonEvaporation::BreakUpFragment " << products->size() 
	   << " gammas from ContinuumDeexcitation " << G4endl;
    G4cout << "   Residual: " << nucleus << G4endl;
  }
  // Products from discrete gamma transitions
  G4FragmentVector* discrProducts = discrDeexcitation->DoChain();
  if(discrProducts) {
    eOccupancy = discrDeexcitation->GetEO();
    vShellNumber = discrDeexcitation->GetVacantSN();

    // not sure if the following line is needed!
    discrDeexcitation->SetVaccantSN(-1);

    if (verbose > 0) {
      G4cout << "G4PhotonEvaporation::BreakUpFragment " << discrProducts->size() 
	     << " gammas from DiscreteDeexcitation " << G4endl;
      G4cout << "   Residual: " << nucleus << G4endl;
    }
    G4FragmentVector::iterator i;
    for (i = discrProducts->begin(); i != discrProducts->end(); ++i)
      {
	products->push_back(*i);
      }
    delete discrProducts;
  }

  if (verbose > 0) {
    G4cout << "*-*-* Photon evaporation: " << products->size() << G4endl;
  }
  return products;
}

G4FragmentVector* G4PhotonEvaporation::BreakUp(const G4Fragment& aNucleus)
{
  //G4cout << "G4PhotonEvaporation::BreakUp" << G4endl;
  nucleus = new G4Fragment(aNucleus);

  contDeexcitation->SetNucleus(nucleus);
  discrDeexcitation->SetNucleus(nucleus);
  
  // Do one photon emission

  // Products from continuum gamma transitions

  G4FragmentVector* products = contDeexcitation->DoTransition();  
  if( !products ) { products = new G4FragmentVector(); }
  else if(verbose > 0) {
    G4cout << "G4PhotonEvaporation::BreakUp " << products->size() 
	   << " gammas from ContinuesDeexcitation " << G4endl;
    G4cout << "   Residual: " << nucleus << G4endl;
  }

  if (0 == products->size())
    {
      // Products from discrete gamma transitions
      G4FragmentVector* discrProducts = discrDeexcitation->DoTransition();

      if (discrProducts) {
	eOccupancy = discrDeexcitation->GetEO();
	vShellNumber = discrDeexcitation->GetVacantSN();

	// not sure if the following line is needed!
	discrDeexcitation->SetVaccantSN(-1);
	//
	if (verbose > 0) {
	  G4cout << " = BreakUp = " << discrProducts->size() 
		 << " gammas from DiscreteDeexcitation " 
		 << G4endl;
	  G4cout << "   Residual: " << nucleus << G4endl;
	}
	G4FragmentVector::iterator i;
	for (i = discrProducts->begin(); i != discrProducts->end(); ++i)
	  {
	    products->push_back(*i);
	  }
	delete discrProducts;
      }
    }
  
  // Add deexcited nucleus to products
  products->push_back(nucleus);

  if (verbose > 0) {
    G4cout << "*-*-*-* Photon evaporation: " << products->size() << G4endl;
  }

  return products;
}

G4FragmentVector* G4PhotonEvaporation::BreakItUp(const G4Fragment& aNucleus)
{
  // The same pointer of primary nucleus
  nucleus = new G4Fragment(aNucleus);
  contDeexcitation->SetNucleus(nucleus);
  discrDeexcitation->SetNucleus(nucleus);

  //G4cout << "G4PhotonEvaporation::BreakItUp:  " << nucleus << G4endl;

  // Do the whole gamma chain 
  G4FragmentVector* products = contDeexcitation->DoChain();  
  if( !products ) { products = new G4FragmentVector; }

  // Products from continuum gamma transitions
  if (verbose > 0) {
    G4cout << " = BreakItUp = " << products->size()
	   << " gammas from ContinuumDeexcitation " << G4endl;
  }

  // Products from discrete gamma transitions
  G4FragmentVector* discrProducts = discrDeexcitation->DoChain();
  if(discrProducts) {
    eOccupancy = discrDeexcitation->GetEO();
    vShellNumber = discrDeexcitation->GetVacantSN();

    // not sure if the following line is needed!
    discrDeexcitation->SetVaccantSN(-1);

    if (verbose > 0) {
      G4cout << " = BreakItUp = " << discrProducts->size() 
	     << " gammas from DiscreteDeexcitation " << G4endl;
    }
    G4FragmentVector::iterator i;
    for (i = discrProducts->begin(); i != discrProducts->end(); ++i)
      {
	products->push_back(*i);
      }
    delete discrProducts;
  }
  // Add deexcited nucleus to products
  products->push_back(nucleus);

  if (verbose > 0) {
    G4cout << "*-*-* Photon evaporation: " << products->size() << G4endl;
  }
  return products;
}

G4double 
G4PhotonEvaporation::GetEmissionProbability(G4Fragment* theNucleus)
{
  nucleus = theNucleus;
  G4double prob = 
    probAlgorithm->EmissionProbability(*nucleus, nucleus->GetExcitationEnergy());
  return prob;
}

void 
G4PhotonEvaporation::SetEmissionStrategy(G4VEmissionProbability * alg)
{
  // CD - not sure about always wanting to delete this pointer....
  if(myOwnProbAlgorithm) { delete probAlgorithm; }
  probAlgorithm = alg;
  myOwnProbAlgorithm = false;
}


void G4PhotonEvaporation::SetVerboseLevel(G4int verb)
{
  verbose = verb;
  contDeexcitation->SetVerboseLevel(verbose);
  discrDeexcitation->SetVerboseLevel(verbose);
}

void G4PhotonEvaporation::SetICM(G4bool ic)
{
  static_cast<G4DiscreteGammaDeexcitation*>(discrDeexcitation)->SetICM(ic);
}

void G4PhotonEvaporation::SetMaxHalfLife(G4double hl)
{
  static_cast<G4DiscreteGammaDeexcitation*>(discrDeexcitation)->SetHL(hl);
}

void G4PhotonEvaporation::SetTimeLimit(G4double val)
{
  discrDeexcitation->SetTimeLimit(val);
}

void G4PhotonEvaporation::RDMForced(G4bool fromRDM)
{
  static_cast<G4DiscreteGammaDeexcitation*>(discrDeexcitation)->SetRDM(fromRDM);
}

void G4PhotonEvaporation::SetEOccupancy(G4ElectronOccupancy eo)
{
  discrDeexcitation->SetEO(eo);
}




