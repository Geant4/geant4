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
// $Id: G4PhotonEvaporationOLD.cc 91873 2015-08-07 17:50:24Z vnivanch $
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4PhotonEvaporationOLD
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
// 22 October 2015, V.Ivanchenko renamed to G4PhotonEvaporationOLD
//
// -------------------------------------------------------------------
//

#include "G4PhotonEvaporationOLD.hh"

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
#include "G4NuclearLevelStore.hh"

static const G4double tolerance = 0.1*CLHEP::keV;

G4PhotonEvaporationOLD::G4PhotonEvaporationOLD(const G4String & aName)
  : G4VEvaporationChannel(aName),
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

  // Time for short-cut simulation of photon evaporationOLD
  p->SetTimeLimit(timeLimit);

  // Default time limit for isomere production 
  SetMaxHalfLife(DBL_MAX);
  //  SetMaxHalfLife(1.e-6*second);

  nucleus = 0;
}

G4PhotonEvaporationOLD::~G4PhotonEvaporationOLD()
{ 
  if(myOwnProbAlgorithm) delete probAlgorithm;
  delete discrDeexcitation;
  delete contDeexcitation;
}

G4Fragment* G4PhotonEvaporationOLD::EmittedFragment(G4Fragment* aNucleus)
{
  //G4cout << "G4PhotonEvaporationOLD::EmittedFragment" << G4endl;
  vShellNumber = -1;
  G4Fragment* gamma = 0;
  if(contDeexcitation->CanDoTransition(aNucleus)) { 
    gamma = contDeexcitation->GenerateGamma(aNucleus);
    if(gamma && verbose > 1) {
      G4cout << "G4PhotonEvaporationOLD::EmittedFragment continium deex: "   
	     << gamma << G4endl;
      G4cout << "   Residual: " << aNucleus << G4endl;
    }
  } else if(discrDeexcitation->CanDoTransition(aNucleus)) {

    // Do one photon emission by the discrete deexcitation 
    gamma = discrDeexcitation->GenerateGamma(aNucleus);
    if(gamma) { 
      vShellNumber = discrDeexcitation->GetVacantSN();
      if (verbose > 1) {
	G4cout << "G4PhotonEvaporationOLD::EmittedFragment discrete deex: "   
	       << gamma << G4endl;
        G4cout << "   Residual: " << aNucleus << G4endl;
      }
    }
  }
  return gamma; 
}

G4bool G4PhotonEvaporationOLD::BreakUpChain(G4FragmentVector* products,
					 G4Fragment* aNucleus)
{
  //G4cout << "G4PhotonEvaporationOLD::BreakUpChain" << G4endl;
  G4Fragment* gamma = 0;
  // one continues emission is not excluded
  if(contDeexcitation->CanDoTransition(aNucleus)) { 
    gamma = contDeexcitation->GenerateGamma(aNucleus);
    if(gamma) { 
      if (verbose > 1) {
	G4cout << "G4PhotonEvaporationOLD::EmittedFragment continium deex: "   
	       << gamma << G4endl;
	G4cout << "   Residual: " << aNucleus << G4endl;
      }
      products->push_back(gamma);
    }
  }
  // main emissions are discrete
  discrDeexcitation->DoChain(products, aNucleus);
  return false;
}

G4FragmentVector* G4PhotonEvaporationOLD::BreakUpFragment(G4Fragment* aNucleus)
{
  //G4cout << "G4PhotonEvaporationOLD::BreakUpFragment" << G4endl;
  G4FragmentVector* products = new G4FragmentVector();
  BreakUpChain(products, aNucleus);
  return products;
}

G4FragmentVector* G4PhotonEvaporationOLD::BreakUp(const G4Fragment& theNucleus)
{
  //G4cout << "G4PhotonEvaporationOLD::BreakUp" << G4endl;
  G4Fragment* aNucleus = new G4Fragment(theNucleus);
  G4FragmentVector* products = new G4FragmentVector();
  //discrDeexcitation->DoChain(products, aNucleus);
  BreakUpChain(products, aNucleus);
  products->push_back(aNucleus);
  return products;
}

G4FragmentVector* G4PhotonEvaporationOLD::BreakItUp(const G4Fragment& theNucleus)
{
  //G4cout << "G4PhotonEvaporationOLD::BreakItUp" << G4endl;
  G4Fragment* aNucleus = new G4Fragment(theNucleus);
  G4FragmentVector* products = new G4FragmentVector();
  BreakUpChain(products, aNucleus);
  products->push_back(aNucleus);
  return products;
}

G4double 
G4PhotonEvaporationOLD::GetEmissionProbability(G4Fragment* theNucleus)
{
  G4double prob = 0.0;
  G4int Z = theNucleus->GetZ_asInt();
  G4int A = theNucleus->GetA_asInt();
  G4double eexc = theNucleus->GetExcitationEnergy();
  if(0 < Z && Z < A && eexc > tolerance) {
    prob = probAlgorithm->EmissionProbability(*theNucleus, eexc);
  }
  return prob;
}

void 
G4PhotonEvaporationOLD::SetEmissionStrategy(G4VEmissionProbability * alg)
{
  if(myOwnProbAlgorithm) { delete probAlgorithm; }
  probAlgorithm = alg;
  myOwnProbAlgorithm = false;
}

void G4PhotonEvaporationOLD::SetVerboseLevel(G4int verb)
{
  verbose = verb;
  contDeexcitation->SetVerboseLevel(verbose);
  discrDeexcitation->SetVerboseLevel(verbose);
}

void G4PhotonEvaporationOLD::SetICM(G4bool ic)
{
  static_cast<G4DiscreteGammaDeexcitation*>(discrDeexcitation)->SetICM(ic);
}

void G4PhotonEvaporationOLD::SetMaxHalfLife(G4double hl)
{
  static_cast<G4DiscreteGammaDeexcitation*>(discrDeexcitation)->SetHL(hl);
}

void G4PhotonEvaporationOLD::SetTimeLimit(G4double val)
{
  discrDeexcitation->SetTimeLimit(val);
}

void G4PhotonEvaporationOLD::RDMForced(G4bool fromRDM)
{
  static_cast<G4DiscreteGammaDeexcitation*>(discrDeexcitation)->SetRDM(fromRDM);
}

void G4PhotonEvaporationOLD::SetEOccupancy(G4ElectronOccupancy eo)
{
  discrDeexcitation->SetEO(eo);
}




