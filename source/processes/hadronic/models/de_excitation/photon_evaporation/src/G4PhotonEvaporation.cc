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
// $Id: G4PhotonEvaporation.cc,v 1.16 2010-11-23 18:03:21 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "Randomize.hh"
#include "G4Gamma.hh"
#include "G4LorentzVector.hh"
#include "G4VGammaTransition.hh"
#include "G4Fragment.hh"
#include "G4FragmentVector.hh"
#include "G4ContinuumGammaDeexcitation.hh"
#include "G4DiscreteGammaDeexcitation.hh"
#include "G4E1Probability.hh"

G4PhotonEvaporation::G4PhotonEvaporation()
  :_verbose(0),_myOwnProbAlgorithm (true),
   _eOccupancy(0), _vShellNumber(-1),_gammaE(0.)
{ 
  _probAlgorithm = new G4E1Probability;
  G4DiscreteGammaDeexcitation* p = new G4DiscreteGammaDeexcitation();
  p->SetICM(false);
  _discrDeexcitation = p;
  _contDeexcitation = new G4ContinuumGammaDeexcitation;
  _nucleus = 0;
}

G4PhotonEvaporation::~G4PhotonEvaporation()
{ 
  if(_myOwnProbAlgorithm) delete _probAlgorithm;
  delete _discrDeexcitation;
  delete _contDeexcitation;
}

void G4PhotonEvaporation::Initialize(const G4Fragment& fragment)
{
  _nucleus = const_cast<G4Fragment*>(&fragment);
}

G4Fragment* G4PhotonEvaporation::EmittedFragment(G4Fragment* nucleus)
{
  //G4cout << "G4PhotonEvaporation::EmittedFragment" << G4endl;
  _nucleus = nucleus;
  
  // Do one photon emission by the continues deexcitation  
  _contDeexcitation->SetNucleus(_nucleus);
  _contDeexcitation->Initialize();

  if(_contDeexcitation->CanDoTransition()) {  
    G4Fragment* gamma = _contDeexcitation->GenerateGamma();
    if(gamma) { 
      if (_verbose > 0) {
	G4cout << "G4PhotonEvaporation::EmittedFragment continium deex: "   
	       << gamma << G4endl;
        G4cout << "   Residual: " << nucleus << G4endl;
      }
      return gamma; 
    }
  }

  // Do one photon emission by the discrete deexcitation 
  _discrDeexcitation->SetNucleus(_nucleus);
  _discrDeexcitation->Initialize();

  if(_discrDeexcitation->CanDoTransition()) {  
    G4Fragment* gamma = _discrDeexcitation->GenerateGamma();
    if(gamma) { 
      if (_verbose > 0) {
	G4cout << "G4PhotonEvaporation::EmittedFragment discrete deex: "   
	       << gamma << G4endl;
        G4cout << "   Residual: " << nucleus << G4endl;
      }
      return gamma; 
    }
  }

  if (_verbose > 0) {
    G4cout << "G4PhotonEvaporation unable emit gamma: " 
	   << nucleus << G4endl;
  }
  return 0;
}

G4FragmentVector* G4PhotonEvaporation::BreakUpFragment(G4Fragment* nucleus)
{
  //G4cout << "G4PhotonEvaporation::BreakUpFragment" << G4endl;
  // The same pointer of primary nucleus
  _nucleus = nucleus;
  _contDeexcitation->SetNucleus(_nucleus);
  _discrDeexcitation->SetNucleus(_nucleus);

  // Do the whole gamma chain 
  G4FragmentVector* products = _contDeexcitation->DoChain();  
  if( !products ) { products = new G4FragmentVector(); }

  if (_verbose > 0) {
    G4cout << "G4PhotonEvaporation::BreakUpFragment " << products->size() 
	   << " gammas from ContinuumDeexcitation " << G4endl;
    G4cout << "   Residual: " << nucleus << G4endl;
  }
  // Products from discrete gamma transitions
  G4FragmentVector* discrProducts = _discrDeexcitation->DoChain();
  if(discrProducts) {
    _eOccupancy = _discrDeexcitation->GetEO();
    _vShellNumber = _discrDeexcitation->GetVacantSN();

    // not sure if the following line is needed!
    _discrDeexcitation->SetVaccantSN(-1);

    if (_verbose > 0) {
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

  if (_verbose > 0) {
    G4cout << "*-*-* Photon evaporation: " << products->size() << G4endl;
  }
  return products;
}

G4FragmentVector* G4PhotonEvaporation::BreakUp(const G4Fragment& nucleus)
{
  //G4cout << "G4PhotonEvaporation::BreakUp" << G4endl;
  _nucleus = new G4Fragment(nucleus);

  _contDeexcitation->SetNucleus(_nucleus);
  _discrDeexcitation->SetNucleus(_nucleus);
  
  // Do one photon emission

  // Products from continuum gamma transitions

  G4FragmentVector* products = _contDeexcitation->DoTransition();  
  if( !products ) { products = new G4FragmentVector(); }
  else if(_verbose > 0) {
    G4cout << "G4PhotonEvaporation::BreakUp " << products->size() 
	   << " gammas from ContinuesDeexcitation " << G4endl;
    G4cout << "   Residual: " << nucleus << G4endl;
  }

  if (0 == products->size())
    {
      // Products from discrete gamma transitions
      G4FragmentVector* discrProducts = _discrDeexcitation->DoTransition();

      if (discrProducts) {
	_eOccupancy = _discrDeexcitation->GetEO();
	_vShellNumber = _discrDeexcitation->GetVacantSN();

	// not sure if the following line is needed!
	_discrDeexcitation->SetVaccantSN(-1);
	//
	if (_verbose > 0) {
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
  products->push_back(_nucleus);

  if (_verbose > 0) {
    G4cout << "*-*-*-* Photon evaporation: " << products->size() << G4endl;
  }

  return products;
}

G4FragmentVector* G4PhotonEvaporation::BreakItUp(const G4Fragment& nucleus)
{
  // The same pointer of primary nucleus
  _nucleus = new G4Fragment(nucleus);
  _contDeexcitation->SetNucleus(_nucleus);
  _discrDeexcitation->SetNucleus(_nucleus);

  //G4cout << "G4PhotonEvaporation::BreakItUp:  " << nucleus << G4endl;

  // Do the whole gamma chain 
  G4FragmentVector* products = _contDeexcitation->DoChain();  
  if( !products ) { products = new G4FragmentVector; }

  // Products from continuum gamma transitions
  if (_verbose > 0) {
    G4cout << " = BreakItUp = " << products->size()
	   << " gammas from ContinuumDeexcitation " << G4endl;
  }

  // Products from discrete gamma transitions
  G4FragmentVector* discrProducts = _discrDeexcitation->DoChain();
  if(discrProducts) {
    _eOccupancy = _discrDeexcitation->GetEO();
    _vShellNumber = _discrDeexcitation->GetVacantSN();

    // not sure if the following line is needed!
    _discrDeexcitation->SetVaccantSN(-1);

    if (_verbose > 0) {
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
  products->push_back(_nucleus);

  if (_verbose > 0) {
    G4cout << "*-*-* Photon evaporation: " << products->size() << G4endl;
  }
  return products;
}

G4double G4PhotonEvaporation::GetEmissionProbability() const
{
  G4double prob = 
    _probAlgorithm->EmissionProbability(*_nucleus,_nucleus->GetExcitationEnergy());
  return prob;
}


void G4PhotonEvaporation::SetEmissionStrategy(G4VEmissionProbability * probAlgorithm)
{

  // CD - not sure about always wanting to delete this pointer....

  if(_myOwnProbAlgorithm) delete _probAlgorithm;

  _probAlgorithm = probAlgorithm;

  _myOwnProbAlgorithm = false;
}


void G4PhotonEvaporation::SetVerboseLevel(G4int verbose)
{
  _verbose = verbose;
  _contDeexcitation->SetVerboseLevel(verbose);
  _discrDeexcitation->SetVerboseLevel(verbose);
}

void G4PhotonEvaporation::SetICM(G4bool ic)
{
 (dynamic_cast <G4DiscreteGammaDeexcitation*> (_discrDeexcitation))->SetICM(ic);
}

void G4PhotonEvaporation::SetMaxHalfLife(G4double hl)
{
 (dynamic_cast <G4DiscreteGammaDeexcitation*> (_discrDeexcitation))->SetHL(hl);
}

void G4PhotonEvaporation::RDMForced(G4bool fromRDM)
{
 (dynamic_cast <G4DiscreteGammaDeexcitation*> (_discrDeexcitation))->SetRDM(fromRDM);
}

void G4PhotonEvaporation::SetEOccupancy(G4ElectronOccupancy eo)
{
  _discrDeexcitation->SetEO(eo);
}

#ifdef debug
void G4PhotonEvaporation::CheckConservation(const G4Fragment & theInitialState,
					    G4FragmentVector * Result) const
{
  G4double ProductsEnergy =0;
  G4ThreeVector ProductsMomentum;
  G4int ProductsA = 0;
  G4int ProductsZ = 0;
  G4FragmentVector::iterator h;
  for (h = Result->begin(); h != Result->end(); h++) {
    G4LorentzVector tmp = (*h)->GetMomentum();
    ProductsEnergy += tmp.e();
    ProductsMomentum += tmp.vect();
    ProductsA += static_cast<G4int>((*h)->GetA());
    ProductsZ += static_cast<G4int>((*h)->GetZ());
  }

  if (ProductsA != theInitialState.GetA()) {
    G4cout << "!!!!!!!!!! Baryonic Number Conservation Violation !!!!!!!!!!" << G4endl;
    G4cout << "G4PhotonEvaporation.cc: Barionic Number Conservation test for evaporation fragments" 
	   << G4endl; 
    G4cout << "Initial A = " << theInitialState.GetA() 
	   << "   Fragments A = " << ProductsA << "   Diference --> " 
	   << theInitialState.GetA() - ProductsA << G4endl;
  }
  if (ProductsZ != theInitialState.GetZ()) {
    G4cout << "!!!!!!!!!! Charge Conservation Violation !!!!!!!!!!" << G4endl;
    G4cout << "G4PhotonEvaporation.cc: Charge Conservation test for evaporation fragments" 
	   << G4endl; 
    G4cout << "Initial Z = " << theInitialState.GetZ() 
	   << "   Fragments Z = " << ProductsZ << "   Diference --> " 
	   << theInitialState.GetZ() - ProductsZ << G4endl;
  }
  if (std::abs(ProductsEnergy-theInitialState.GetMomentum().e()) > 1.0*keV) {
    G4cout << "!!!!!!!!!! Energy Conservation Violation !!!!!!!!!!" << G4endl;
    G4cout << "G4PhotonEvaporation.cc: Energy Conservation test for evaporation fragments" 
	   << G4endl; 
    G4cout << "Initial E = " << theInitialState.GetMomentum().e()/MeV << " MeV"
	   << "   Fragments E = " << ProductsEnergy/MeV  << " MeV   Diference --> " 
	   << (theInitialState.GetMomentum().e() - ProductsEnergy)/MeV << " MeV" << G4endl;
  } 
  if (std::abs(ProductsMomentum.x()-theInitialState.GetMomentum().x()) > 1.0*keV || 
      std::abs(ProductsMomentum.y()-theInitialState.GetMomentum().y()) > 1.0*keV ||
      std::abs(ProductsMomentum.z()-theInitialState.GetMomentum().z()) > 1.0*keV) {
    G4cout << "!!!!!!!!!! Momentum Conservation Violation !!!!!!!!!!" << G4endl;
    G4cout << "G4PhotonEvaporation.cc: Momentum Conservation test for evaporation fragments" 
	   << G4endl; 
    G4cout << "Initial P = " << theInitialState.GetMomentum().vect() << " MeV"
	   << "   Fragments P = " << ProductsMomentum  << " MeV   Diference --> " 
	   << theInitialState.GetMomentum().vect() - ProductsMomentum << " MeV" << G4endl;
  }
  return;
}
#endif



