//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//

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
//      Modifications: 
//      
//        8 March 2002, Fan Lei (flei@space.qinetiq.com)
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
// -------------------------------------------------------------------

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
#include "G4AtomicDeexcitation.hh"

G4PhotonEvaporation::G4PhotonEvaporation()
  :_verbose(0),_myOwnProbAlgorithm (true),
   _eOccupancy(0), _vShellNumber(-1),_gammaE(0.),_applyARM(false)
{ 
  _probAlgorithm = new G4E1Probability;
  _discrDeexcitation = new G4DiscreteGammaDeexcitation;
  _contDeexcitation = new G4ContinuumGammaDeexcitation;
  (dynamic_cast <G4DiscreteGammaDeexcitation*> (_discrDeexcitation))->SetICM(false);
}


G4PhotonEvaporation::~G4PhotonEvaporation()
{ 
  if(_myOwnProbAlgorithm) delete _probAlgorithm;
  delete _discrDeexcitation;
  delete _contDeexcitation;
}


void G4PhotonEvaporation::Initialize(const G4Fragment& fragment)
{
  _nucleus = fragment;
  return;
}


G4FragmentVector* G4PhotonEvaporation::BreakUp(const G4Fragment& nucleus)
{
  _nucleus = nucleus;

  G4FragmentVector* products = new G4FragmentVector;
 
  _contDeexcitation->SetNucleus(nucleus);
  _discrDeexcitation->SetNucleus(nucleus);
  
  // Do one photon emission

  // Products from continuum gamma transitions

  G4FragmentVector* contProducts = _contDeexcitation->DoTransition();  

  G4int nCont = 0;
  if (contProducts != 0) nCont = contProducts->size();

  G4FragmentVector::iterator i;
  if (nCont > 0)
    {
      G4Fragment modifiedNucleus = _contDeexcitation->GetNucleus();
      _discrDeexcitation->SetNucleus(modifiedNucleus);
      for (i = contProducts->begin(); i != contProducts->end(); i++)
	{
	  products->push_back(*i);
	}
    }
  else
    {
      // Products from discrete gamma transitions
      G4FragmentVector* discrProducts = _discrDeexcitation->DoTransition();
      
      G4int nDiscr = 0;
      if (discrProducts != 0) nDiscr = discrProducts->size();
      
      if (_verbose > 0)
	G4cout << " = BreakUp = " << nDiscr 
	       << " gammas from DiscreteDeexcitation " 
	       << G4endl;
      
      for (i = discrProducts->begin(); i != discrProducts->end(); i++)
	{
	  products->push_back(*i);
	}
      discrProducts->clear();
      delete discrProducts;
    }


  _gammaE = 0.;
  if (products->size() > 0)
    {
      _gammaE = (*(products->begin()))->GetMomentum().e();
    }

  contProducts->clear();
  delete contProducts;  // delete vector, not fragments 

  // Add deexcited nucleus to products
  G4Fragment* finalNucleus = new G4Fragment(_discrDeexcitation->GetNucleus());
  products->push_back(finalNucleus);


  if (_verbose > 0)
    G4cout << "*-*-*-* Photon evaporation: " << products->size() << G4endl;

  return products;
}


G4FragmentVector* G4PhotonEvaporation::BreakItUp(const G4Fragment& nucleus)
{
  _nucleus = nucleus;

  G4FragmentVector* products = new G4FragmentVector;

  _contDeexcitation->SetNucleus(nucleus);
  _discrDeexcitation->SetNucleus(nucleus);

  // Do the whole gamma chain 

  G4FragmentVector* contProducts = _contDeexcitation->DoChain();  

  // Products from continuum gamma transitions
  G4int nCont = 0;
  if (contProducts != 0) nCont = contProducts->size();

  if (_verbose > 0)
    G4cout << " = BreakItUp = " << nCont 
	   << " gammas from ContinuumDeexcitation " << G4endl;

  G4FragmentVector::iterator i;
  if (nCont > 0)
    {
      G4Fragment modifiedNucleus = _contDeexcitation->GetNucleus();
      _discrDeexcitation->SetNucleus(modifiedNucleus);
      for (i = contProducts->begin(); i != contProducts->end(); i++)
	{
	  products->push_back(*i);
	}
    }

  // Products from discrete gamma transitions
  G4FragmentVector* discrProducts = _discrDeexcitation->DoChain();
  _eOccupancy = _discrDeexcitation->GetEO();
  _vShellNumber = _discrDeexcitation->GetVacantSN();
  // _eOccupancy.DumpInfo() ;
  _discrDeexcitation->SetVaccantSN(-1);

  G4int nDiscr = 0;
  if (discrProducts != 0) nDiscr = discrProducts->size();

  if (_verbose > 0)
    G4cout << " = BreakItUp = " << nDiscr 
	   << " gammas from DiscreteDeexcitation " << G4endl;

  for (i = discrProducts->begin(); i != discrProducts->end(); i++)
    {
      products->push_back(*i);
    }

  // now to see if apply Atomic relaxation model or not
  // only when there is a vacant orbital electron
  
  if (_applyARM && _vShellNumber != -1) 
    {
      G4int aZ = static_cast<G4int>(_discrDeexcitation->GetNucleus().GetZ());
      G4int eShell = _vShellNumber+1;
      if ( eShell > 0 ) {
	G4AtomicDeexcitation* atomDeex = new G4AtomicDeexcitation();
	// no auger electron for now
	atomDeex->ActivateAugerElectronProduction(0);
	std::vector<G4DynamicParticle*>* armProducts = atomDeex->GenerateParticles(aZ,eShell);
	G4DynamicParticle* aParticle;	
	if (_verbose > 0)
	  G4cout << " = BreakItUp = " << armProducts->size()
		 << " particles from G4AtomicDeexcitation " << G4endl;
	for (size_t i = 0;  i < armProducts->size(); i++)
	  {
	    aParticle = (*armProducts)[i] ;
	    G4LorentzVector lParticle = aParticle->Get4Momentum();
	    G4Fragment* aFragment = new
	      G4Fragment(lParticle,aParticle->GetDefinition());
	    aFragment->SetCreationTime(aParticle->GetProperTime());
	    products->push_back(aFragment);
	  }
	delete armProducts;
	delete atomDeex;
      }
    }
  // Add deexcited nucleus to products
  G4Fragment* finalNucleus = new G4Fragment(_discrDeexcitation->GetNucleus());
  products->push_back(finalNucleus);
  if (_verbose > 0)
    G4cout << " = BreakItUp = Nucleus added to products" << G4endl;

  if (_verbose > 0)
    G4cout << "*-*-* Photon evaporation: " << products->size() << G4endl;

#ifdef debug
  if ( armProducts->size() == 0)
    CheckConservation(nucleus,products);
#endif
  contProducts->clear();
  discrProducts->clear();
  delete contProducts;  // delete vector, not fragments 
  delete discrProducts;
  return products;
}


G4double G4PhotonEvaporation::GetEmissionProbability() const
{
  G4double prob = 0.;
  if (_probAlgorithm != 0) prob = _probAlgorithm->EmissionProbability(_nucleus,_gammaE);
  return prob;
}


void G4PhotonEvaporation::SetEmissionStrategy(G4VEmissionProbability * probAlgorithm)
{

  // CD - not sure about always wanting to delete this pointer....

  if(_myOwnProbAlgorithm) delete _probAlgorithm;

  _probAlgorithm = probAlgorithm;

  _myOwnProbAlgorithm = false;

  return;
}


void G4PhotonEvaporation::SetVerboseLevel(G4int verbose)
{
  _verbose = verbose;
  _contDeexcitation->SetVerboseLevel(verbose);
  _discrDeexcitation->SetVerboseLevel(verbose);

  return;
}

void G4PhotonEvaporation::SetICM(G4bool ic)
{
 (dynamic_cast <G4DiscreteGammaDeexcitation*> (_discrDeexcitation))->SetICM(ic);
 return;
}

void G4PhotonEvaporation::SetMaxHalfLife(G4double hl)
{
 (dynamic_cast <G4DiscreteGammaDeexcitation*> (_discrDeexcitation))->SetHL(hl);
 return;
}

void G4PhotonEvaporation::RDMForced(G4bool fromRDM)
{
 (dynamic_cast <G4DiscreteGammaDeexcitation*> (_discrDeexcitation))->SetRDM(fromRDM);
 return;
}

void G4PhotonEvaporation::SetEOccupancy(G4ElectronOccupancy eo)
{
  _discrDeexcitation->SetEO(eo);
  return;
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
  if (abs(ProductsEnergy-theInitialState.GetMomentum().e()) > 1.0*keV) {
    G4cout << "!!!!!!!!!! Energy Conservation Violation !!!!!!!!!!" << G4endl;
    G4cout << "G4PhotonEvaporation.cc: Energy Conservation test for evaporation fragments" 
	   << G4endl; 
    G4cout << "Initial E = " << theInitialState.GetMomentum().e()/MeV << " MeV"
	   << "   Fragments E = " << ProductsEnergy/MeV  << " MeV   Diference --> " 
	   << (theInitialState.GetMomentum().e() - ProductsEnergy)/MeV << " MeV" << G4endl;
  } 
  if (abs(ProductsMomentum.x()-theInitialState.GetMomentum().x()) > 1.0*keV || 
      abs(ProductsMomentum.y()-theInitialState.GetMomentum().y()) > 1.0*keV ||
      abs(ProductsMomentum.z()-theInitialState.GetMomentum().z()) > 1.0*keV) {
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



