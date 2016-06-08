
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//      
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
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


G4PhotonEvaporation::G4PhotonEvaporation()
{ 
  _probAlgorithm = new G4E1Probability;
  _myOwnProbAlgorithm = true;
  _discrDeexcitation = new G4DiscreteGammaDeexcitation;
  _contDeexcitation = new G4ContinuumGammaDeexcitation;
  _verbose = 0;
  _gammaE = 0.;
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
  if (contProducts != 0) nCont = contProducts->entries();

  G4int i;
  if (nCont > 0)
    {
      G4Fragment modifiedNucleus = _contDeexcitation->GetNucleus();
      _discrDeexcitation->SetNucleus(modifiedNucleus);
      for (i=0; i<nCont; i++)
	{
	  products->insert(contProducts->at(i));
	}
    }
  else
    {
      // Products from discrete gamma transitions
      G4FragmentVector* discrProducts = _discrDeexcitation->DoTransition();
      
      G4int nDiscr = 0;
      if (discrProducts != 0) nDiscr = discrProducts->entries();
      
      if (_verbose > 0)
	G4cout << " = BreakUp = " << nDiscr 
	       << " gammas from DiscreteDeexcitation " 
	       << endl;
      
      for (i=0; i<nDiscr; i++)
	{
	  products->insert(discrProducts->at(i));
	}
      delete discrProducts;
    }


  _gammaE = 0.;
  if (products->entries() > 0)
    {
      _gammaE = products->at(0)->GetMomentum().e();
    }

  delete contProducts;  // delete vector, not fragments 

  // Add deexcited nucleus to products
  G4Fragment* finalNucleus = new G4Fragment(_discrDeexcitation->GetNucleus());
  products->insert(finalNucleus);


  if (_verbose > 0)
    G4cout << "*-*-*-* Photon evaporation: " << products->entries() << endl;

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
  if (contProducts != 0) nCont = contProducts->entries();

  if (_verbose > 0)
    G4cout << " = BreakItUp = " << nCont 
	   << " gammas from ContinuumDeexcitation " << endl;

  G4int i;
  if (nCont > 0)
    {
      G4Fragment modifiedNucleus = _contDeexcitation->GetNucleus();
      _discrDeexcitation->SetNucleus(modifiedNucleus);
      for (i=0; i<nCont; i++)
	{
	  products->insert(contProducts->at(i));
	}
    }

  // Products from discrete gamma transitions
  G4FragmentVector* discrProducts = _discrDeexcitation->DoChain();

  G4int nDiscr = 0;
  if (discrProducts != 0) nDiscr = discrProducts->entries();

  if (_verbose > 0)
    G4cout << " = BreakItUp = " << nDiscr 
	   << " gammas from DiscreteDeexcitation " << endl;

  for (i=0; i<nDiscr; i++)
    {
      products->insert(discrProducts->at(i));
    }

  // Add deexcited nucleus to products
  G4Fragment* finalNucleus = new G4Fragment(_discrDeexcitation->GetNucleus());
  products->insert(finalNucleus);

  if (_verbose > 0)
    G4cout << " = BreakItUp = Nucleus added to products" << endl;

  delete contProducts;  // delete vector, not fragments 
  delete discrProducts;

  if (_verbose > 0)
    G4cout << "*-*-* Photon evaporation: " << products->entries() << endl;

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
