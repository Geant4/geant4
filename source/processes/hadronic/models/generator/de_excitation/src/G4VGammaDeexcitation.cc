// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
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
//      File name:     G4VGammaDeexcitation
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 23 October 1998
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#include "G4VGammaDeexcitation.hh"

#include "globals.hh"
#include "Randomize.hh"
#include "G4Gamma.hh"
#include "G4LorentzVector.hh"
#include "G4VGammaTransition.hh"
#include "G4Fragment.hh"
#include "G4FragmentVector.hh"


G4VGammaDeexcitation::G4VGammaDeexcitation(): _verbose(0), _transition(0)
{ }


G4VGammaDeexcitation::~G4VGammaDeexcitation()
{ 
  //  if (_transition != 0) delete _transition;
}


G4FragmentVector* G4VGammaDeexcitation::DoTransition()
{
  // Template method

 Initialize();
 G4FragmentVector* products = new G4FragmentVector;
  
  if (CanDoTransition())
    {
     G4Fragment* gamma = GenerateGamma();
     if (gamma != 0)
       {
	 products->append(gamma);
	 UpdateNucleus(gamma);
	 Update(gamma);
       }
    }

  if (_verbose > 1)
    G4cout << "G4VGammaDeexcitation::DoTransition - Transition deleted " << endl;

  if (_transition != 0) delete _transition;

  return products;
}

G4FragmentVector* G4VGammaDeexcitation::DoChain()
{
  Initialize();
  G4FragmentVector* products = new G4FragmentVector;

  while (CanDoTransition())
    {
      if (_verbose > 5) G4cout << "G4VGammaDeexcitation::DoChain -  Looping" << endl;

      G4Fragment* gamma = GenerateGamma();
      if (gamma != 0) 
	{
	  products->append(gamma);
	  UpdateNucleus(gamma);
	}
     Update(gamma);
    } 

  if (_verbose > 1)
      G4cout << "G4VGammaDeexcitation::DoChain - Transition deleted, end of chain " << endl;

  if (_transition != 0) delete _transition;
  
  return products;
}


const G4Fragment& G4VGammaDeexcitation::GetNucleus() const
{
  return _nucleus; 
}


void G4VGammaDeexcitation::SetNucleus(const G4Fragment& nucleus)
{
  _nucleus = G4Fragment(nucleus);
}


G4Fragment* G4VGammaDeexcitation::GenerateGamma()
{
  G4double eGamma = 0.;

  if (_transition != 0) eGamma = _transition->GammaEnergy();

  if (_verbose > 1 && _transition != 0) 
    {
      G4cout << "G4VGammaDeexcitation::GenerateGamma - Gamma energy " << eGamma 
	     << " ** New excitation " << _transition->GetEnergyTo() << endl;
    }

  // Photon momentum isotropically generated 
 
    if (eGamma > 0.) 
      {
	G4double cosTheta = 1. - 2. * G4UniformRand();
	G4double sinTheta = sqrt(1. - cosTheta * cosTheta);
	G4double phi = twopi * G4UniformRand();
	
	G4ThreeVector pGamma( eGamma * sinTheta * cos(phi),
			      eGamma * sinTheta * sin(phi),
			      eGamma * cosTheta );
	
	G4LorentzVector gamma(pGamma, eGamma);
	//	gamma.boost(_nucleus.GetMomentum().boostVector() );
	G4Fragment* gammaFragment = new G4Fragment(gamma,G4Gamma::GammaDefinition());

        if (_verbose > 1)
	  G4cout << "G4VGammaDeexcitation::GenerateGamma -  Gamma fragment generated " << endl;

	return gammaFragment;
      }
    else
      {
	return 0;
      }
}


void G4VGammaDeexcitation::UpdateNucleus(const G4Fragment*  gamma)
{
  G4LorentzVector p4Gamma = gamma->GetMomentum();
  G4ThreeVector pGamma(p4Gamma);
  G4double eGamma = gamma->GetMomentum().e();
  
  G4LorentzVector p4Nucleus(_nucleus.GetMomentum() );
  //  p4Nucleus.boost(-_nucleus.GetMomentum().boostVector() );
  
  G4LorentzVector p4Residual(p4Nucleus - pGamma, p4Nucleus.e() - eGamma);
  //  G4LorentzVector p4Residual(-pGamma, p4Nucleus.e() - eGamma);
  //  p4Residual.boost( _nucleus.GetMomentum().boostVector() );
  
  // Update excited nucleus parameters

  _nucleus.SetMomentum(p4Residual);

  if (_transition != 0) 
  {
    G4double excitation =_transition->GetEnergyTo();
    if (excitation < 0.) excitation = 0.0;
    _nucleus.SetExcitationEnergy(excitation);
  }

  return;
}


void G4VGammaDeexcitation::Update(const G4Fragment*  gamma)
{
  if (_transition !=  0) 
    { 
      delete _transition;
      _transition = 0;
      if (_verbose > 1)
	G4cout << "G4VGammaDeexcitation::Update - Transition deleted " << endl;
    }

  _transition = CreateTransition();
  if (_transition != 0) 
    {
      _transition->SetEnergyFrom(_nucleus.GetExcitationEnergy());
    }

  return;
}


void G4VGammaDeexcitation::Initialize()
{
  _transition = CreateTransition();
  if (_transition != 0) _transition->SetEnergyFrom(_nucleus.GetExcitationEnergy());
  return;
}


void G4VGammaDeexcitation::SetVerboseLevel(G4int verbose)
{
  _verbose = verbose;
  return;
}



