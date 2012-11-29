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
// $Id$
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
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
//        21 Nov 2001, Fan Lei (flei@space.qinetiq.com)
//           Modified GenerateGamma() and UpdateUncleus() for implementation
//           of Internal Conversion processs 
//      
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added creation time evaluation for products of evaporation
//
//        19 April 2010, J. M. Quesada calculations in CM system
//              pending final boost to lab system  (not critical)
//
//        23 April 2010, V.Ivanchenko rewite kinematic part using PDG formula
//                                    for 2-body decay
//
//        07 May   2011, V.Ivanchenko implement check ICM flag - produce or not e-
//
// -------------------------------------------------------------------

#include "G4VGammaDeexcitation.hh"

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4LorentzVector.hh"
#include "G4VGammaTransition.hh"
#include "G4Fragment.hh"
#include "G4FragmentVector.hh"

#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

#include "G4DiscreteGammaTransition.hh"


G4VGammaDeexcitation::G4VGammaDeexcitation(): _transition(0), _verbose(0),
					      _electronO (0), _vSN(-1)
{ 
  _nucleus = 0;
  fTimeLimit = DBL_MAX;
}

G4VGammaDeexcitation::~G4VGammaDeexcitation()
{ 
  if (_transition != 0) { delete _transition; }
}

G4FragmentVector* G4VGammaDeexcitation::DoTransition()
{
  Initialize();
  G4FragmentVector* products = new G4FragmentVector();
 
  if (CanDoTransition())
    {
      G4Fragment* gamma = GenerateGamma();
      if (gamma != 0) { products->push_back(gamma); }
    }
 
  if (_verbose > 1) {
    G4cout << "G4VGammaDeexcitation::DoTransition - Transition deleted " << G4endl;
  }
 
  return products;
}

G4FragmentVector* G4VGammaDeexcitation::DoChain()
{
  if (_verbose > 1) { G4cout << "G4VGammaDeexcitation::DoChain" << G4endl; }
  const G4double tolerance = CLHEP::keV;

  Initialize();
  G4FragmentVector* products = new G4FragmentVector();
  
  while (CanDoTransition())
    {      
      _transition->SetEnergyFrom(_nucleus->GetExcitationEnergy());
      G4Fragment* gamma = GenerateGamma();
      if (gamma != 0) 
	{
	  products->push_back(gamma);
	  //G4cout << "Eex(keV)= " << _nucleus->GetExcitationEnergy()/keV << G4endl;
	  if(_nucleus->GetExcitationEnergy() <= tolerance) { break; }
	  Update();
	}
    } 
  
  if (_verbose > 1) {
    G4cout << "G4VGammaDeexcitation::DoChain - Transition deleted, end of chain " << G4endl;
  }
  
  return products;
}

G4Fragment* G4VGammaDeexcitation::GenerateGamma()
{
  // 23/04/10 V.Ivanchenko rewrite complitely
  G4double eGamma = 0.;
  
  if (_transition) {
    _transition->SelectGamma();  // it can be conversion electron too
    eGamma = _transition->GetGammaEnergy(); 
    //G4cout << "G4VGammaDeexcitation::GenerateGamma - Egam(MeV)= " 
    //	   << eGamma << G4endl; 
    if(eGamma <= 0.0) { return 0; }
  } else { return 0; }

  G4double excitation = _nucleus->GetExcitationEnergy() - eGamma;
  if(excitation < 0.0) { excitation = 0.0; } 
  if (_verbose > 1) 
    {
      G4cout << "G4VGammaDeexcitation::GenerateGamma - Edeexc(MeV)= " << eGamma 
	     << " ** left Eexc(MeV)= " << excitation
	     << G4endl;
    }

  G4double gammaTime = _transition->GetGammaCreationTime();
  
  // Do complete Lorentz computation 
  G4LorentzVector lv = _nucleus->GetMomentum();
  G4double Mass = _nucleus->GetGroundStateMass() + excitation;

  // select secondary
  G4ParticleDefinition* gamma = G4Gamma::Gamma();

  G4DiscreteGammaTransition* dtransition = 
    dynamic_cast <G4DiscreteGammaTransition*> (_transition);

  G4bool eTransition = false;
  if (dtransition && !( dtransition->IsAGamma()) ) {
    eTransition = true; 
    gamma = G4Electron::Electron(); 
    _vSN = dtransition->GetOrbitNumber();   
    _electronO.RemoveElectron(_vSN);
    lv += G4LorentzVector(0.0,0.0,0.0,CLHEP::electron_mass_c2 - dtransition->GetBondEnergy());
  }

  G4double cosTheta = 1. - 2. * G4UniformRand(); 
  G4double sinTheta = std::sqrt(1. - cosTheta * cosTheta);
  G4double phi = twopi * G4UniformRand();

  G4double eMass = gamma->GetPDGMass();
  G4LorentzVector Gamma4P;

  //G4cout << "Egamma= " << eGamma << " Mass= " << eMass << " t= " << gammaTime
  //	 << " tlim= " << fTimeLimit << G4endl;

  if(gammaTime > fTimeLimit) {
    // shortcut for long lived levels
    // not correct position of stopping ion gamma emission
    // 4-momentum balance is breaked
    G4double e = eGamma + eMass;
    G4double mom = std::sqrt(eGamma*(eGamma + 2*eMass));
    Gamma4P.set(mom * sinTheta * std::cos(phi),
		mom * sinTheta * std::sin(phi),
		mom * cosTheta, e); 
    lv -= Gamma4P;
    e = lv.e();
    if(e < Mass) { e = Mass; }
    mom = std::sqrt((e - Mass)*(e + Mass));
    G4ThreeVector v = lv.vect().unit();
    lv.set(mom*v.x(), mom*v.y(), mom*v.z(), e); 

  } else {
    // 2-body decay in rest frame
    G4double Ecm       = lv.mag();
    G4ThreeVector bst  = lv.boostVector();

    G4double GammaEnergy = 0.5*((Ecm - Mass)*(Ecm + Mass) + eMass*eMass)/Ecm;
    if(GammaEnergy <= eMass) { return 0; }

    G4double mom = std::sqrt((GammaEnergy - eMass)*(GammaEnergy + eMass));
    Gamma4P.set(mom * sinTheta * std::cos(phi),
		mom * sinTheta * std::sin(phi),
		mom * cosTheta,
		GammaEnergy);

    Gamma4P.boost(bst);  
    lv -= Gamma4P;
  }

  // modified primary fragment 
  gammaTime += _nucleus->GetCreationTime();

  _nucleus->SetMomentum(lv);
  _nucleus->SetCreationTime(gammaTime);

  // e- is not produced
  if(eTransition && !dtransition->GetICM()) { return 0; }

  // gamma or e- are produced
  G4Fragment * thePhoton = new G4Fragment(Gamma4P,gamma);
  thePhoton->SetCreationTime(gammaTime);

  //G4cout << "G4VGammaDeexcitation::GenerateGamma : " << thePhoton << G4endl;
  //G4cout << "       Left nucleus: " << _nucleus << G4endl;
  return thePhoton;
}

void G4VGammaDeexcitation::Update()
{
  if (_transition !=  0) 
    { 
      delete _transition;
      _transition = 0;
      if (_verbose > 1) {
	G4cout << "G4VGammaDeexcitation::Update - Transition deleted " << G4endl;
      }
    }
  
  _transition = CreateTransition();
  if (_transition != 0) 
    {
      _transition->SetEnergyFrom(_nucleus->GetExcitationEnergy());
      // if ( _vSN != -1) (dynamic_cast <G4DiscreteGammaTransition*> (_transition))->SetICM(false);
      // the above line is commented out for bug fix #952. It was intruduced for reason that
      // the k-shell electron is most likely one to be kicked out and there is no time for 
      // the atom to deexcite before the next IC. But this limitation is causing other problems as 
      // reported in #952
    }
  
  return;
}
