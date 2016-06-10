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
// $Id: G4PromptPhotonEvaporation.cc 85841 2014-11-05 15:35:06Z gcosmo $
//
// -------------------------------------------------------------------
//
//      GEANT4 class file
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4PromptPhotonEvaporation
//
//      Author:        Vladimir Ivantchenko
//
//      Creation date: 20 December 2011
//
//Modifications:
//
// 
// -------------------------------------------------------------------
//

#include "G4PromptPhotonEvaporation.hh"

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Gamma.hh"
#include "G4NuclearLevelManager.hh"
#include "G4NuclearLevelStore.hh"

#include "G4LorentzVector.hh"
#include "G4VGammaTransition.hh"
#include "G4Fragment.hh"
#include "G4FragmentVector.hh"
#include "G4ContinuumGammaDeexcitation.hh"
#include "G4DiscreteGammaDeexcitation.hh"
#include "G4E1Probability.hh"

G4PromptPhotonEvaporation::G4PromptPhotonEvaporation()
  :fVerbose(0), fICM(true), fRDM(false), fMaxHalfTime(DBL_MAX),
   fEmissionProbability(0.0),levelManager(0),nucleus(0)
{
  fNuclearLevelStore = G4NuclearLevelStore::GetInstance(); 
  theA = theZ = 0;
  fEnergyFermi = fExcEnergyMax = gammaE = 0.0;
}

G4PromptPhotonEvaporation::~G4PromptPhotonEvaporation()
{ 
}

G4double 
G4PromptPhotonEvaporation::GetEmissionProbability(G4Fragment* theNucleus) 
{
  fEmissionProbability = 0.0;
  nucleus = theNucleus;
  G4double ex = nucleus->GetExcitationEnergy();

  if(nucleus->GetZ_asInt() != theZ || nucleus->GetA_asInt() != theA) {
    G4int Z = nucleus->GetZ_asInt();
    G4int A = nucleus->GetA_asInt();
    fExcEnergyMax = -1.0;
    if(1 < A && ex > keV) {
      fEnergyFermi = G4NucleiProperties::GetNuclearMass(A-1, Z) 
	+ neutron_mass_c2 - nucleus->GetGroundStateMass();
      fExcEnergyMax = fEnergyFermi + 15*MeV;
    }
    if(ex < fExcEnergyMax) {

      theZ = Z;
      theA = A;
      levelManager = fNuclearLevelStore->GetManager(Z,A);

      // continium transition
      if(ex >= fEnergyFermi) {

	// discrete transition
      } else {
      }
    }
  }
  return fEmissionProbability;
}

G4Fragment* 
G4PromptPhotonEvaporation::EmittedFragment(G4Fragment* theNucleus)
{
  //G4cout << "G4PromptPhotonEvaporation::EmittedFragment" << G4endl;
  
  G4Fragment* gamma = 0;
  if(theNucleus->GetExcitationEnergy() <= keV) { return gamma; }
  if(GetEmissionProbability(theNucleus) <= 0.0){ return gamma; }

  //G4cout << "G4PromptPhotonEvaporation::EmittedFragment" << G4endl;
  /*
    G4Fragment* gamma = _contDeexcitation->GenerateGamma();
    if(gamma) { 
      if (_verbose > 0) {
	G4cout << "G4PromptPhotonEvaporation::EmittedFragment continium deex: "   
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
	G4cout << "G4PromptPhotonEvaporation::EmittedFragment discrete deex: "   
	       << gamma << G4endl;
        G4cout << "   Residual: " << nucleus << G4endl;
      }
      return gamma; 
    }
  }
 
  if (_verbose > 0) {
    G4cout << "G4PromptPhotonEvaporation unable emit gamma: " 
	   << nucleus << G4endl;
  }
  */
  return gamma;
}

G4bool 
G4PromptPhotonEvaporation::BreakUpChain(G4FragmentVector*, 
					G4Fragment*)
{
  return false;
}

G4FragmentVector* 
G4PromptPhotonEvaporation::BreakUpFragment(G4Fragment* theNucleus)
{
  //G4cout << "G4PromptPhotonEvaporation::BreakUpFragment" << G4endl;
  G4FragmentVector* products = new G4FragmentVector();
  BreakUpChain(products, theNucleus);
  return products;
}

G4FragmentVector* 
G4PromptPhotonEvaporation::BreakUp(const G4Fragment& theNucleus)
{
  //G4cout << "G4PromptPhotonEvaporation::BreakUp" << G4endl;
  G4Fragment* aNucleus = new G4Fragment(theNucleus);
  G4FragmentVector* products = new G4FragmentVector();
  BreakUpChain(products, aNucleus);
  products->push_back(aNucleus);
  return products;
}


