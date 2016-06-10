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
// $Id: G4ContinuumGammaDeexcitation.cc 88987 2015-03-17 10:39:50Z gcosmo $
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4ContinuumGammaDeexcitation
//
//      Authors:       Carlo Dallapiccola (dallapiccola@umdhep.umd.edu)
//                     Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 23 October 1998
//
//      Modifications: 
//
//      02 May 2003,   Vladimir Ivanchenko change interface to G4NuclearlevelManager
//
//      19 April 2010 J. M. Quesada: smaller value of tolerance parameter
//      
// -------------------------------------------------------------------
//
// Class G4ContinuumGammaDeexcitation.cc
//
//       Concrete class derived from G4VGammaDeexcitation
//
//
#include "G4ContinuumGammaDeexcitation.hh"

#include "G4Gamma.hh"
#include "G4ContinuumGammaTransition.hh"
#include "G4NuclearLevelManager.hh"
#include "G4NuclearLevelStore.hh"
#include "G4Fragment.hh"      
#include "G4ConstantLevelDensityParameter.hh"

G4ContinuumGammaDeexcitation::G4ContinuumGammaDeexcitation()
  : nucleusZ(0), nucleusA(0), levelManager(0), ctransition(0)
{
  store = G4NuclearLevelStore::GetInstance();
}

G4ContinuumGammaDeexcitation::~G4ContinuumGammaDeexcitation() 
{}

G4bool G4ContinuumGammaDeexcitation::CanDoTransition(G4Fragment* nucleus) 
{
  G4int Z = nucleus->GetZ_asInt();
  G4int A = nucleus->GetA_asInt();
  G4double excitation = nucleus->GetExcitationEnergy();

  if (excitation <= _tolerance) {
    if (_verbose > 1) { 
      G4cout << "G4ContinuumGammaDeexcitation::CanDoTransition fails; Z= " << Z
	     << " A= " << A << " Eex(meV)= " << excitation << G4endl;
    }
    return false;
  }

  if (nucleusA != A || nucleusZ != Z) {
    levelManager = store->GetManager(Z,A);
    nucleusA = A;
    nucleusZ = Z;
  }

  if (_verbose > 1) {
    G4cout << "G4ContinuumGammaDeexcitation: "
	   << " Z= " << Z << "  A= " << A << " Eex= " << excitation
	   << G4endl;
  }
  if (levelManager && excitation <= levelManager->MaxLevelEnergy() + _tolerance) {
    if (_verbose > 0) {
      G4cout << "G4ContinuumGammaDeexcitation::CanDoTransition -  Excitation " 
	     << excitation << " below max discrete level " 
	     << levelManager->MaxLevelEnergy() << G4endl;
    }
    return false;
  }

  if(!ctransition) {
    ctransition = 
      new G4ContinuumGammaTransition(levelManager,Z,A,excitation,_verbose);
    _transition = ctransition;
  } else {
    ctransition->Update(levelManager,Z,A,excitation);
  }
  /*
    G4cout <<"G4ContinuumGammaDeexcitation::CanDoTransition: " 
	   << " Eex(MeV)= " << excitation 
	   << " Emax(MeV)= " << levelManager->MaxLevelEnergy()
	   << " Z= " << Z << " A= " << A << G4endl;
  */
  return true;
}



