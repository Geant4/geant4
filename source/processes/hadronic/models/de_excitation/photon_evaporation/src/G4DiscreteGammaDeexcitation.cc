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
// $Id: G4DiscreteGammaDeexcitation.cc 88987 2015-03-17 10:39:50Z gcosmo $
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4DiscreteGammaDeexcitation
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 23 October 1998
//
//      Modifications: 
//      8 March 2001, Fan Lei (flei@space.qinetiq.com)
//         Added the following as part if the IC implementation
//            void SetICM(G4bool hl) { _icm = hl; };
//            void SetRDM(G4bool hl) { _rdm = hl; };
//            void SetHL(G4double hl) { _max_hl = hl; };
//         Changed in CreateTransition() from 
//                     return new G4DiscreteGammaTransition(*level);  
//                 To
//                     return new G4DiscreteGammaTransition(*level,Z);
//         Added in CanDoTransition
//                if (level->HalfLife() > _max_hl && !_rdm ) canDo = false;
//
//		3 November 2011 L. Desorgher
//       		remove the   Z<= 98  limit
//      
// -------------------------------------------------------------------

#include "G4DiscreteGammaDeexcitation.hh"
#include "G4DiscreteGammaTransition.hh"
#include "G4NuclearLevelManager.hh"
#include "G4NuclearLevelStore.hh"
#include "G4SystemOfUnits.hh"

G4DiscreteGammaDeexcitation::G4DiscreteGammaDeexcitation(): 
  nucleusZ(0), nucleusA(0), maxhl(1e-6*CLHEP::second), icm(false),
  rdm(false), levelManager(0), dtransition(0)
{
  store = G4NuclearLevelStore::GetInstance();
}

G4DiscreteGammaDeexcitation::~G4DiscreteGammaDeexcitation() 
{}

G4bool G4DiscreteGammaDeexcitation::CanDoTransition(G4Fragment* nucleus) 
{
  G4int Z = nucleus->GetZ_asInt();
  G4int A = nucleus->GetA_asInt();
  G4double excitation = nucleus->GetExcitationEnergy();

  if(excitation <= _tolerance) {
    if (_verbose > 1) { 
      G4cout << "G4DiscreteGammaDeexcitation::CanDoTransition fails; Z= " << Z
	     << " A= " << A << " Eex(meV)= " << excitation/MeV << G4endl;
    }
    return false;
  }

  if (nucleusA != A || nucleusZ != Z) {
    levelManager = store->GetManager(Z,A);
    nucleusA = A;
    nucleusZ = Z;
  }
  if(!levelManager ||
     excitation > levelManager->MaxLevelEnergy() + _tolerance) { return false; }

  if (_verbose > 1) {
    G4cout << "G4DiscreteGammaDeexcitation::CanDoTransition "
	   << " Z= " << Z << "  A= " << A << " Eex= " << excitation
	   << G4endl;
  }

  const G4NuclearLevel* level = levelManager->NearestLevel(excitation);

  // long lived level	
  if (!level || level->HalfLife() > maxhl) { return false; }
  //  if (!level || (level->HalfLife() > maxhl && !rdm) ) { return false; }

  if (_verbose > 1) {
    G4cout << "G4DiscreteGammaDeexcitation: Elevel(MeV)= " 
	   << level->Energy()/MeV << ", Eex(MeV)= " << excitation << G4endl;
  }
  if(!_transition) {
    dtransition = new G4DiscreteGammaTransition(level,Z,_verbose);
    dtransition->SetICM(icm);  
    _transition = dtransition;
  } else {
    dtransition->Update(level,Z);
  }
  // control on ICM
  if(level->HalfLife() > _timeLimit) {
    dtransition->SetICM(true);
  } else {
    dtransition->SetICM(icm);  
  }  

  dtransition->SetEnergyFrom(excitation);
  
  if (_verbose > 1) {
    G4cout << "CanDoTransition done: Eex(MeV)= " 
	   << excitation/MeV << ", level enrgies: Emin(MeV)= " 
	   << levelManager->MinLevelEnergy()/MeV << " Emax(MeV)= "
	   << levelManager->MaxLevelEnergy()/MeV << G4endl;
  }
  return true;
}

