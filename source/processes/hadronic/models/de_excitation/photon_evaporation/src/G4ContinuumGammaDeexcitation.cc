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
// $Id: G4ContinuumGammaDeexcitation.cc,v 1.8 2010-11-17 19:17:17 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

//
// Constructor
//

G4ContinuumGammaDeexcitation::G4ContinuumGammaDeexcitation()
  : _nucleusZ(0), _nucleusA(0), _levelManager(0)
{}

G4ContinuumGammaDeexcitation::~G4ContinuumGammaDeexcitation() 
{}

G4VGammaTransition* G4ContinuumGammaDeexcitation::CreateTransition()
{
  G4Fragment* nucleus = GetNucleus();
  G4int Z = nucleus->GetZ_asInt();
  G4int A = nucleus->GetA_asInt();
  G4double excitation = nucleus->GetExcitationEnergy();

  if (_nucleusA != A || _nucleusZ != Z)
    {
      _levelManager = G4NuclearLevelStore::GetInstance()->GetManager(Z,A);
      _nucleusA = A;
      _nucleusZ = Z;
    }

  if (_verbose > 1) {
    G4cout << "G4ContinuumGammaDeexcitation::CreateTransition - Created" << G4endl;
  }
  G4VGammaTransition* gt =  
    new G4ContinuumGammaTransition(_levelManager,Z,A,excitation,_verbose );

  return gt;
}
    

G4bool G4ContinuumGammaDeexcitation::CanDoTransition() 
{
  //JMQ: far too small, creating sometimes continuum gammas instead 
  //     of the right discrete ones (when excitation energy is slightly 
  //     over maximum discrete  energy): changed
  //  G4double tolerance = 10*eV;
  const G4double tolerance = CLHEP::keV;

  if (_transition == 0) 
    {
      if (_verbose > 0) {
	G4cout << "G4ContinuumGammaDeexcitation::CanDoTransition - Null transition "
	       << G4endl;
      }
      return false;
    }

  G4Fragment* nucleus = GetNucleus();
  G4double excitation = nucleus->GetExcitationEnergy();

  if (_nucleusZ < 2 || _nucleusA < 3)
    {
      if (_verbose > 1) { 
	G4cout << "G4ContinuumGammaDeexcitation::CanDoTransition - n/p/H"
	       << G4endl;
      }
      return false;
    }

  if (excitation <= tolerance) 
    {
      if (_verbose > 1) { 
	G4cout << "G4ContinuumGammaDeexcitation::CanDoTransition -  Excitation "
	       << excitation/CLHEP::keV << " keV is too small"
	       << G4endl;
      }
      return false;
    }
  if (excitation <= (_levelManager->MaxLevelEnergy() + tolerance)) 
    {  
      if (_verbose > 0) {
	G4cout << "G4ContinuumGammaDeexcitation::CanDoTransition -  Excitation " 
	       << excitation << " below max discrete level " 
	       << _levelManager->MaxLevelEnergy() << G4endl;
      }
      return false;
    }
  
  if (_verbose > 1) {
    G4cout <<"G4ContinuumGammaDeexcitation::CanDoTransition - CanDo" 
	   << " Eex(keV)= " << excitation/CLHEP::keV 
	   << " Emax(keV)= " << _levelManager->MaxLevelEnergy()/CLHEP::keV 
	   << " Z= " << _nucleusZ << " A= " << _nucleusA
	   << G4endl;
  }
  return true;
}



