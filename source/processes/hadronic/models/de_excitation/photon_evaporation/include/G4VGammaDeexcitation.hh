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
// $Id: G4VGammaDeexcitation.hh,v 1.9 2010-11-17 16:50:53 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added creation time evaluation for products of evaporation
//
//        21 Nov 2001, Fan Lei (flei@space.qinetiq.com)
//           Modified GenerateGamma() and UpdateUncleus() for implementation
//           of Internal Conversion processs
// 
//        8 March 2002, Fan Lei (flei@space.qinetiq.com)
//          Added  SetEO () , GetEO(), UpdateElectrons() to allow the assignment 
//          and modification of electron configuration.
//
//        18 October 2002, F. Lei
//          Added GetVaccantSN() and _vSN in order to link to ARM in low-e em
//          _vSN is updated in UpdateElectron()
//          Added SetVaccantSN(). It is need to to re-set _vSN after each 
//          IC happened.
//
//        28 April 2010, V.Ivanchenko cleanup methods
//
// -------------------------------------------------------------------
//

#ifndef G4VGAMMADEEXCITATION_HH
#define G4VGAMMADEEXCITATION_HH 1

#include "globals.hh"
#include "G4VGammaTransition.hh"
#include "G4Fragment.hh"
#include "G4FragmentVector.hh"
#include "G4ElectronOccupancy.hh"

class G4VGammaDeexcitation {

public:

  G4VGammaDeexcitation();

  virtual ~G4VGammaDeexcitation();

  virtual G4VGammaTransition * CreateTransition() = 0;
  virtual G4bool CanDoTransition() = 0;

  // Single gamma transition
  G4FragmentVector * DoTransition();

  // Chain of gamma transitions
  G4FragmentVector * DoChain();

  G4Fragment * GenerateGamma();

  inline G4Fragment* GetNucleus();

  inline void SetNucleus(G4Fragment* nucleus);

  inline void SetVerboseLevel(G4int verbose);

  inline void Initialize();

  inline void SetEO(G4ElectronOccupancy eo) { _electronO = eo; };
  inline void SetVaccantSN( G4int val ) { _vSN = val;};
  
  inline G4ElectronOccupancy GetEO() { return _electronO; };    
  inline G4int GetVacantSN() {return _vSN;};

protected:

  void Update();

  G4VGammaTransition* _transition; // Owned pointer
  G4int _verbose;

private:

  G4Fragment* _nucleus;
  G4ElectronOccupancy _electronO;
  G4int _vSN;

  G4VGammaDeexcitation(const G4VGammaDeexcitation & right);
  const G4VGammaDeexcitation & operator = (const G4VGammaDeexcitation & right);
  G4bool operator == (const G4VGammaDeexcitation & right) const;
  G4bool operator != (const G4VGammaDeexcitation & right) const;

};

inline G4Fragment* G4VGammaDeexcitation::GetNucleus() 
{
  return _nucleus; 
}

inline void G4VGammaDeexcitation::SetNucleus(G4Fragment* nucleus)
{
  _nucleus = nucleus;
}

inline void G4VGammaDeexcitation::SetVerboseLevel(G4int verbose)
{
  _verbose = verbose;
}

inline void G4VGammaDeexcitation::Initialize()
{
  if (_transition != 0) { delete _transition; }
  _transition = CreateTransition();
  if (_transition != 0) {
    _transition->SetEnergyFrom(_nucleus->GetExcitationEnergy());
  }
}

#endif






