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
// $Id: G4VFermiBreakUp.hh 98577 2016-07-25 13:05:12Z vnivanch $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)
//
// Modifications:
// 01.04.2011 General cleanup by V.Ivanchenko

#ifndef G4VFermiBreakUp_h
#define G4VFermiBreakUp_h 1

#include "globals.hh"
#include "G4FragmentVector.hh"

class G4VFermiBreakUp 
{
public:

  explicit G4VFermiBreakUp() {};
  virtual ~G4VFermiBreakUp() {};

  virtual void Initialise() = 0;

  // check if the Fermi Break Up model can be used 
  // mass is an effective mass of a fragment
  virtual G4bool IsApplicable(G4int Z, G4int A, G4double mass) const = 0;

  // vector of products is added to the provided vector
  // if no decay channel is found out for the primary fragment 
  // then it is added to the results vector
  // if primary decays then it is deleted 
  virtual void BreakFragment(G4FragmentVector* results, 
			     G4Fragment* theNucleus) = 0;

private:

  G4VFermiBreakUp(const G4VFermiBreakUp &right) = delete;  
  const G4VFermiBreakUp & operator=(const G4VFermiBreakUp &right) = delete;
  G4bool operator==(const G4VFermiBreakUp &right) const = delete;
  G4bool operator!=(const G4VFermiBreakUp &right) const = delete;
};

#endif


