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
// $Id: G4VFermiBreakUp.hh 85443 2014-10-29 14:35:57Z gcosmo $
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

  G4VFermiBreakUp();
  virtual ~G4VFermiBreakUp();

  // primary fragment is copied to the new instance, the copy is deleted 
  // or is added to the list of products 
  virtual G4FragmentVector * BreakItUp(const G4Fragment &theNucleus) = 0;

  // new interface - vector of products is added to the provided vector
  // primary fragment is deleted or is modified and added to the list
  // of products
  virtual void BreakFragment(G4FragmentVector*, G4Fragment* theNucleus) = 0;

private:

  G4VFermiBreakUp(const G4VFermiBreakUp &right);  
  const G4VFermiBreakUp & operator=(const G4VFermiBreakUp &right);
  G4bool operator==(const G4VFermiBreakUp &right) const;
  G4bool operator!=(const G4VFermiBreakUp &right) const;
};

#endif


