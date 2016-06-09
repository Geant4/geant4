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
//
// $Id$
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)
//
// Modifications:
// 01.04.2011 General cleanup by V.Ivanchenko

#ifndef G4Be8FermiFragment_h
#define G4Be8FermiFragment_h 1

#include "G4UnstableFermiFragment.hh"

class G4Be8FermiFragment : public G4UnstableFermiFragment
{
public:

  G4Be8FermiFragment(G4int anA, G4int aZ, G4int Pol, G4double ExE);

  virtual ~G4Be8FermiFragment();
  
private:
  
  G4Be8FermiFragment(const G4Be8FermiFragment &right);
  const G4Be8FermiFragment & operator=(const G4Be8FermiFragment &right);
  G4bool operator==(const G4Be8FermiFragment &right) const;
  G4bool operator!=(const G4Be8FermiFragment &right) const;
  
};


#endif


