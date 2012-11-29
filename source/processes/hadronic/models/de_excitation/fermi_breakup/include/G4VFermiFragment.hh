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
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)
//
// Modifications:
// 01.04.2011 General cleanup by V.Ivanchenko

#ifndef G4VFermiFragment_h
#define G4VFermiFragment_h 1

#include "G4FragmentVector.hh"
#include "globals.hh"

class G4VFermiFragment 
{
public:

  G4VFermiFragment(G4int anA, G4int aZ, G4int Pol, G4double ExE);

  virtual ~G4VFermiFragment();
  
private:

  G4VFermiFragment(const G4VFermiFragment &right);  
  const G4VFermiFragment & operator=(const G4VFermiFragment &right);
  G4bool operator==(const G4VFermiFragment &right) const;
  G4bool operator!=(const G4VFermiFragment &right) const;
  
public:

  virtual G4FragmentVector * GetFragment(const G4LorentzVector & aMomentum) const = 0;

  inline G4int GetA(void) const 
  {
    return A;
  }
  
  inline G4int GetZ(void) const 
  {
    return Z;
  }
  
  inline G4int GetPolarization(void) const 
  {
    return Polarization;
  }

  inline G4double GetExcitationEnergy(void) const 
  {
    return ExcitEnergy;
  }

  inline G4double GetFragmentMass(void) const
  {
    return fragmentMass;
  }

  inline G4double GetTotalEnergy(void) const
  {
    return (GetFragmentMass() + GetExcitationEnergy());
  }

  inline G4bool IsStable() const
  {
    return isStable;
  }

protected:

  G4bool isStable;

  G4int A;

  G4int Z;

  G4int Polarization;

  G4double ExcitEnergy;

  G4double fragmentMass;

};


#endif


