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
// $Id: G4FermiPair.hh 85677 2014-11-03 17:44:12Z vnivanch $
//
// FermiBreakUp de-excitation model
// by V. Ivanchenko (July 2016)
//

#ifndef G4FermiPair_h
#define G4FermiPair_h 1

#include "globals.hh"
#include "G4FermiFragment.hh"
#include "G4Fragment.hh"

class G4FermiPair 
{
public:

  explicit G4FermiPair(const G4FermiFragment* f1, const G4FermiFragment* f2);

  inline G4int GetA() const;
  inline G4int GetZ() const;
  inline G4double GetMass() const;
  inline G4double GetExcitationEnergy() const;
  inline G4double GetTotalEnergy() const;
  inline G4double GetDynamicMinMass(G4double Eex) const;
  inline const G4FermiFragment* GetFragment1() const;
  inline const G4FermiFragment* GetFragment2() const;
  
private:

  inline G4FermiPair(const G4FermiPair &) = delete;
  inline const G4FermiPair & operator=(const G4FermiPair &) = delete;
  inline G4bool operator==(const G4FermiPair &) const = delete;
  inline G4bool operator!=(const G4FermiPair &) const = delete;
  
  G4int totalZ;
  G4int totalA;

  G4double mass;
  G4double excitEnergy;

  const G4FermiFragment* fragment1;
  const G4FermiFragment* fragment2;

};

inline G4int G4FermiPair::GetA() const
{
  return totalA;
}

inline G4int G4FermiPair::GetZ() const
{
  return totalZ;
}

inline G4double G4FermiPair::GetMass() const
{
  return mass;
}

inline G4double G4FermiPair::GetExcitationEnergy() const
{
  return excitEnergy;
}

inline G4double G4FermiPair::GetTotalEnergy() const
{
  return mass + excitEnergy;
}

inline const G4FermiFragment* G4FermiPair::GetFragment1() const
{
  return fragment1;
}

inline const G4FermiFragment* G4FermiPair::GetFragment2() const
{
  return fragment2;
}

inline G4double G4FermiPair::GetDynamicMinMass(G4double Eex) const
{
  return fragment1->GetTotalEnergy() + fragment2->GetTotalEnergy()
    + fragment1->GetCoulombBarrier(fragment2->GetA(), fragment2->GetZ(), Eex);
}

#endif


