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
// FermiBreakUp de-excitation model
// by V. Ivanchenko (July 2016)
//

#ifndef G4FermiFragment_h
#define G4FermiFragment_h 1

#include "globals.hh"
#include "G4FragmentVector.hh"
#include "G4VCoulombBarrier.hh"

class G4FermiFragment 
{
public:

  explicit G4FermiFragment(G4int anA, G4int aZ, G4int sp, G4double exc);

  ~G4FermiFragment();

  inline G4int GetA() const { return A; }
  
  inline G4int GetZ() const { return Z; }
  
  inline G4int GetSpin() const { return spin; }

  inline G4double GetExcitationEnergy() const { return excitEnergy; }

  inline G4double GetFragmentMass() const { return fragmentMass; }

  inline G4double GetTotalEnergy(void) const 
  {
    return (fragmentMass + excitEnergy);
  }

  inline G4double GetCoulombBarrier(G4int Ares, G4int Zres, G4double Eex) const
  {
    return cBarrier->GetCoulombBarrier(Ares, Zres, Eex);
  }

  inline G4bool operator==(const G4FermiFragment &right) const
  {
    return (A == right.A && Z == right.Z &&
	    std::abs(excitEnergy - right.excitEnergy) < 0.0001);
  }
  
  G4FermiFragment(const G4FermiFragment &right) = delete;
  const G4FermiFragment & operator=(const G4FermiFragment &right) = delete;

private:

  G4double excitEnergy;
  G4double fragmentMass;

  G4VCoulombBarrier* cBarrier;

  G4int A;
  G4int Z;
  G4int spin;
};


#endif


