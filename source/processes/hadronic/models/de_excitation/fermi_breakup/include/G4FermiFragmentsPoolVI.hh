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

#ifndef G4FermiFragmentsPoolVI_hh 
#define G4FermiFragmentsPoolVI_hh 1

#include "globals.hh"
#include "G4FermiFragment.hh"
#include "G4FermiPair.hh"
#include "G4FermiChannels.hh"

#include <vector>

class G4FermiFragmentsPoolVI
{
public:

  G4FermiFragmentsPoolVI();

  ~G4FermiFragmentsPoolVI();

  void Initialise();

  const G4FermiChannels* ClosestChannels(const G4int Z, const G4int A,
                                         const G4double mass) const;

  void DumpFragment(const G4FermiFragment*) const;

  void Dump() const;

  G4bool HasDecay(const G4int Z, const G4int A, const G4double eexc) const;

  G4bool IsInitialized() const { return isInitialized; };

  G4FermiFragmentsPoolVI(const G4FermiFragmentsPoolVI &right) = delete;  
  const G4FermiFragmentsPoolVI & operator=
  (const G4FermiFragmentsPoolVI &right) = delete;
  G4bool operator==(const G4FermiFragmentsPoolVI &right) const = delete;
  G4bool operator!=(const G4FermiFragmentsPoolVI &right) const = delete;
  
private:

  G4bool IsInThePool(const G4int Z, const G4int A, const G4double exc) const;

  G4double fTolerance{0.0};
  G4double fElim{0.0};

  const G4int maxZ{9};
  const G4int maxA{17};

  G4bool isInitialized{false};

  // pool 
  std::vector<const G4FermiFragment*> fragment_pool;

  // list of channels sorted by Z and A
  std::vector<G4FermiChannels*>* list_c[9][17] = {{nullptr}};
};

#endif
