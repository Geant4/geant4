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
// $Id: G4FermiFragmentsPoolVI.hh,v 1.5 2006-06-29 20:13:13 gunter Exp $
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
#include "G4FermiDecayProbability.hh"

#include <vector>

class G4FermiFragmentsPoolVI
{
public:

  explicit G4FermiFragmentsPoolVI();

  ~G4FermiFragmentsPoolVI();

  const G4FermiChannels* ClosestChannels(G4int Z, G4int A, G4double mass) const;

  void DumpFragment(const G4FermiFragment*) const;

  void Dump() const;
 
  G4bool IsApplicable(G4int ZZ, G4int AA, G4double etot) const;

  inline const G4FermiDecayProbability* FermiDecayProbability() const;

  inline G4int GetMaxZ() const;

  inline G4int GetMaxA() const;

  inline G4double GetEnergyLimit() const;

  inline G4double GetTolerance() const;
  
private:

  void Initialise();

  G4bool IsInThePool(G4int Z, G4int A, G4double exc) const;

  G4bool IsInPhysPairs(const G4FermiFragment* f1, const G4FermiFragment* f2,
		       G4double exc) const;

  G4bool IsInUnphysPairs(const G4FermiFragment* f1, const G4FermiFragment* f2,
			 G4double exc) const;

  G4int maxZ;
  G4int maxA;
  G4double tolerance;
  G4double elim;
  G4double elim_unstable;

  G4FermiDecayProbability theDecay;

  // pool 
  std::vector<const G4FermiFragment*> fragment_pool;
  std::vector<const G4FermiFragment*> funstable;

  // lists of configurations sorted by A 

  // "stable" fragments
  std::vector<const G4FermiFragment*> list_f[17]; 
  // list of channels for "stable" fragments
  std::vector<G4FermiChannels*> list_c[17];
  // pairs of stable fragments 
  std::vector<const G4FermiPair*>     list_p[17]; 

  // "unstable" fragments
  std::vector<const G4FermiFragment*> list_g[17]; 
  // list of channels of stable and unstable fragments 
  std::vector<G4FermiChannels*> list_d[17];
  // pairs of stable and unstable fragments 
  std::vector<const G4FermiPair*>     list_u[17]; 

};

inline G4int G4FermiFragmentsPoolVI::GetMaxZ() const
{
  return maxZ;
}

inline G4int G4FermiFragmentsPoolVI::GetMaxA() const
{
  return maxA;
}

inline const G4FermiDecayProbability* 
G4FermiFragmentsPoolVI::FermiDecayProbability() const
{
  return &theDecay;
}

inline G4double G4FermiFragmentsPoolVI::GetEnergyLimit() const
{
  return elim;
}

inline G4double G4FermiFragmentsPoolVI::GetTolerance() const
{
  return tolerance;
}

#endif

