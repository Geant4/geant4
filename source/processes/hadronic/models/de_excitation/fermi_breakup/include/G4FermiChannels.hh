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
// $Id: G4FermiChannels.hh 85677 2014-11-03 17:44:12Z vnivanch $
//
// FermiBreakUp de-excitation model
// by V. Ivanchenko (July 2016)
//

#ifndef G4FermiChannels_h
#define G4FermiChannels_h 1

#include "globals.hh"
#include "G4FermiFragment.hh"
#include "G4FermiPair.hh"
#include <vector>

class G4FermiChannels 
{
public:

  explicit G4FermiChannels(size_t nmax, G4double ex, G4double gmass) 
    : nch(0), excitation(ex), ground_mass(gmass) 
  { fvect.reserve(nmax); cum_prob.reserve(nmax); };

  inline size_t GetNumberOfChannels() const;
  inline const std::vector<const G4FermiPair*>& GetChannels() const;
  inline const G4FermiPair* GetPair(size_t idx) const;
  inline const G4FermiPair* SamplePair(G4double rand) const;

  inline void AddChannel(const G4FermiPair*);
  inline std::vector<G4double>& GetProbabilities();

  inline G4double GetExcitation() const;
  inline G4double GetMass() const;

private:

  inline const G4FermiChannels& operator=(const G4FermiChannels&) = delete;
  inline G4FermiChannels(const G4FermiChannels &) = delete;
  inline G4bool operator==(const G4FermiChannels &) const = delete;
  inline G4bool operator!=(const G4FermiChannels &) const = delete;

  size_t nch;
  G4double excitation;
  G4double ground_mass;
  std::vector<const G4FermiPair*> fvect;
  std::vector<G4double> cum_prob;

};

inline size_t G4FermiChannels::GetNumberOfChannels() const
{
  return nch;
}

inline const std::vector<const G4FermiPair*>& G4FermiChannels::GetChannels() const
{
  return fvect;
}

inline const G4FermiPair* G4FermiChannels::GetPair(size_t idx) const
{
  return (idx < nch) ? fvect[idx] : nullptr;
}

inline const G4FermiPair* G4FermiChannels::SamplePair(G4double rand) const
{
  const G4FermiPair* ptr = nullptr;
  for(size_t i=0; i<nch; ++i) { 
    if(rand <= cum_prob[i]) { ptr = fvect[i]; break; }
  }
  return ptr;
}

inline void G4FermiChannels::AddChannel(const G4FermiPair* ptr)
{
  fvect.push_back(ptr);
  cum_prob.push_back(1.0);
  ++nch;
}

inline std::vector<G4double>& G4FermiChannels::GetProbabilities()
{
  return cum_prob;
}

inline G4double G4FermiChannels::GetExcitation() const
{
  return excitation;
}

inline G4double G4FermiChannels::GetMass() const
{
  return excitation + ground_mass;
}

#endif


