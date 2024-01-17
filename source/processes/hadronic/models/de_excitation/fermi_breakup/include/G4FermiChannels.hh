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

#ifndef G4FermiChannels_h
#define G4FermiChannels_h 1

#include "globals.hh"
#include "G4FermiFragment.hh"
#include "G4FermiPair.hh"
#include <vector>

class G4FermiChannels
{
public:

  explicit G4FermiChannels(const G4FermiFragment* ptr) : frag(ptr) {}
  ~G4FermiChannels() { for (auto const & p : fvect) { delete p; } }

  G4int GetZ() const { return frag->GetZ(); }
  const G4FermiFragment* GetFragment() const { return frag; }

  std::size_t NumberPairs() const { return nch; }
  const std::vector<G4FermiPair*>& GetChannels() const { return fvect; }
  const G4FermiPair* GetPair(std::size_t idx) const {
    return (idx < nch) ? fvect[idx] : nullptr;
  }

  void AddChannel(G4FermiPair* ptr) {
    fvect.push_back(ptr);
    ++nch;
  }

  G4double GetExcitation() const { return frag->GetExcitationEnergy(); }
  G4double GetMass() const { return frag->GetTotalEnergy(); }

  inline const G4FermiChannels& operator=(const G4FermiChannels&) = delete;
  inline G4FermiChannels(const G4FermiChannels &) = delete;
  inline G4bool operator==(const G4FermiChannels &) = delete;
  inline G4bool operator!=(const G4FermiChannels &) = delete;

private:

  const G4FermiFragment* frag;
  std::size_t nch{0};
  std::vector<G4FermiPair*> fvect;
};

#endif
