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
// Class description:
//
// This is a stack class used by G4StackManager. This class object
// stores G4StackedTrack class objects in the form of bi-directional
// linked list.

// Author: Makoto Asai (SLAC) - 09/Dec/96
// --------------------------------------------------------------------
#ifndef G4TrackStack_hh
#define G4TrackStack_hh 1

#include <vector>

#include "G4StackedTrack.hh"
#include "G4Types.hh"

class G4SmartTrackStack;

class G4TrackStack : public std::vector<G4StackedTrack>
{
  public:

    G4TrackStack() = default;
    explicit G4TrackStack(std::size_t n)
      : safetyValue1(G4int(4*n/5)),
        safetyValue2(G4int(4*n/5-100)), nstick(100) { reserve(n); }
   ~G4TrackStack();
  
    G4TrackStack& operator=(const G4TrackStack&) = delete;
    G4bool operator==(const G4TrackStack&) const = delete;
    G4bool operator!=(const G4TrackStack&) const = delete;
  
    inline void PushToStack(const G4StackedTrack& aStackedTrack)
      { push_back(aStackedTrack); }
    inline G4StackedTrack PopFromStack()
      { G4StackedTrack st = back(); pop_back(); return st; }
    void TransferTo(G4TrackStack* aStack);
    void TransferTo(G4SmartTrackStack* aStack);
  
    void clearAndDestroy();

    inline std::size_t GetNTrack() const { return size(); }
    inline std::size_t GetMaxNTrack() const { return max_size(); }
    inline G4int GetSafetyValue1() const { return safetyValue1; }
    inline G4int GetSafetyValue2() const { return safetyValue2; }
    inline G4int GetNStick() const { return nstick; }
  
    G4double getTotalEnergy() const;
    inline void SetSafetyValue2(G4int x) { safetyValue2 = x  < 0 ? 0 : x; }
  
  private:

    G4int safetyValue1{0};
    G4int safetyValue2{0};
    G4int nstick{0};
};

#endif
