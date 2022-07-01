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
// G4SmartTrackStack
//
// Class description:
//
// This is a 'smart' stack class used by G4StackManager. This class object
// stores G4StackedTrack class objects in various dedicated stacks

// Author: S.Kamperis - 4 October 2012
// --------------------------------------------------------------------
#ifndef G4SmartTrackStack_hh
#define G4SmartTrackStack_hh 1

#include <array>

#include "G4StackedTrack.hh"
#include "G4TrackStack.hh"
#include "globals.hh"

class G4SmartTrackStack
{
  public:

    G4SmartTrackStack();
   ~G4SmartTrackStack();

    G4SmartTrackStack & operator=(const G4SmartTrackStack&) = delete;
    G4bool operator==(const G4SmartTrackStack&) const = delete;
    G4bool operator!=(const G4SmartTrackStack&) const = delete;

    void PushToStack(const G4StackedTrack& aStackedTrack);
    G4StackedTrack PopFromStack();
    void clear();
    void clearAndDestroy();
    void TransferTo(G4TrackStack* aStack);
    G4double getEnergyOfStack(G4TrackStack* aTrackStack);
    void dumpStatistics();

    inline G4int GetNTrack() const { return nTracks; }
    inline G4int GetMaxNTrack() const { return maxNTracks; }

  private:
    G4int fTurn = 0;
    static constexpr G4int nTurn{5};
    std::array<G4double,nTurn> energies;
    std::array<G4TrackStack*,nTurn> stacks;
    // = 0 : all primaries and secondaries except followings
    // = 1 : secondary neutrons
    // = 2 : secondary electrons
    // = 3 : secondary gammas
    // = 4 : secondary positrons
    G4int maxNTracks{0};
    G4int nTracks{0};
};

#endif
