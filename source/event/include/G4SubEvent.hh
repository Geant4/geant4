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
// This is a stack class used for transferring a sub-event.

// Author: Makoto Asai (JLAB) - 23/Aug/23
// --------------------------------------------------------------------
#ifndef G4SubEvent_hh
#define G4SubEvent_hh 1

#include <vector>

class G4Event;

#include "G4StackedTrack.hh"
#include "G4Types.hh"

class G4SubEvent : public std::vector<G4StackedTrack>
{
  public:

    G4SubEvent() = default;
    explicit G4SubEvent(G4int ty, std::size_t maxEnt)
      : fSubEventType(ty), fMaxEnt(maxEnt) {;}
   ~G4SubEvent();
  
    G4SubEvent& operator=(const G4SubEvent&) = delete;
    G4bool operator==(const G4SubEvent&) const = delete;
    G4bool operator!=(const G4SubEvent&) const = delete;
  
    inline void PushToStack(const G4StackedTrack& aStackedTrack)
    { push_back(aStackedTrack); }
    G4StackedTrack PopFromStack();
  
    void clearAndDestroy();

    inline G4int GetSubEventType() const { return fSubEventType; }
    inline std::size_t GetNTrack() const { return size(); }
    inline std::size_t GetMaxNTrack() const { return fMaxEnt; }
  
    G4double getTotalEnergy() const;

    // pointer to G4Event where this sub-event belongs to
    inline void SetEvent(G4Event* evt) { fpEvent = evt; }
    inline G4Event* GetEvent() const { return fpEvent; }
  
  private:

    G4int fSubEventType = -1;
    std::size_t fMaxEnt = 1000;
    G4Event* fpEvent = nullptr;

};

#endif
