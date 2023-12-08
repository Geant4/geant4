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
// stores G4StackedTrack class objects into G4SubEvent

// Author: Makoto Asai (JLab) - 23/Aug/23
// --------------------------------------------------------------------
#ifndef G4SubEventTrackStack_hh
#define G4SubEventTrackStack_hh 1

#include "G4Types.hh"
#include "G4SubEvent.hh"

class G4Event;
class G4StackedTrack;

class G4SubEventTrackStack 
{
  public:

    G4SubEventTrackStack() = default;
    explicit G4SubEventTrackStack(G4int ty, std::size_t maxEnt)
      : fSubEventType(ty), fMaxEnt(maxEnt) {;}
   ~G4SubEventTrackStack();
  
    G4SubEventTrackStack& operator=(const G4SubEventTrackStack&) = delete;
    G4bool operator==(const G4SubEventTrackStack&) const = delete;
    G4bool operator!=(const G4SubEventTrackStack&) const = delete;
  
    void PushToStack(const G4StackedTrack& aStackedTrack);
    void ReleaseSubEvent();
    //G4StackedTrack PopFromStack();
  
    void PrepareNewEvent(G4Event* ev);
    void clearAndDestroy();

    inline G4int GetSubEventType() const { return fSubEventType; }
    inline G4int GetNTrack() const
    { 
      if(fCurrentSE==nullptr) return 0;
      return (G4int)fCurrentSE->size();
    }
    inline std::size_t GetMaxNTrack() const { return fMaxEnt; }
    inline void SetVerboseLevel(G4int val) { verboseLevel = val; }
  
  private:

    G4int fSubEventType = -1;
    std::size_t fMaxEnt = 1000;
    G4SubEvent* fCurrentSE = nullptr;
    G4Event* fCurrentEvent = nullptr;
    G4int verboseLevel = 0;
};

#endif
