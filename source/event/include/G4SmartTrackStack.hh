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
// $Id: G4SmartTrackStack.hh 66892 2013-01-17 10:57:59Z gunter $
//
//
//  Last Modification : 04/Oct/12 S. Kamperis
//


#ifndef G4SmartTrackStack_h
#define G4SmartTrackStack_h 1

#include "G4StackedTrack.hh"
#include "G4TrackStack.hh"
#include "globals.hh"

// class description:
//
// This is a 'smart' stack class used by G4StackManager. This class object
// stores G4StackedTrack class objects in various dedicated stacks

class G4SmartTrackStack
{
  public:
      G4SmartTrackStack();
      ~G4SmartTrackStack();

  private:
      const G4SmartTrackStack & operator=
                          (const G4SmartTrackStack &right);
      G4int operator==(const G4SmartTrackStack &right) const;
      G4int operator!=(const G4SmartTrackStack &right) const;

  public:
      void PushToStack(const G4StackedTrack& aStackedTrack);
      G4StackedTrack PopFromStack();
      void clear();
      void clearAndDestroy();
      void TransferTo(G4TrackStack* aStack);
      G4double getEnergyOfStack(G4TrackStack* aTrackStack);
      void dumpStatistics();

  private:
      G4int fTurn;
      G4int nTurn; // should be 5
      G4double energies[5];
      G4TrackStack* stacks[5];
      // = 0 : all primaries and secondaries except followings
      // = 1 : secondary neutrons
      // = 2 : secondary electrons
      // = 3 : secondary gammas
      // = 4 : secondary positrons
      G4int maxNTracks;
      G4int nTracks;

  public:
      G4int GetNTrack() const { return nTracks; }
      G4int GetMaxNTrack() const { return maxNTracks; }

  private:
      inline G4int n_stackedTrack() const
      {
	      return stacks[0]->GetNTrack() +
		     stacks[1]->GetNTrack() +
		     stacks[2]->GetNTrack() +
		     stacks[3]->GetNTrack() +
		     stacks[4]->GetNTrack();
      }
};

#endif

