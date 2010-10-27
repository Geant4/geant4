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
// $Id: G4StackedTrack.hh,v 1.12 2010-10-27 07:21:13 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//  Last Modification : 02/Feb/96 M.Asai
//


#ifndef G4StackedTrack_h
#define G4StackedTrack_h 1

#include "G4Track.hh"
#include "globals.hh"
#include "G4Allocator.hh"
class G4VTrajectory;

// class description:
//
//  This class is exclusively used by G4StackManager and G4TrackStack
// classes for storing a G4Track object in the form of bi-directional
// linked list.

class G4StackedTrack 
{
  public:
      inline void *operator new(size_t);
      inline void operator delete(void *aStackedTrack);

      G4StackedTrack();
      G4StackedTrack(G4Track * aTrack, G4VTrajectory * aTrajectory = 0);
      ~G4StackedTrack();

      const G4StackedTrack & operator=(const G4StackedTrack &right);
      G4int operator==(const G4StackedTrack &right) const;
      G4int operator!=(const G4StackedTrack &right) const;

  private:
      G4double priorityWeight;
      G4Track * track;
      G4VTrajectory * trajectory;
      G4StackedTrack * previousStackedTrack;
      G4StackedTrack * nextStackedTrack;

  public:
      inline G4double GetPriorityWeight() const
      { return priorityWeight; }
      inline void SetPriorityWeight(const G4double value)
      { priorityWeight = value; }
      inline G4Track * GetTrack() const
      { return track; }
      inline G4VTrajectory * GetTrajectory() const
      { return trajectory; }
      inline G4StackedTrack * GetPrevious() const
      { return previousStackedTrack; }
      inline G4StackedTrack * GetNext() const
      { return nextStackedTrack; }
      inline void SetPrevious(G4StackedTrack * value)
      { previousStackedTrack = value; }
      inline void SetNext(G4StackedTrack * value)
      { nextStackedTrack = value; }
};

#if defined G4EVENT_ALLOC_EXPORT
  extern G4DLLEXPORT G4Allocator<G4StackedTrack> aStackedTrackAllocator;
#else
  extern G4DLLIMPORT G4Allocator<G4StackedTrack> aStackedTrackAllocator;
#endif

inline void * G4StackedTrack::operator new(size_t)
{
  void * aStackedTrack;
  aStackedTrack = (void *) aStackedTrackAllocator.MallocSingle();
  return aStackedTrack;
}

inline void G4StackedTrack::operator delete(void * aStackedTrack)
{
  aStackedTrackAllocator.FreeSingle((G4StackedTrack *) aStackedTrack);
}


#endif

