// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4StackedTrack.hh,v 1.4 2000-01-26 06:42:15 asaim Exp $
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
      int operator==(const G4StackedTrack &right) const;
      int operator!=(const G4StackedTrack &right) const;

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

extern G4Allocator<G4StackedTrack> aStackedTrackAllocator;

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

