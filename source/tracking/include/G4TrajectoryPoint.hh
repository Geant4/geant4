// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TrajectoryPoint.hh,v 1.6 2000-11-11 06:34:10 tsasaki Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
// G4TrajectoryPoint.hh
//
// class description:
//   This class represents the trajectory of a particle tracked.
//   It includes information of 
//     1) List of trajectory points which compose the trajectory,
//     2) static information of particle which generated the 
//        trajectory, 
//     3) trackID and parent particle ID of the trajectory,
//     4) termination condition of the trajectory.
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
// ---------------------------------------------------------------

#ifndef G4TrajectoryPoint_h
#define G4TrajectoryPoint_h 1

#include "G4VTrajectoryPoint.hh"
#include "globals.hh"                // Include from 'global'
#include "G4ThreeVector.hh"          // Include from 'geometry'
#include "G4Allocator.hh"            // Include from 'particle+matter'


////////////////////////
class G4TrajectoryPoint : public G4VTrajectoryPoint
//////////////////////// 
{

//--------
public: // without description
//--------

// Constructor/Destructor
   G4TrajectoryPoint();
   G4TrajectoryPoint(G4ThreeVector pos);
   G4TrajectoryPoint(const G4TrajectoryPoint &right);
   virtual ~G4TrajectoryPoint();

// Operators
   inline void *operator new(size_t);
   inline void operator delete(void *aTrajectoryPoint);
   inline int operator==(const G4TrajectoryPoint& right) const
   { return (this==&right); };

// Get/Set functions
   inline const G4ThreeVector GetPosition() const
   { return fPosition; };


//---------
   private:
//---------

// Member data
   G4ThreeVector fPosition;

};


extern G4Allocator<G4TrajectoryPoint> aTrajectoryPointAllocator;

inline void* G4TrajectoryPoint::operator new(size_t)
{
   void *aTrajectoryPoint;
   aTrajectoryPoint = (void *) aTrajectoryPointAllocator.MallocSingle();
   return aTrajectoryPoint;
}

inline void G4TrajectoryPoint::operator delete(void *aTrajectoryPoint)
{
   aTrajectoryPointAllocator.FreeSingle((G4TrajectoryPoint *) aTrajectoryPoint);
}

#endif

