// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Trajectory.hh,v 1.1 1999-01-07 16:14:30 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
// G4Trajectory.hh
//
// Description:
//   This class represents the trajectory of a particle tracked.
//   It includes information of 
//     1) List of trajectory points which compose the trajectory,
//     2) static information of particle which generated the 
//        trajectory, 
//     3) trackID and parent particle ID of the trajectory,
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Makoto  Asai   (e-mail: asai@kekvax.kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
// ---------------------------------------------------------------

class G4Trajectory;

#ifndef G4Trajectory_h
#define G4Trajectory_h 1

#include <stdlib.h>                 // Include from 'system'
#include "G4ios.hh"               // Include from 'system'
#include <rw/tvordvec.h>            // RWTValOrderedVector
#include "globals.hh"               // Include from 'global'
#include "G4ParticleDefinition.hh"  // Include from 'particle+matter'
#include "G4TrajectoryPoint.hh"     // Include from 'tracking'
#include "G4Track.hh"
#include "G4Step.hh"

class G4Polyline;                   // Forward declaration.

///////////////////
class G4Trajectory
///////////////////
{

//--------
   public:
//--------

// Constructor/Destrcutor
   G4Trajectory(G4Track* aTrack);
   G4Trajectory(G4Trajectory &);
   ~G4Trajectory();

// Operators
   inline int operator == (const G4Trajectory& right){return (this==&right);}; 

// Get/Set functions 
   inline G4int GetTrackID() const
   { return fTrackID; }
   inline G4int GetParentID() const
   { return fParentID; }
   inline G4String GetParticleName() const
   { return ParticleName; }
   inline G4double GetCharge() const
   { return PDGCharge; }
   inline G4int GetPDGEncoding() const
   { return PDGEncoding; }

// Other member functions
   void ShowTrajectory();
     // Print all information of the trajectory to stdout

   void DrawTrajectory(G4int i_mode=0);

   void AppendStep(G4Step* aStep);

   G4ParticleDefinition* GetParticleDefinition();

//---------
   private:
//---------

  RWTValOrderedVector<G4TrajectoryPoint> positionRecord;
  G4int fTrackID;
  G4int fParentID;
  G4String ParticleName;
  G4double PDGCharge;
  G4int    PDGEncoding;

};

#endif










