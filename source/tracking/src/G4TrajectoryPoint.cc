// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TrajectoryPoint.cc,v 1.3 1999-10-14 05:39:53 tsasaki Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// ---------------------------------------------------------------
//
// G4TrajectoryPoint.cc
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
// ---------------------------------------------------------------

#include "G4TrajectoryPoint.hh"

G4Allocator<G4TrajectoryPoint> aTrajectoryPointAllocator;

//////////////////////////////////////
G4TrajectoryPoint::G4TrajectoryPoint()
//////////////////////////////////////
{
  fPosition = G4ThreeVector(0.,0.,0.);
}

///////////////////////////////////////////////////////
G4TrajectoryPoint::G4TrajectoryPoint(G4ThreeVector pos)
///////////////////////////////////////////////////////
{
  fPosition = pos;
}

////////////////////////////////////////////////////////////////////
G4TrajectoryPoint::G4TrajectoryPoint(const G4TrajectoryPoint &right)
////////////////////////////////////////////////////////////////////
 : fPosition(right.fPosition)
{
}

///////////////////////////////////////
G4TrajectoryPoint::~G4TrajectoryPoint()
///////////////////////////////////////
{
}


