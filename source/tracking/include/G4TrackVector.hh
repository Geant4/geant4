// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TrackVector.hh,v 1.3 1999-10-14 05:39:44 tsasaki Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
//  G4TrackVector.hh
//
//  Description:
//    This class keeps a List of G4Track objects. It is implemented 
//    as a RougeWave pointer ordered vector.
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
//---------------------------------------------------------------

#ifndef G4TrackVector_h
#define G4TrackVector_h 1

#include <rw/tpordvec.h>           // Include from 'RogueWave'
#include "G4Track.hh"              // Include form 'tracking'

///////////////////////////////////////////////////
typedef RWTPtrOrderedVector<G4Track> G4TrackVector;
///////////////////////////////////////////////////

#endif


