// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TrackVector.hh,v 1.6 2000-11-11 06:34:10 tsasaki Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
//  G4TrackVector.hh
//
//  class description:
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

#include "g4rw/tpordvec.h"           // Include from 'RogueWave'
#include "G4Track.hh"              // Include form 'tracking'

///////////////////////////////////////////////////
typedef G4RWTPtrOrderedVector<G4Track> G4TrackVector;
///////////////////////////////////////////////////

#endif


