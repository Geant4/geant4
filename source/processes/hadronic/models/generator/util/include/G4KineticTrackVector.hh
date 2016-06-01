// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4KineticTrackVector.hh,v 1.2 1998/10/09 13:46:50 maxim Exp $
// GEANT4 tag $Name: geant4-00 $

// Modified at 8-Oct-1998 by Maxim Komogorov. Methods BoostBeam,Boost,Shift
// were added.

#ifndef G4KineticTrackVector_h
#define G4KineticTrackVector_h 1

#include "globals.hh"
#include "G4KineticTrack.hh"
#include <rw/tpordvec.h>

class G4KineticTrackVector : public RWTPtrOrderedVector<G4KineticTrack>
    {
public:
  G4KineticTrackVector();

public:
    void BoostBeam(G4ThreeVector& BeamMom);
    void Boost(G4ThreeVector& Velocity);
    void Shift(G4ThreeVector& Pos);
    };

#endif

