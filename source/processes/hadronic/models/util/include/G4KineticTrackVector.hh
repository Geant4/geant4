//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//

// Modified at 8-Oct-1998 by Maxim Komogorov. Methods BoostBeam,Boost,Shift
// were added.

#ifndef G4KineticTrackVector_h
#define G4KineticTrackVector_h 1

#include "globals.hh"
#include "G4KineticTrack.hh"
#include <vector>

class G4KineticTrackVector : public std::vector<G4KineticTrack *>
    {
public:
  G4KineticTrackVector();

public:
    void BoostBeam(G4ThreeVector& BeamMom);
    void Boost(G4ThreeVector& Velocity);
    void Shift(G4ThreeVector& Pos);
    };

struct DeleteKineticTrack{void operator()(G4KineticTrack * aT){delete aT;}};

#endif

