// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DiffractiveHHScatterer.hh,v 1.1.10.1 1999/12/07 20:51:42 gunter Exp $

#ifndef G4DiffractiveHHScatterer_h
#define G4DiffractiveHHScatterer_h 1
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4DiffractiveExcitation --------------
//             by Gunter Folger, October 1998.
//      diffractive Excitation used by strings models
//	Take a projectile and a target
//	excite the projectile and target
// ------------------------------------------------------------

#include "globals.hh"
class G4DiffractiveExcitation;
class G4LundStringFragmentation;
class G4KineticTrack;
#include "G4KineticTrackVector.hh"

class G4DiffractiveHHScatterer
{
public:

   G4DiffractiveHHScatterer();
   
   G4KineticTrackVector * Scatter(const G4KineticTrack & aTrack, const G4KineticTrack & bTrack);

private:

const G4DiffractiveExcitation * theExcitation;
G4LundStringFragmentation * theStringFragmentation;
};

#endif
