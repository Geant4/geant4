// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QResonance.hh,v 1.2 1999-12-15 14:52:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4QResonance_h
#define G4QResonance_h 1

// -------------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QResonance ----------------
//             by Mikhail Kossov, Sept 1999.
//  class header for mesonic P-wave Resonances used by the CHIPS Model
// -------------------------------------------------------------------

//#include <iostream.h>
#include "Randomize.hh"
#include "G4QHadron.hh"

class G4QResonance : public G4QHadron
{
public:
 G4QResonance(G4int PDGcode, G4double maxMass);
 ~G4QResonance();
 virtual G4double CalculateMass(G4double maxMass,G4int PDG); // Calculate the Resonance Mass
};

#endif
