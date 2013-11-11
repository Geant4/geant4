/*
 * File:   G4FFGVerboseMacros.hh
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on August 24, 2012, 16:23
 */

#ifndef G4FFGVERBOSEMACROS_HH
#define G4FFGVERBOSEMACROS_HH

#include "globals.hh"

/** G4FFG_DEPTH is used to track the depth of the function calls in the
 *  fission fragment generator code.
 */
extern G4long G4FFG_DEPTH;

/** G4FFG_LOCATION__ outputs the current location in the code */
#define G4FFG_LOCATION__ \
G4String debugOutput(__FILE__); \
debugOutput = debugOutput.substr(debugOutput.find_last_of('/') + 1); \
G4cout << G4FFG_FUNCTION_SIGNATURE__ << " at " << debugOutput << ":" << __LINE__;

/** G4FFG_SPACING__ indents the debug messages according to the debugging depth */
#define G4FFG_SPACING__ \
for(G4int depth = 0; depth < G4FFG_DEPTH; depth++) \
{ \
    G4cout << "  "; \
}

#endif /* G4FFGVERBOSEMACROS_HH */

