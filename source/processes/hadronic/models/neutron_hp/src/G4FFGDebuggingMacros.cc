/*
 * File:   G4FFGDebuggingMacros.cc
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on August 17, 2012, 17:46
 */
 
#include "G4FFGDebuggingMacros.hh"

// Only set the values if G4DEBUG_VERBOSE is defined.
#ifdef G4DEBUG_VERBOSE
    G4long G4FFGDEBUG_RECURSIVE_REPRESSION = 0;
    G4long G4FFGDEBUG_DATA_STRUCTURE_REPRESSION = 0;
    G4long G4FFGDEBUG_RECURSIVE_REPRESSION_COUNTER = 0;
    G4long G4FFGDEBUG_DATA_STRUCTURE_REPRESSION_COUNTER = 0;
#endif /* G4DEBUG_VERBOSE */

