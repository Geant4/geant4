//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
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

