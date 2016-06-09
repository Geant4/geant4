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
// $Id: G4ListOfCollisions.hh,v 1.1.4.1 2004/03/24 13:18:28 hpw Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// $Id: G4ListOfCollisions.hh,v 1.0 1998/06/30
// -----------------------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, A. Feliciello, 30th June 1998
//               16 Nov 1999  M.G. Pia  STL migration
// -----------------------------------------------------------------------------

#ifndef G4ListOfCollisions_h
#define G4ListOfCollisions_h 1

#include "globals.hh"
#include <vector>
#include "G4CollisionInitialState.hh"

typedef std::vector<G4CollisionInitialState*> G4ListOfCollisions;
struct DeleteCollisionInitialState { void operator()(G4CollisionInitialState *aC) {delete aC;} };

#endif
