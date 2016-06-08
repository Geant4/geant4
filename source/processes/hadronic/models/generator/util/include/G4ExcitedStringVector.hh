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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4ExcitedStringVector.hh,v 1.6 2001/10/04 20:00:32 hpw Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
#ifndef G4ExcitedStringVector_h
#define G4ExcitedStringVector_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      ---------------- G4ExcitedStringVector ----------------
//             by Gunter Folger, June 1998.
// ------------------------------------------------------------

#include "G4ExcitedString.hh"
#include "g4std/vector"

typedef G4std::vector<G4ExcitedString *> G4ExcitedStringVector;
struct DeleteString { void operator()(G4ExcitedString* aS){delete aS;} };

#endif
