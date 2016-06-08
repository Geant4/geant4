// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: rstream.h,v 1.4 1999/11/25 10:14:45 gcosmo Exp $
// GEANT4 tag $Name: geant4-01-01 $
//
// 
//---------------------------------------------------------------
//  GEANT 4 utility file
//
//  Definition of utility function rwEatwhite.
//---------------------------------------------------------------

#ifndef __rstream
#define __rstream

#include "g4rw/defs.h"
#include "G4Types.hh"
#include "g4std/iostream"

inline G4std::istream& rwEatwhite(G4std::istream& stream)
{
  return stream >> G4std::ws;
}

#endif
