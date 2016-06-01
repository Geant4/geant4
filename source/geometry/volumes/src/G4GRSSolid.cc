// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GRSSolid.cc,v 2.0 1998/07/02 17:06:27 gunter Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// class G4GRSSolid Implementation

#include "G4GRSSolid.hh"

G4GRSSolid::~G4GRSSolid()
{
  delete frot;			// safe if null
}
