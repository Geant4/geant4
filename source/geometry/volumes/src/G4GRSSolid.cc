// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GRSSolid.cc,v 1.2 1999-12-15 14:50:26 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// class G4GRSSolid Implementation

#include "G4GRSSolid.hh"

G4GRSSolid::~G4GRSSolid()
{
  delete frot;			// safe if null
}
