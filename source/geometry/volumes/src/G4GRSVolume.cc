// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4GRSVolume.cc,v 2.0 1998/07/02 17:06:28 gunter Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// class G4GRSVolume Implementation

#include "G4GRSVolume.hh"

G4GRSVolume::~G4GRSVolume()
{
  delete frot;			// safe if null
}
