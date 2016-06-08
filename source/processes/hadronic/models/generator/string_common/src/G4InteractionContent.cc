// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4InteractionContent.cc,v 1.1.10.1 1999/12/07 20:51:54 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4InteractionContent----------------
//             by Gunter Folger, June 1998.
//       class for a storing two colliding particles in PartonString Models
// ------------------------------------------------------------

#include "G4InteractionContent.hh"

G4InteractionContent::G4InteractionContent(G4VSplitableHadron *aPrimaryParticipant)
      : theNumberOfHard(NULL), theNumberOfSoft(NULL)
{
	theProjectile=aPrimaryParticipant;
}

G4InteractionContent::G4InteractionContent(const G4InteractionContent &right)
     : theNumberOfHard(0), theNumberOfSoft(0)
{}


G4InteractionContent::~G4InteractionContent()
{}

int G4InteractionContent::operator==(const G4InteractionContent &right) const
{
	return this==&right;
}

int G4InteractionContent::operator!=(const G4InteractionContent &right) const
{
        return this!=&right;
}


const G4InteractionContent & G4InteractionContent::operator=(const G4InteractionContent &right)
{
	G4Exception("G4InteractionContent::operator= meant to not be accessable");
	return *this;
}


