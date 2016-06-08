// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VParticipants.cc,v 1.1.8.1.2.1 1999/12/08 17:34:43 gunter Exp $
// GEANT4 tag $Name: geant4-01-01 $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4VParticipants ----------------
//             by Gunter Folger, May 1998.
//      abstract class finding participants in a hadron Nucleus collision
//       in Parton String Models.
// ------------------------------------------------------------

#include "G4VParticipants.hh"
#include "Randomize.hh"

G4VParticipants::G4VParticipants() : theNucleus(NULL)
{}



G4VParticipants::~G4VParticipants()
{
// G4cout << "G4VParticipants::~G4VParticipants()" << G4endl;
	if ( theNucleus != NULL ) delete theNucleus;
}







