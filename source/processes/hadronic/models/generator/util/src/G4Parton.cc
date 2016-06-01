// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Parton.cc,v 1.1 1998/08/22 08:58:11 hpw Exp $
// GEANT4 tag $Name: geant4-00 $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4Parton ----------------
//             by Gunter Folger, June 1998.
//       class for Parton (inside a string) used by Parton String Models
// ------------------------------------------------------------

#include "G4Parton.hh"

G4Parton::G4Parton(G4int PDGcode)
{
	PDGencoding=PDGcode;
}

G4Parton::~G4Parton()
{}
