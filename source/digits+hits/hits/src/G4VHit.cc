// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VHit.cc,v 2.1 1998/07/12 02:53:44 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//

// G4VHit
#include "G4VHit.hh"
#include "globals.hh"

G4VHit::G4VHit()
{;}

G4VHit::~G4VHit()
{;}

int G4VHit::operator==(const G4VHit &right) const
{ return false; }

void G4VHit::Draw()
{;}

void G4VHit::Print()
{;}

