// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VDigi.cc,v 1.3.4.1 1999/12/07 20:47:47 gunter Exp $
// GEANT4 tag $Name: geant4-03-00 $
//

// G4VDigi
#include "G4VDigi.hh"
#include "globals.hh"

G4VDigi::G4VDigi()
{;}

G4VDigi::~G4VDigi()
{;}

int G4VDigi::operator==(const G4VDigi &right) const
{ return false; }

void G4VDigi::Draw() 
{;}

void G4VDigi::Print() 
{;}

