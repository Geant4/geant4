// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gsatt.cc,v 1.2 1999-05-12 08:09:59 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include <rw/cstring.h>
#include "G3toG4.hh"
#include "G3VolTable.hh"

void PG4gsatt(RWCString tokens[])
{
    // fill the parameter containers
    G3fillParams(tokens,PTgsatt);

    // interpret the parameters
    G4String name = Spar[0];
    G4String attr = Spar[1];
    G4int ival = Ipar[0];

    G4gsatt(name, attr, ival);
}

void G4gsatt(G4String name, G4String attr, G4int ival)
{
    // get logical volume pointer
    G4LogicalVolume *lvol = G3Vol.GetLV(name);
    G4cerr << "G4gsatt not implemented" << endl;
    // apply attribute
// $$$    lvol->ApplyAttribute(attr, ival);
}
