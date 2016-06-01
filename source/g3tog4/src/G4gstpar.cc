// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gstpar.cc,v 2.0 1998/07/02 16:17:21 gunter Exp $
// GEANT4 tag $Name: geant4-00 $
//
#include "G3toG4.hh"
#include "G3VolTable.hh"

void PG4gstpar(RWCString tokens[])
{
    // fill the parameter containers
    G3fillParams(tokens,PTgstpar);

    // interpret the parameters
    G4String chpar = Spar[0];
    G4int itmed = Ipar[0];
    G4double parval = Rpar[0];

    G4gstpar(itmed,chpar,parval);
}

void G4gstpar(G4int itmed, G4String chpar, G4double parval)
{
    // set special tracking medium parameter. Apply to all logical
    // volumes making use of the specified tracking medium.
    G3Vol.MatchTmed(itmed);
    G4LogicalVolume *lvol;
    while (lvol = G3Vol.NextTmed()) {
        // $$$ apply tracking parameter
    }
}
