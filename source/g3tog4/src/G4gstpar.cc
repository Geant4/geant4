// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gstpar.cc,v 1.5 1999-12-15 14:49:44 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G3toG4.hh"
#include "G3VolTable.hh"

void PG4gstpar(G4String tokens[])
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
  G4cerr << "G4gstpar: not implemented." << G4endl;
}
