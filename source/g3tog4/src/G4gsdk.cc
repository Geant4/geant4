// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gsdk.cc,v 1.2 1999-11-15 10:39:38 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4Decay.hh"
#include "G3toG4.hh"
#include "G3PartTable.hh"

void PG4gsdk(G4String tokens[])
{
    // fill the parameter containers
    G3fillParams(tokens,PTgsdk);

    // interpret the parameters
    G4int ipart = Ipar[0];
    G4int inum = Ipar[1];
    G4int *mode = &Ipar[3];
    G4double *bratio = Rpar;

    G4gsdk(ipart,bratio,mode);
}

void G4gsdk(G4int ipart, G4double bratio[], G4int mode[])
{
    // create decay object for the particle
    G4Decay *decay = new G4Decay();
    // add decay modes
    for (G4int i=0; i<6; i++) {
        if (mode[i] != 0) {
// $$$            decay->AddMode(mode[i],bratio[i]);
        }
    }
    // associate decay object with particle ipart
    G4ParticleDefinition *part = G3Part.get(ipart);
// $$$    part->SetDecay(decay);
}
