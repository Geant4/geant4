// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gsdvx.cc,v 1.3 1999-12-05 17:50:13 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by I.Hrivnacova, V.Berejnoi, 27 Sep 99

#include "globals.hh"
#include "G3toG4.hh"

void G4gsdvn2(G4String name, G4String moth, G4int ndiv, G4int iaxis,
              G4double c0, G4int numed);

void G4gsdvt2(G4String name, G4String moth, G4double Step, G4int iaxis,
              G4double c0, G4int numed, G4int ndvmx);

void PG4gsdvx(G4String tokens[])
{
    // fill the parameter containers
    G3fillParams(tokens,PTgsdvx);

    // interpret the parameters
    G4String name = Spar[0];
    G4String moth = Spar[1];
    G4int ndiv = Ipar[0];
    G4int iaxis = Ipar[1];
    G4int numed = Ipar[2];
    G4int ndvmx = Ipar[3];
    G4double Step = Rpar[0];
    G4double c0 = Rpar[1];

    G4gsdvx(name,moth,ndiv,iaxis,Step,c0,numed,ndvmx);
}

void G4gsdvx(G4String name, G4String moth, G4int ndiv, G4int iaxis,
             G4double Step, G4double c0, G4int numed, G4int ndvmx)
{
    // pass to gsdvn2 or gsdvt2
    if (Step > 0.) {
        G4gsdvt2(name,moth,Step,iaxis,c0,numed,ndvmx);
    } else if (ndiv > 0) {
        G4gsdvn2(name,moth,ndiv,iaxis,c0,numed);
    }
}
