// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gsdet.cc,v 1.1 1999-01-07 16:06:48 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "globals.hh"
#include "G3toG4.hh"

void PG4gsdet(RWCString tokens[])
{
    // fill the parameter containers
    G3fillParams(tokens,PTgsdet);

    // interpret the parameters
    G4String chset = Spar[0];
    G4String chdet = Spar[1];
    G4int nv = Ipar[0];
    G4String chnmsv[100];
    for (G4int i=0; i<=nv; i++ ) chnmsv[i] = Spar[2+i].data();
    G4int *nbits = &Ipar[1];
    G4int idtyp = Ipar[1+nv];
    G4int nwhi = Ipar[2+nv];
    G4int nwdi = Ipar[3+nv];

    G4gsdet(chset,chdet,nv,chnmsv,nbits,idtyp,nwhi,nwdi);
}

void G4gsdet(G4String chset, G4String chdet, G4int, G4String*,
             G4int*, G4int idtyp, G4int nwhi, G4int nwdi)
{ 
    G4gsdetv(chset, chdet, idtyp, nwhi, nwdi);
}
