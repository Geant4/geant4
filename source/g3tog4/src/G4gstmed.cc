// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gstmed.cc,v 1.3 1999-11-15 10:39:40 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4LogicalVolume.hh"
#include "G3toG4.hh"
#include "G3MatTable.hh"
#include "G3MedTable.hh"

void PG4gstmed(G4String tokens[])
{
    // fill the parameter containers
    G3fillParams(tokens,PTgstmed);

    // interpret the parameters
    G4String name = Spar[0];
    G4int itmed = Ipar[0];
    G4int nmat = Ipar[1];
    G4int isvol = Ipar[2];
    G4int ifield = Ipar[3];
    G4int nwbuf = Ipar[4];
    G4double fieldm = Rpar[0];
    G4double tmaxfd = Rpar[1];
    G4double stemax = Rpar[2];
    G4double deemax = Rpar[3];
    G4double epsil = Rpar[4];
    G4double stmin = Rpar[5];
    G4double *ubuf = &Rpar[6];

    G4gstmed(itmed,name,nmat,isvol,ifield,fieldm,tmaxfd,stemax,
             deemax,epsil,stmin,ubuf,nwbuf);
}

void G4gstmed(G4int itmed, G4String, G4int nmat, G4int isvol,
              G4int ifield, G4double fieldm, G4double tmaxfd,
              G4double stemax, G4double deemax, G4double epsil,
              G4double stmin, G4double* , G4int )
{
    // get the pointer to material nmat
    G4Material* material = G3Mat.get(nmat);
    // Store this medium in the G3Med structure
    G3Med.put(itmed,material);
}
