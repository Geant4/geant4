// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gstmed.cc,v 1.1 1999-01-07 16:06:51 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4LogicalVolume.hh"
#include "G3toG4.hh"
#include "G3MatTable.hh"
#include "G3MedTable.hh"

void PG4gstmed(RWCString tokens[])
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
    // NB. there is the possibility for redundancy in the mag field
    //     and user limits objects. Who cares.
    // Generate the mag field object
    G4MagneticField *field = NULL;
// $$$     G4MagneticField* field = new G4MagneticField(ifield, fieldm, tmaxfd);
    // Generate the user limits object
    G4UserLimits *limits = NULL;
// $$$    G4UserLimits* limits = new G4UserLimits(stmin, stemax);
    // Store this medium in the G3Med structure
    G3Med.put(itmed,material,field,limits,isvol,deemax,epsil);
}
