// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gsvolu.cc,v 1.2 1999-05-06 04:27:23 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4ios.hh"
#include "G3toG4.hh"
#include "G3VolTable.hh"
#include "G3MedTable.hh"

G4LogicalVolume* G4makevol(G4String vname, G4String shape, G4int nmed,
                           G4double Rpar[], G4int npar);
G4bool G3CalcParamsFn(G4double *Rpar, G4int npar, G4double *Rparm,
                       G4String shape, G4String shapem);

void PG4gsvolu(RWCString tokens[]) {
    // fill the parameter containers
    G3fillParams(tokens,PTgsvolu);

    // interpret the parameters
    G4String vname = Spar[0];
    G4String shape = Spar[1];
    G4int nmed = Ipar[0];
    G4int npar = Ipar[1];
    G4double *pars = Rpar;

    G4gsvolu(vname, shape, nmed, pars, npar);
}

void G4gsvolu(G4String vname, G4String shape, G4int nmed, G4double* Rpar,
              G4int npar)
{
    G4double rangehi[3];
    G4double rangelo[3];
    G4bool isIndexed = npar == 0;
    G4LogicalVolume *lvol = 0;
    G4int nmax = 11;

    if (npar>11) nmax = npar;
    G4double *param = new G4double[nmax]; 
    for (G4int i=0; i<npar; i++) {
        param[i] = Rpar[i];
    }

    // If volume is indexed (npar=0, gsposp used, or negative length parameters)
    // cannot Build a G4LogicalVolume now.
    // Build a G3Vol table entry, with 0 pointer, and defer LV creation
    // to the gsposp call.
    if ( npar == 0 ) {
      isIndexed = true;
    } else {
      // Negative length parameters?
      if (G3CalcParamsFn(param, npar, 0, shape, shape)) {
        isIndexed = true;
        G4cout << "G4gsvolu: Volume " << vname << " type " << shape <<
          " is indexed via negative parameters" << endl;
        G4cout << "  Params: ";
        for (G4int i=0; i<npar; i++) G4cout << param[i] << " ";
        G4cout << endl;
      }
    }
    if ( isIndexed ) {
      G4VSolid *solid = 0;
      G3Vol.PutLV(&vname, lvol, nmed, shape, param, npar, solid);
    } else {
      // Else proceed with volume creation.
      lvol = G4makevol(vname, shape, nmed, param, npar);
    }
    delete [] param;
}



