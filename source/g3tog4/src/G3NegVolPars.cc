// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3NegVolPars.cc,v 1.3 1999-05-12 08:09:51 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "globals.hh"
#include "G4VSolid.hh"
#include "G3toG4.hh"
#include "G3VolTable.hh"
#include <math.h>

G4bool G3CalcParamsFn(G4double* Rpar, G4int npar, G4double* Rparm,
		      G4String shape, G4String shapem);

// pars in case of gsposp are those passed to gsposp; for other cases,
// they come from the logical volume
G4bool G3NegVolPars(G4double pars[], G4int *nparpt,
		    G4String name, G4String moth,
		    char routine[]){
  G4bool NegPresent = FALSE;
  // Retrieve the parameters of the volume
  /*
  G4String shape;
  G4int nmed, npar;
  G4double *Rpar;
  G4VSolid *solid = NULL;
  G3Vol.GetLVPars(&name, &shape, &nmed, &Rpar, &npar, solid);
  if (pars == NULL) {
    pars = Rpar;
    *nparpt = npar;
  } else {
    Rpar = pars;
  }
  // and of the mother
  G4String shapem;
  G4int nparm;
  G4double *Rparm;
  G4VSolid *solidm = NULL;
  G3Vol.GetLVPars(&moth, &shapem, &nmed, &Rparm, &nparm, solidm);

  if (strcmp(routine,"GSPOS") == 0 || strcmp(routine,"GSVOLU") == 0) {
    NegPresent = G3CalcParamsFn(Rpar, npar, Rparm, shape, shapem);
  }
  if (strcmp(routine,"GSDVN") == 0) {
    // just set the flag. The parametrization function figures out
    // what to do.
    for (G4int i=0;i<npar;i++) {
      if (Rpar[i] < 0) {
	NegPresent = TRUE;
      }
    }
  }
  */
  return NegPresent;
}



