// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3NegVolPars.cc,v 1.6 2000-03-07 10:51:40 stesting Exp $ 
// GEANT4 tag $Name: not supported by cvs2svn $
//
// modified by I. Hrivnacova, 13.10.99 

#include "globals.hh"
#include "G3VolTable.hh"
#include "G4VSolid.hh"
#include "G3toG4.hh"
#include <math.h>

G4bool G3CalcParamsFn(G4double* Rpar, G4int npar, G4double* Rparm,
                                  G4String shape, G4String shapem);

G4bool G3NegVolPars(G4double pars[], G4int *nparpt, 
          G3VolTableEntry* vte, G3VolTableEntry* mvte,  const char routine[])
{
  G4bool NegPresent = FALSE;

  // retrieve parameters 

  // the volume
  G4String shape = vte->GetShape();
  G4double* Rpar = vte->GetRpar();
  G4int npar = vte->GetNpar();  
  if (npar ==0) {
    // no solid parameters are defined in vte
    npar = *nparpt;
    Rpar = pars;
  }
  else {  
    // solid parameters are already defined in vte
    // pars[], nparpt are ignored
    // TO DO: check if g3 ignores them too or resets
    // vte parameters according to this new ones !!
  }  
      
  // mother
  G4String shapem = mvte->GetShape();
  G4double* Rparm = mvte->GetRpar();
  G4int nparm = mvte->GetNpar();    

  if (strcmp(routine,"GSPOS") == 0 || strcmp(routine,"GSVOLU") == 0) {
    NegPresent = G3CalcParamsFn(Rpar,npar,Rparm,shape,shapem);
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
  return NegPresent;
}
