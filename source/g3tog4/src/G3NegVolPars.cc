// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3NegVolPars.cc,v 1.7 2000-11-24 09:50:12 gcosmo Exp $ 
// GEANT4 tag $Name: not supported by cvs2svn $
//
// modified by I. Hrivnacova, 13.10.99 

#include "globals.hh"
#include "G3VolTable.hh"
#include "G4VSolid.hh"
#include "G3toG4.hh"
#include <math.h>

G4bool G3CalcParamsFn(G4double *Rpar, G4int npar, G4double *Rparm,
                      G4String shape, G4String shapem)
  // Returns true only in case the parameters *after* processing
  // this method remain negative.
{
  G4bool NegPresent = FALSE;
  // for normal single volume positioning, just substitute for the
  // negative parameters
  // treat only the legal cases
  if (shapem == shape) {
    if (shape == "BOX" || shape == "TRD1" || shape == "TRD2" || 
        shape == "ELTU") {
      for (G4int i=0;i<npar;i++) {
        if (Rpar[i] < 0) {
          if (Rparm != 0) Rpar[i] = Rparm[i];
	  if (Rpar[i] < 0) NegPresent = TRUE;
        }
      }
    }
    if (shape == "TRAP") {
      for (G4int i=0;i<11;i++) {
        if (i != 1 && i != 2 && i != 6 && i != 10) {
          if (Rpar[i]<0) {
            if (Rparm != 0) Rpar[i] = Rparm[i];
	    if (Rpar[i] < 0) NegPresent = TRUE;
	  }  
        }
      }
    }
    if (shape == "TUBE" || shape == "TUBS" || shape == "PARA") {
      for (G4int i=0;i<3;i++) {
        if (Rpar[i] < 0) {
          if (Rparm != 0) Rpar[i] = Rparm[i];
	  if (Rpar[i] < 0) NegPresent = TRUE;
        }
      }
    }
    if (shape == "CONE" || shape == "CONS") {
      for (G4int i=0;i<5;i++) {
        if (Rpar[i] < 0) {
          if (Rparm != 0) Rpar[i] = Rparm[i];
	  if (Rpar[i] < 0) NegPresent = TRUE;
        }
      }
    }
    if (shape == "SPHE") {
      for (G4int i=0;i<2;i++) {
        if (Rpar[i] < 0) {
          if (Rparm != 0) Rpar[i] = Rparm[i];
	  if (Rpar[i] < 0) NegPresent = TRUE;
        }
      }
    }
    if (shape == "PGON") {
      G4int nz = int(Rpar[3]);
      G4int ipl;
      for (G4int i=0;i<nz;i++) {
        ipl = 5 + i*3;
        if (Rpar[ipl] < 0) {
          if (Rparm != 0) Rpar[ipl] = Rparm[ipl];
	  if (Rpar[ipl] < 0) NegPresent = TRUE;
        }
        if (Rpar[ipl+1] < 0) {
          if (Rparm != 0)  Rpar[ipl] = Rparm[ipl];
	  if (Rpar[ipl] < 0) NegPresent = TRUE;
        }
      }
    }
    if (shape == "PCON") {
      G4int nz = int(Rpar[2]);
      G4int ipl;
      for (G4int i=0;i<nz;i++) {
        ipl = 4 + i*3;
        if (Rpar[ipl] < 0) {
          if (Rparm != 0) Rpar[ipl] = Rparm[ipl];
	  if (Rpar[ipl] < 0) NegPresent = TRUE;
        }
        if (Rpar[ipl+1] < 0) {
          // TO DO
	  // check - folowing argument might be ipl+1
          if (Rparm != 0) Rpar[ipl] = Rparm[ipl];
	  if (Rpar[ipl] < 0) NegPresent = TRUE;
        }
      }
    }
  }

  if (shape == "BOX") {
    if (shapem == "TRD1") {
      if (Rpar[1] < 0) {
        if (Rparm != 0) Rpar[1] = Rparm[2];
	if (Rpar[1] < 0) NegPresent = TRUE;
      }
      if (Rpar[2] < 0) {
        if (Rparm != 0) Rpar[2] = Rparm[3];
	if (Rpar[2] < 0) NegPresent = TRUE;
      }
      if (Rpar[0] < 0) {
        if (Rparm != 0) Rpar[0] = G4std::min(Rparm[0],Rparm[1]) +
                                  abs(Rparm[0]-Rparm[1])*.5*Rpar[2]/Rparm[3];
	if (Rpar[0] < 0) NegPresent = TRUE;
      }
    }
    if (shapem == "TRD2") {
      if (Rpar[2] < 0) {
        if (Rparm != 0) Rpar[2] = Rparm[4];
        if (Rpar[2]<0) NegPresent = TRUE;
      }
      if (Rpar[0] < 0) {
        if (Rparm != 0) Rpar[0] = G4std::min(Rparm[0],Rparm[1]) +
                                  abs(Rparm[0]-Rparm[1])*.5*Rpar[2]/Rparm[4];
        if (Rpar[0]<0) NegPresent = TRUE;
      }
      if (Rpar[1] < 0) {
        if (Rparm != 0) Rpar[1] = G4std::min(Rparm[2],Rparm[3]) +
                                  abs(Rparm[2]-Rparm[3])*.5*Rpar[2]/Rparm[4];
        if (Rpar[1]<0) NegPresent = TRUE;
      }
    }
    if (shapem == "TRAP") {
      if (Rpar[2] < 0) {
        if (Rparm != 0) Rpar[2] = Rparm[0];
        if (Rpar[2] < 0) NegPresent = TRUE;
      }
      if (Rpar[0] < 0) {
        if (Rparm != 0) {
          G4double xlo = G4std::min(Rparm[4],Rparm[8]) +
            abs(Rparm[4]-Rparm[8])*.5*Rpar[2]/Rparm[0];
          G4double xhi = G4std::min(Rparm[5],Rparm[9]) +
            abs(Rparm[5]-Rparm[9])*.5*Rpar[2]/Rparm[0];
          Rpar[0] = G4std::min(xlo,xhi);
        }
        if (Rpar[0] < 0) NegPresent = TRUE;
      }
      if (Rpar[1] < 0) {
        if (Rparm != 0) Rpar[1] = G4std::min(Rparm[3],Rparm[7]) +
                        abs(Rparm[3]-Rparm[7])*.5*Rpar[2]/Rparm[0];
        if (Rpar[1] < 0) NegPresent = TRUE;
      }
    }
  }
  return NegPresent;
}

G4bool G3NegVolPars(G4double pars[], G4int *nparpt,
                    G3VolTableEntry* vte,
                    G3VolTableEntry* mvte, const char routine[])
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
