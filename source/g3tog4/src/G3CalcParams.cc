// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3CalcParams.cc,v 1.3 1999-12-15 14:49:42 gunter Exp $ 
// GEANT4 tag $Name: not supported by cvs2svn $ 
//
// Change: G3CalcParamsFn returns true only in case
// the parameters *after* processing this method
// remain negative 
//
// modified by I.Hrivnacova, 27 Sep 99

#include "globals.hh"

G4bool G3CalcParamsFn(G4double* Rpar, G4int npar, G4double* Rparm,
                                  G4String shape, G4String shapem);

G4bool G3CalcParamsFn(G4double *Rpar, G4int npar, G4double *Rparm,
                       G4String shape, G4String shapem)
{
  G4bool NegPresent = FALSE;
  // for normal single volume positioning, just substitute for the
  // negative parameters
  // treat only the legal cases
  if (shapem == shape) {
    if (shape == "BOX" || shape == "TRD1" || shape == "TRD2") {
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
