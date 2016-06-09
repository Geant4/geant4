//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G3NegVolPars.cc,v 1.13 2006-06-29 18:13:08 gunter Exp $ 
// GEANT4 tag $Name: not supported by cvs2svn $
//
// modified by I. Hrivnacova, 13.10.99 

#include "globals.hh"
#include "G3VolTable.hh"
#include "G4VSolid.hh"
#include "G3toG4.hh"
#include <cmath>

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
        if (Rparm != 0) Rpar[0] = std::min(Rparm[0],Rparm[1]) +
                                  std::abs(Rparm[0]-Rparm[1])*.5*Rpar[2]/Rparm[3];
	if (Rpar[0] < 0) NegPresent = TRUE;
      }
    }
    if (shapem == "TRD2") {
      if (Rpar[2] < 0) {
        if (Rparm != 0) Rpar[2] = Rparm[4];
        if (Rpar[2]<0) NegPresent = TRUE;
      }
      if (Rpar[0] < 0) {
        if (Rparm != 0) Rpar[0] = std::min(Rparm[0],Rparm[1]) +
                                  std::abs(Rparm[0]-Rparm[1])*.5*Rpar[2]/Rparm[4];
        if (Rpar[0]<0) NegPresent = TRUE;
      }
      if (Rpar[1] < 0) {
        if (Rparm != 0) Rpar[1] = std::min(Rparm[2],Rparm[3]) +
                                  std::abs(Rparm[2]-Rparm[3])*.5*Rpar[2]/Rparm[4];
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
          G4double xlo = std::min(Rparm[4],Rparm[8]) +
            std::abs(Rparm[4]-Rparm[8])*.5*Rpar[2]/Rparm[0];
          G4double xhi = std::min(Rparm[5],Rparm[9]) +
            std::abs(Rparm[5]-Rparm[9])*.5*Rpar[2]/Rparm[0];
          Rpar[0] = std::min(xlo,xhi);
        }
        if (Rpar[0] < 0) NegPresent = TRUE;
      }
      if (Rpar[1] < 0) {
        if (Rparm != 0) Rpar[1] = std::min(Rparm[3],Rparm[7]) +
                        std::abs(Rparm[3]-Rparm[7])*.5*Rpar[2]/Rparm[0];
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
