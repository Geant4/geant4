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
// $Id: G3NegVolPars.cc 67982 2013-03-13 10:36:03Z gcosmo $ 
//
// modified by I. Hrivnacova, 13.10.99 

#include "globals.hh"
#include "G3VolTable.hh"
#include "G4VSolid.hh"
#include "G3toG4.hh"
#include <cmath>

G4bool G3CalcParamsFn(G4double *rpar, G4int npar, G4double *rparm,
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
        if (rpar[i] < 0) {
          if (rparm != 0) rpar[i] = rparm[i];
	  if (rpar[i] < 0) NegPresent = TRUE;
        }
      }
    }
    if (shape == "TRAP") {
      for (G4int i=0;i<11;i++) {
        if (i != 1 && i != 2 && i != 6 && i != 10) {
          if (rpar[i]<0) {
            if (rparm != 0) rpar[i] = rparm[i];
	    if (rpar[i] < 0) NegPresent = TRUE;
	  }  
        }
      }
    }
    if (shape == "TUBE" || shape == "TUBS" || shape == "PARA") {
      for (G4int i=0;i<3;i++) {
        if (rpar[i] < 0) {
          if (rparm != 0) rpar[i] = rparm[i];
	  if (rpar[i] < 0) NegPresent = TRUE;
        }
      }
    }
    if (shape == "CONE" || shape == "CONS") {
      for (G4int i=0;i<5;i++) {
        if (rpar[i] < 0) {
          if (rparm != 0) rpar[i] = rparm[i];
	  if (rpar[i] < 0) NegPresent = TRUE;
        }
      }
    }
    if (shape == "SPHE") {
      for (G4int i=0;i<2;i++) {
        if (rpar[i] < 0) {
          if (rparm != 0) rpar[i] = rparm[i];
	  if (rpar[i] < 0) NegPresent = TRUE;
        }
      }
    }
    if (shape == "PGON") {
      G4int nz = G4int(rpar[3]);
      G4int ipl;
      for (G4int i=0;i<nz;i++) {
        ipl = 5 + i*3;
        if (rpar[ipl] < 0) {
          if (rparm != 0) rpar[ipl] = rparm[ipl];
	  if (rpar[ipl] < 0) NegPresent = TRUE;
        }
        if (rpar[ipl+1] < 0) {
          if (rparm != 0)  rpar[ipl] = rparm[ipl];
	  if (rpar[ipl] < 0) NegPresent = TRUE;
        }
      }
    }
    if (shape == "PCON") {
      G4int nz = G4int(rpar[2]);
      G4int ipl;
      for (G4int i=0;i<nz;i++) {
        ipl = 4 + i*3;
        if (rpar[ipl] < 0) {
          if (rparm != 0) rpar[ipl] = rparm[ipl];
	  if (rpar[ipl] < 0) NegPresent = TRUE;
        }
        if (rpar[ipl+1] < 0) {
          // TO DO
	  // check - folowing argument might be ipl+1
          if (rparm != 0) rpar[ipl] = rparm[ipl];
	  if (rpar[ipl] < 0) NegPresent = TRUE;
        }
      }
    }
  }

  if (shape == "BOX") {
    if (shapem == "TRD1") {
      if (rpar[1] < 0) {
        if (rparm != 0) rpar[1] = rparm[2];
	if (rpar[1] < 0) NegPresent = TRUE;
      }
      if (rpar[2] < 0) {
        if (rparm != 0) rpar[2] = rparm[3];
	if (rpar[2] < 0) NegPresent = TRUE;
      }
      if (rpar[0] < 0) {
        if (rparm != 0) rpar[0] = std::min(rparm[0],rparm[1]) +
                                  std::abs(rparm[0]-rparm[1])*.5*rpar[2]/rparm[3];
	if (rpar[0] < 0) NegPresent = TRUE;
      }
    }
    if (shapem == "TRD2") {
      if (rpar[2] < 0) {
        if (rparm != 0) rpar[2] = rparm[4];
        if (rpar[2]<0) NegPresent = TRUE;
      }
      if (rpar[0] < 0) {
        if (rparm != 0) rpar[0] = std::min(rparm[0],rparm[1]) +
                                  std::abs(rparm[0]-rparm[1])*.5*rpar[2]/rparm[4];
        if (rpar[0]<0) NegPresent = TRUE;
      }
      if (rpar[1] < 0) {
        if (rparm != 0) rpar[1] = std::min(rparm[2],rparm[3]) +
                                  std::abs(rparm[2]-rparm[3])*.5*rpar[2]/rparm[4];
        if (rpar[1]<0) NegPresent = TRUE;
      }
    }
    if (shapem == "TRAP") {
      if (rpar[2] < 0) {
        if (rparm != 0) rpar[2] = rparm[0];
        if (rpar[2] < 0) NegPresent = TRUE;
      }
      if (rpar[0] < 0) {
        if (rparm != 0) {
          G4double xlo = std::min(rparm[4],rparm[8]) +
            std::abs(rparm[4]-rparm[8])*.5*rpar[2]/rparm[0];
          G4double xhi = std::min(rparm[5],rparm[9]) +
            std::abs(rparm[5]-rparm[9])*.5*rpar[2]/rparm[0];
          rpar[0] = std::min(xlo,xhi);
        }
        if (rpar[0] < 0) NegPresent = TRUE;
      }
      if (rpar[1] < 0) {
        if (rparm != 0) rpar[1] = std::min(rparm[3],rparm[7]) +
                        std::abs(rparm[3]-rparm[7])*.5*rpar[2]/rparm[0];
        if (rpar[1] < 0) NegPresent = TRUE;
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
  G4double* rpar = vte->GetRpar();
  G4int npar = vte->GetNpar();  
  if (npar ==0) {
    // no solid parameters are defined in vte
    npar = *nparpt;
    rpar = pars;
  }
  else {  
    // solid parameters are already defined in vte
    // pars[], nparpt are ignored
    // TO DO: check if g3 ignores them too or resets
    // vte parameters according to this new ones !!
  }  
      
  // mother
  G4String shapem = mvte->GetShape();
  G4double* rparm = mvte->GetRpar();

  if (strcmp(routine,"GSPOS") == 0 || strcmp(routine,"GSVOLU") == 0) {
    NegPresent = G3CalcParamsFn(rpar,npar,rparm,shape,shapem);
  }
  if (strcmp(routine,"GSDVN") == 0) {
    // just set the flag. The parametrization function figures out
    // what to do.
    for (G4int i=0;i<npar;i++) {
      if (rpar[i] < 0) {
        NegPresent = TRUE;
      }
    }
  }
  return NegPresent;
}
