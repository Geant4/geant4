// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3CalcParams.cc,v 1.1 1999-01-07 16:06:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"
#include "G3toG4.hh"
#include "G3CalcParams.hh"
#include <math.h>

// $$$ finish and check

G4bool G3CalcParamsFn(G4double* Rpar, G4int npar, G4double* Rparm,
                                  G4String shape, G4String shapem);

void G3CalcParams::ComputeTransformation(const G4int n, // replication number
                                         G4VPhysicalVolume *pRep) const
{
}

void G3CalcParams::ComputeDimensions(G4Box &pBox,
                                     const G4int, // replication number
                                     const G4VPhysicalVolume *pRep) const
{
  
    G4LogicalVolume* lvol = pRep->GetLogicalVolume();
    G4VPhysicalVolume* pmoth = pRep->GetMother();
    G4LogicalVolume* lmoth = pmoth->GetLogicalVolume();
    G4VSolid* msolid = lmoth->GetSolid();
    pBox.SetXHalfLength(1.*cm); // $$$
    pBox.SetYHalfLength(1.*cm); // $$$
    pBox.SetZHalfLength(1.*cm); // $$$
}

void G3CalcParams::ComputeDimensions(G4Tubs &,
                                     const G4int, // replication number
                                     const G4VPhysicalVolume *) const
{
}

void G3CalcParams::ComputeDimensions(G4Trd &,
                                     const G4int, // replication number
                                     const G4VPhysicalVolume *) const
{
}

void G3CalcParams::ComputeDimensions(G4Trap &,
                                     const G4int, // replication number
                                     const G4VPhysicalVolume *) const
{
}

void G3CalcParams::ComputeDimensions(G4Cons &,
                                     const G4int, // replication number
                                     const G4VPhysicalVolume *) const
{
}

void G3CalcParams::ComputeDimensions(G4Sphere &,
                                     const G4int, // replication number
                                     const G4VPhysicalVolume *) const
{
}

void G3CalcParams::ComputeDimensions(G4Torus &,
                                     const G4int, // replication number
                                     const G4VPhysicalVolume *) const
{
}

void G3CalcParams::ComputeDimensions(G4Para &,
                                     const G4int, // replication number
                                     const G4VPhysicalVolume *) const
{
}

void G3CalcParams::ComputeDimensions(G4Hype &,
                                     const G4int, // replication number
                                     const G4VPhysicalVolume *) const
{
}

G4double *G3CalcParamsFunction(G4double *Rpar, G4int npar, G4double *Rparm,
                               G4String shape, G4String shapem)
{
// will have available the parameters of the mother (and the daughter).
// what should really expect is pointer to mother solid, so extract
// shape and params, and daughter solid, extract shape and params,
// and then return updated params with a G4double* return.
// In absence of methods to obtain params, shapes from solids, assume
// I get what I need as call params.

    G4bool NegPresent = G3CalcParamsFn(Rpar, npar, Rparm, shape, shapem);
    return Rpar;
}

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
          if (Rparm != NULL) Rpar[i] = Rparm[i];
          NegPresent = TRUE;
        }
      }
    }
    if (shape == "TRAP") {
      for (G4int i=0;i<11;i++) {
        if (i != 1 && i != 2 && i != 6 && i != 10) {
          if (Rpar[i]<0) {
            if (Rparm != NULL) Rpar[i] = Rparm[i];
            NegPresent = TRUE;
          }
        }
      }
    }
    if (shape == "TUBE" || shape == "TUBS" || shape == "PARA") {
      for (G4int i=0;i<3;i++) {
        if (Rpar[i] < 0) {
          if (Rparm != NULL) Rpar[i] = Rparm[i];
          NegPresent = TRUE;
        }
      }
    }
    if (shape == "CONE" || shape == "CONS") {
      for (G4int i=0;i<5;i++) {
        if (Rpar[i] < 0) {
          if (Rparm != NULL) Rpar[i] = Rparm[i];
          NegPresent = TRUE;
        }
      }
    }
    if (shape == "SPHE") {
      for (G4int i=0;i<2;i++) {
        if (Rpar[i] < 0) {
          if (Rparm != NULL) Rpar[i] = Rparm[i];
          NegPresent = TRUE;
        }
      }
    }
    if (shape == "PGON") {
      G4int nz = int(Rpar[3]);
      G4int ipl;
      for (G4int i=0;i<nz;i++) {
        ipl = 5 + i*3;
        if (Rpar[ipl] < 0) {
          if (Rparm != NULL) Rpar[ipl] = Rparm[ipl];
          NegPresent = TRUE;
        }
        if (Rpar[ipl+1] < 0) {
          if (Rparm != NULL) Rpar[ipl] = Rparm[ipl];
          NegPresent = TRUE;
        }
      }
    }
    if (shape == "PCON") {
      G4int nz = int(Rpar[2]);
      G4int ipl;
      for (G4int i=0;i<nz;i++) {
        ipl = 4 + i*3;
        if (Rpar[ipl] < 0) {
          if (Rparm != NULL) Rpar[ipl] = Rparm[ipl];
          NegPresent = TRUE;
        }
        if (Rpar[ipl+1] < 0) {
          if (Rparm != NULL) Rpar[ipl] = Rparm[ipl];
          NegPresent = TRUE;
        }
      }
    }
  }

  if (shape == "BOX") {
    if (shapem == "TRD1") {
      if (Rpar[1] < 0) {
        if (Rparm != NULL) Rpar[1] = Rparm[2];
        NegPresent = TRUE;
      }
      if (Rpar[2] < 0) {
        if (Rparm != NULL) Rpar[2] = Rparm[3];
        NegPresent = TRUE;
      }
      if (Rpar[0] < 0) {
        if (Rparm != NULL) Rpar[0] = min(Rparm[0],Rparm[1]) +
          abs(Rparm[0]-Rparm[1])*.5*Rpar[2]/Rparm[3];
        NegPresent = TRUE;
      }
    }
    if (shapem == "TRD2") {
      if (Rpar[2] < 0) {
        if (Rparm != NULL) Rpar[2] = Rparm[4];
        NegPresent = TRUE;
      }
      if (Rpar[0] < 0) {
        if (Rparm != NULL) Rpar[0] = min(Rparm[0],Rparm[1]) +
          abs(Rparm[0]-Rparm[1])*.5*Rpar[2]/Rparm[4];
        NegPresent = TRUE;
      }
      if (Rpar[1] < 0) {
        if (Rparm != NULL) Rpar[1] = min(Rparm[2],Rparm[3]) +
          abs(Rparm[2]-Rparm[3])*.5*Rpar[2]/Rparm[4];
        NegPresent = TRUE;
      }
    }
    if (shapem == "TRAP") {
      if (Rpar[2] < 0) {
        if (Rparm != NULL) Rpar[2] = Rparm[0];
        NegPresent = TRUE;
      }
      if (Rpar[0] < 0) {
        if (Rparm != NULL) {
          G4double xlo = min(Rparm[4],Rparm[8]) +
            abs(Rparm[4]-Rparm[8])*.5*Rpar[2]/Rparm[0];
          G4double xhi = min(Rparm[5],Rparm[9]) +
            abs(Rparm[5]-Rparm[9])*.5*Rpar[2]/Rparm[0];
          Rpar[0] = min(xlo,xhi);
        }
        NegPresent = TRUE;
      }
      if (Rpar[1] < 0) {
        if (Rparm != NULL) Rpar[1] = min(Rparm[3],Rparm[7]) +
          abs(Rparm[3]-Rparm[7])*.5*Rpar[2]/Rparm[0];
        NegPresent = TRUE;
      }
    }
  }
  return NegPresent;
}
