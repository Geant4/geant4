// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gsdvt2.cc,v 1.5 1999-12-05 17:50:13 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by I.Hrivnacova, V.Berejnoi, 29 Oct 99

#include "G3Division.hh"
#include "G3VolTableEntry.hh"
#include "G3VolTable.hh"
#include "globals.hh"
#include "G3toG4.hh"

void G4CreateCloneVTEWithDivision(G4String vname, G3VolTableEntry* mvte,
               G3DivType divType, G4int nofDivisions, G4int iaxis, G4int nmed, 
     	       G4double c0, G4double step);

void PG4gsdvt2(G4String tokens[])
{
  // fill the parameter containers
  G3fillParams(tokens,PTgsdvt2);
  
  // interpret the parameters
  G4String vname = Spar[0];
  G4String vmoth = Spar[1];
  G4int iaxis = Ipar[0];
  G4int numed = Ipar[1];
  G4int ndvmx = Ipar[2];
  G4double Step = Rpar[0];
  G4double c0 = Rpar[1];
  
  G4gsdvt2(vname,vmoth,Step,iaxis,c0,numed,ndvmx);
}

void G4gsdvt2(G4String vname, G4String vmoth, G4double step, G4int iaxis,
               G4double c0, G4int numed, G4int ndvmx)
{
  // find mother VTE
  G3VolTableEntry* mvte = G3Vol.GetVTE(vmoth);
  if (mvte == 0) {
    G4Exception("G4gsdtn2:'" + vmoth + "' has no VolTableEntry");
  }    
  else {
    // a new vte clone copy with division is created
    // for each mother (clone copy)
    
    G4CreateCloneVTEWithDivision(vname, mvte, 
                                  kDvt2, ndvmx, iaxis, numed, c0, step); 
  }  
}
