//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4gsdvt.cc,v 1.6 2001-07-11 09:59:01 gunter Exp $
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

void PG4gsdvt(G4String tokens[])
{
  // fill the parameter containers
  G3fillParams(tokens,PTgsdvt);
  
  // interpret the parameters
  G4String vname = Spar[0];
  G4String vmoth = Spar[1];
  G4int iaxis = Ipar[0];
  G4int numed = Ipar[1];
  G4int ndvmx = Ipar[2];
  G4double Step = Rpar[0];
  
  G4gsdvt(vname,vmoth,Step,iaxis,numed,ndvmx);
}

void G4gsdvt(G4String vname, G4String vmoth, G4double step, G4int iaxis,
             G4int numed, G4int ndvmx)
{
  // find mother VTE
  G3VolTableEntry* mvte = G3Vol.GetVTE(vmoth);
  if (mvte == 0) {
    G4Exception("G4gsdvt:'" + vmoth + "' has no VolTableEntry");
  }    
  else {
    // a new vte clone copy with division is created
    // for each mother (clone copy)
    
    G4CreateCloneVTEWithDivision(vname, mvte, 
                                  kDvt, ndvmx, iaxis, numed, 0., step); 
  }  
}
