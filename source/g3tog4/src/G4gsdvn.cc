// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gsdvn.cc,v 1.5 1999-12-05 17:50:13 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// by I.Hrivnacova, V.Berejnoi, 29 Oct 99

#include "G3G4Interface.hh"
#include "G3Division.hh"
#include "G3VolTableEntry.hh"
#include "G3VolTable.hh"
#include "globals.hh"
#include "G3toG4.hh"

void PG4gsdvn(G4String tokens[])
{
  // fill the parameter containers
  G3fillParams(tokens,PTgsdvn);
  
  // interpret the parameters
  G4String vname = Spar[0];
  G4String vmoth = Spar[1];
  G4int ndiv = Ipar[0];
  G4int iaxis = Ipar[1];
  
  G4gsdvn(vname, vmoth, ndiv, iaxis);
}

void G4CreateCloneVTEWithDivision(G4String vname, G3VolTableEntry* mvte,
        G3DivType divType, G4int nofDivisions, G4int iaxis, G4int nmed, 
	G4double c0, G4double step)
{	
  G3VolTableEntry* vte;

  // loop over all mothers
  for (G4int i=0; i<mvte->GetNoClones(); i++) {
    G3VolTableEntry* mvteClone =  mvte->GetClone(i);
    G4String shape = mvteClone->GetShape();
    G4int    nmed  = mvteClone->GetNmed();
    G4String mvteName = mvteClone->GetName();

    G4String newName = vname;
    if (i>0) {
      char index[4]; sprintf(index, "%d", i);
      newName.append(gSeparator); newName = newName + index;
    }	

    // create new VTE with 0 solid
    // and let vol table know about it
    G3VolTableEntry* vteClone
      = new G3VolTableEntry(newName, shape, 0, 0, nmed, 0, true);
    G3Vol.PutVTE(vteClone);

    // set mother <-> daughter
    // (mother/daughter are reset in case an envelope
    //  needs to be created in G3Division::UpdateVTE)
    mvteClone->AddDaughter(vteClone);
    vteClone->AddMother(mvteClone);

    // create new G3Division 
    G3Division* division
      = new G3Division(divType, vteClone, mvteClone,
                       nofDivisions, iaxis, nmed, c0, step);

    // set division to vte and update it
    vteClone->SetDivision(division);
    division->UpdateVTE();
      
    if (i == 0) {
      // keep the first clone copy
      vte = vteClone;
    }
    else {		
      // let vte know about this clone copy
      vte->AddClone(vteClone);     
    }
  }
}    

void G4gsdvn(G4String vname, G4String vmoth, G4int ndiv, G4int iaxis)
{
  // find mother VTE
  G3VolTableEntry* mvte = G3Vol.GetVTE(vmoth);
 
  if (mvte == 0) {
    G4Exception("G4gsdvn:'" + vmoth + "' has no VolTableEntry");
  }  
  else {
    // a new vte clone copy with division is created
    // for each mother (clone copy)
    
    G4CreateCloneVTEWithDivision(vname, mvte, kDvn, ndiv, iaxis, 0, 0., 0.); 
  }  
}
