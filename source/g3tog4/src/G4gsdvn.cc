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
// $Id: G4gsdvn.cc 104024 2017-05-08 14:44:25Z gcosmo $
//
// by I.Hrivnacova, V.Berejnoi, 29 Oct 99

#include "G3G4Interface.hh"
#include "G3Division.hh"
#include "G3VolTableEntry.hh"
#include "G3VolTable.hh"
#include "globals.hh"
#include "G3toG4.hh"

void PG4gsdvn(G4String *tokens)
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
        G3DivType divType, G4int nofDivisions, G4int iaxis, G4int, 
	G4double c0, G4double step)
{	
  G3VolTableEntry* vte=0;

  // loop over all mothers
  for (G4int i=0; i<mvte->GetNoClones(); i++) {
    G3VolTableEntry* mvteClone =  mvte->GetClone(i);
    G4String shape = mvteClone->GetShape();
    G4int    nmed  = mvteClone->GetNmed();
    G4String mvteName = mvteClone->GetName();
    
    G4String newName = vname;
    if (i>0) {
      char index[12]; sprintf(index, "%d", i);
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
    G4String text = "G4gsdvn:'" + vmoth + "' has no VolTableEntry";
    G4Exception("G4gsdvn()", "G3toG40013", FatalException, text);
    return;
  }  
  else {
    // a new vte clone copy with division is created
    // for each mother (clone copy)
    
    G4CreateCloneVTEWithDivision(vname, mvte, kDvn, ndiv, iaxis, 0, 0., 0.); 
  }  
}
