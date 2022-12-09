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
//
// by I.Hrivnacova, 13.10.99

#include "G3G4Interface.hh"
#include "G3toG4.hh"
#include "G3VolTable.hh"
#include "G3toG4MakeSolid.hh"
#include "G3Division.hh"
#include "G4SystemOfUnits.hh"
#include "G4VSolid.hh"

G4bool G3NegVolPars(G4double pars[], G4int *nparpt, 
          G3VolTableEntry* vte, G3VolTableEntry* mvte, const char routine[]);		       

void PG4gsposp(G4String *tokens){
  // fill the parameter containers
  G3fillParams(tokens,PTgsposp);
  
  // interpret the parameters
  G4String name = Spar[0];
  G4String moth = Spar[1];
  G4String only = Spar[2];
  G4int num = Ipar[0];
  G4int irot = Ipar[1];
  G4int npar = Ipar[2];
  // all parameters are passed to G4gsxxx methods
  // in G3 default units         
  //G4double x = Rpar[0]*cm;
  //G4double y = Rpar[1]*cm;
  //G4double z = Rpar[2]*cm;
  G4double x = Rpar[0];
  G4double y = Rpar[1];
  G4double z = Rpar[2];
  G4double *pars = &Rpar[3];
  
  G4gsposp(name, num, moth, x, y, z, irot, only, pars, npar);
}

void G4ProcessDaughters(G3VolTableEntry* vte)
// restore negative volume parameters and create solid for all
// vte daughters
{
  if (vte->HasNegPars()) {
    G4cerr << " Warning:" << G4endl;
    G4cerr << " G4ProcessDaughters: Ignored (vte has negative parameters)." 
           << G4endl;
  }  
  else {  
    for (G4int i=0; i<vte->GetNoDaughters(); i++) {
   
      G3VolTableEntry* dvte = vte->GetDaughter(i);

      if (dvte->HasNegPars()) {
        if (dvte->GetDivision()) {
           // call division method for creating solid and updating
           // dvte parameters
           dvte->GetDivision()->UpdateVTE();
        }
        else {
          // update negative parameters
          G4double* pars = dvte->GetRpar();
          G4int     npar = dvte->GetNpar();
          G4bool negpars 
	    = G3NegVolPars(pars,&npar, dvte, vte, "GSPOS");

          if (negpars) {
            G4String text = "G3NegVolPars still returns negative parameters!";
            G4Exception("G4ProcessDaughters()", "G3toG40019",
                        FatalException, text);
            return;
  	  }  

          // create solid
          G4bool hasNegPars;
          G4bool deferred;   
          G4bool okAxis[3];
          G4VSolid* solid
            = G3toG4MakeSolid(dvte->GetName(), dvte->GetShape(), pars, npar, 
	                      hasNegPars, deferred, okAxis);  
          if (hasNegPars) {
            G4String text = "G3toG4MakeSolid still returns negative parameters!";
            G4Exception("G4ProcessDaughters()", "G3toG40020",
                        FatalException, text);
            return;
	  }  

          // update dvte			  
          dvte->SetNRpar(npar, pars);
          dvte->SetSolid(solid);
          dvte->SetHasNegPars(hasNegPars);
        }

        // process daughters
	G4ProcessDaughters(dvte);  
      }		
    }
  }     
}
    
void G4CloneDaughters(G3VolTableEntry* vte, G3VolTableEntry* vteClone)
// copy vte daughters to vteClone
// (in case of daughters with negative parameters
// or with divisions new clone copies have to be created)
{
  G4int nofDaughters = vte->GetNoDaughters();
  if (nofDaughters>0)
    for (G4int id=0; id<nofDaughters; id++) {	  
      G3VolTableEntry* dvte = vte->GetDaughter(id);

      if (dvte->HasNegPars() || dvte->GetDivision()){
        // create new dvteClone with Position/Division
	// and set it to vteClone       

        // get master of dvte
        G3VolTableEntry* dvteMaster = dvte->GetMasterClone();

        // generate vteClone name
        G4int cloneNo = dvteMaster->GetNoClones();
        G4String newName = dvteMaster->GetName();
        newName += gSeparator;
        newName = newName + std::to_string(cloneNo);
        
        // create dvteClone
        G4String  dvteShape = dvte->GetShape();
	G4double* dvteRpar  = dvte->GetRpar();
	G4int     dvteNpar  = dvte->GetNpar();
	G4int     dvteNmed  = dvte->GetNmed();
	G4bool   hasNegPars = dvte->HasNegPars();
        G3VolTableEntry* dvteClone
          = new G3VolTableEntry(newName, dvteShape, dvteRpar, dvteNpar,
	                         dvteNmed, 0, hasNegPars);
				 
        // let dvte master and vol table know about it
        G3Vol.PutVTE(dvteClone);
        dvteMaster->AddClone(dvteClone);

        // set mother daughter
        vteClone->AddDaughter(dvteClone);
        dvteClone->AddMother(vteClone);
	
        // copy positions
        G4int nofPositions = dvte->NPCopies();
	for (G4int ip=0; ip<nofPositions; ip++)
	  dvteClone->AddG3Pos(dvte->GetG3PosCopy(ip));
	  
        // copy division
	G3Division* dvteDivision = dvte->GetDivision();
	if (dvteDivision) {          
	  G3Division* dvteCloneDivision
	    = new G3Division(dvteClone, vteClone, *dvteDivision);
          dvteClone->SetDivision(dvteCloneDivision);
          dvteCloneDivision->UpdateVTE();
	}  
                                
	// clone daughters recursively
	G4CloneDaughters(dvte, dvteClone);
      }	
      else {     
        // set dvte to vteClone
        vteClone->AddDaughter(dvte);
        dvte->AddMother(vteClone);
      }	
    }
}

void G4CreateCloneVTE(G3VolTableEntry* vte, G3VolTableEntry* mvte,
              G4double pars[], G4int npar, G4int num,
              G4double x, G4double y, G4double z, G4int irot, G4String vonly)
//
// create a new vte clone copy for each mother
// and derive its parameters from the mother if possible
{
   // create a G3Pos
   G4ThreeVector* offset = new G4ThreeVector(x*cm, y*cm, z*cm);
   G3Pos* aG3Pos = new G3Pos(mvte->GetName(), num, offset, irot, vonly);
            
   // loop over all mothers
   for (G4int i=0; i<mvte->GetNoClones(); i++) {
                    // mvte was retrieved from its "master" name
                    // -> there is no need to call GetMasterClone()
      G3VolTableEntry* mvteClone = mvte->GetClone(i);
        
      G4String tmpName = "TRY";
      G4String vteShape = vte->GetShape();
      G3VolTableEntry* vteClone
        = new G3VolTableEntry(tmpName, vteShape, pars, npar, vte->GetNmed(),
                               0, true);
          
      // negative parameters will be updated only 
      // for vteClone, pars are unchanged
      G4double* clonePars = vteClone->GetRpar();
      G4int     cloneNpar = vteClone->GetNpar();
      G4bool negpars
        = G3NegVolPars(clonePars, &cloneNpar, vteClone, mvteClone, "GSPOS");
      vteClone->SetHasNegPars(negpars);

      G3VolTableEntry* vteSameClone = 0;
      G4VSolid* solid = 0;
      if (!negpars) {
        // check if vteClone with the same parameters exist
        for (G4int ic=0; ic<vte->GetNoClones(); ic++) {
          G3VolTableEntry* checkClone = vte->GetClone(ic);
	  G4int checkNpar = checkClone->GetNpar();
	  G4double* checkPars = checkClone->GetRpar();
	  
          G4bool isSame;
          if (checkNpar != cloneNpar) 
	    isSame = false;
	  else {  
  	    isSame = true;
	    for (G4int ip=0; ip<cloneNpar; ip++) 
	      if (checkPars[ip] != clonePars[ip]) {
	        isSame = false;
		break;  
	      }
	  } 
	  if (isSame) { vteSameClone = checkClone; break; }
        }
       
        if (vteSameClone) {
          delete vteClone;
	  
          // add aG3Pos to the vteClone  
          vteSameClone->AddG3Pos(aG3Pos);
          mvteClone->AddDaughter(vteSameClone);  
          vteSameClone->AddMother(mvteClone);
        }
	else {
          // create the solid
          G4bool hasNegPars;
          G4bool deferred;
          G4bool okAxis[3];
          G4String vteName = vte->GetName();
          G4String cloneShape = vteClone->GetShape();
          solid = G3toG4MakeSolid(vteName, cloneShape, clonePars, cloneNpar,
                                  hasNegPars, deferred, okAxis);
        }				  
      }
   
      if ( negpars || !(vteSameClone)) {
        // generate vteClone name	  
        G4int cloneNo = vte->GetNoClones();
        G4String newName = vte->GetName();
        newName += gSeparator;
        newName = newName + std::to_string(cloneNo);
        
        // update vteClone
        vteClone->SetName(newName);
        vteClone->SetSolid(solid);
        vteClone->SetHasNegPars(negpars);
                                
        // let vte and vol table know about it
        G3Vol.PutVTE(vteClone);
        vte->AddClone(vteClone);
          
        // add aG3Pos to the vteClone  
        vteClone->AddG3Pos(aG3Pos);
        mvteClone->AddDaughter(vteClone);
        vteClone->AddMother(mvteClone);
           
        // copy all daughters
        G4CloneDaughters(vte, vteClone);

        // retrieve daughters parameters
        if (!negpars) G4ProcessDaughters(vteClone);
      }
  }
}

void G4gsposp(G4String vname, G4int num, G4String vmoth, G4double x,
              G4double y, G4double z, G4int irot, G4String vonly,
              G4double pars[], G4int npar)
{
  // find VTEs
  G3VolTableEntry* vte = G3Vol.GetVTE(vname);
  G3VolTableEntry* mvte = G3Vol.GetVTE(vmoth);

  if (vte == 0) {
    G4String err_message1 = "G4gsposp: '" + vname + "' has no VolTableEntry";
    G4Exception("G4psposp()", "G3toG40021", FatalException, err_message1);
    return;
  } 
  if (mvte == 0) {
    G4String err_message2 = "G4gsposp: '" + vmoth + "' has no VolTableEntry";
    G4Exception("G4psposp()", "G3toG40022", FatalException, err_message2);
    return;
  } 
  else { 
    // a new vte clone copy is created for each mother (clone copy)  
    // and its parameters are derived from it if possible
    
    G4CreateCloneVTE(vte, mvte, pars, npar, num, x, y, z, irot, vonly);
  }
}
