// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4gsposp.cc,v 1.8 1999-05-22 06:51:00 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G3toG4.hh"
#include "G3VolTable.hh"
#include "G4makevol.hh"

void PG4gsposp(RWCString tokens[]){
  // fill the parameter containers
  G3fillParams(tokens,PTgsposp);
  
  // interpret the parameters
  G4String name = Spar[0];
  G4String moth = Spar[1];
  G4String only = Spar[2];
  G4int num = Ipar[0];
  G4int irot = Ipar[1];
  G4int npar = Ipar[2];
  G4double x = Rpar[0]*cm;
  G4double y = Rpar[1]*cm;
  G4double z = Rpar[2]*cm;
  G4double *pars = &Rpar[3];
  
  G4gsposp(name, num, moth, x, y, z, irot, only, pars, npar);
}

void G4gsposp(const G4String& vname, const G4int num, const G4String& vmoth, 
	      G4double x, G4double y, G4double z, const G4int irot, 
	      const G4String& vonly, G4double* pars, G4int npar){
  
  // VTable entry
  VolTableEntry* VTE = G3Vol.GetVTE(vname);

  if (VTE != 0) {
    // get some parameters from the current entry
    G4String Shape = VTE->GetShape();
    G4int Nmed = VTE->GetNmed();

    // get pointer to key for current entry
    const G4String* Key = G3Vol.GetVTD()->find(&vname);

    // remove the current entry from the dictionary
    G3Vol.GetVTD()->remove(Key);

    // insert a new entry
    G4makevol(vname, Shape, Nmed, pars, npar);

    G4gspos(vname, num, vmoth, x, y, z, irot, vonly);

  } else {
    G4cerr << "G4gsposp: no G3VolTable entry for logical volume '"
	   << vname << "'" << endl;
  }
}






