// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4makevol.cc,v 1.13 1999-05-22 07:30:33 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "globals.hh"
#include "G4Material.hh"
#include "G4makevol.hh"
#include "G3VolTable.hh"
#include "G3MedTable.hh"
#include "G3toG4MakeSolid.hh"
        
void G4makevol(const G4String& vname, const G4String& shape, const G4int nmed,
	       const G4double* Rpar, const G4int npar){
    
  G4bool NegVolPars = false;
  G4bool Deferred   = false;
  G4VSolid* solid = G3toG4MakeSolid(vname, shape, Rpar, npar, NegVolPars,
				    Deferred);

  // get the material corresponding to the tracking medium

  G4Material* material=0;
  if (nmed>0) {
    material = G3Med.get(nmed);
  } else {
    G4cerr << "No material found for tracking medium index " << nmed << endl;
    assert(material != 0);
  }

  // external store is needed to handle deferred LV creation 
  G3Vol.PutVTE(vname, shape, Rpar, npar, nmed, material, solid, Deferred, 
	      NegVolPars);
}












































