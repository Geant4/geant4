// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4makevol.cc,v 1.15 1999-05-28 02:19:38 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "globals.hh"
#include "G4Material.hh"
#include "G4makevol.hh"
#include "G3VolTable.hh"
#include "G3MedTable.hh"
#include "G3toG4MakeSolid.hh"
        
void G4makevol(G4String& vname, G4String& shape, G4int nmed,
	       G4double* Rpar, G4int npar){
  G4bool OKAxis[3];
  G4bool NegVolPars;
  G4bool Deferred;
  G4VSolid* solid = G3toG4MakeSolid(vname, shape, Rpar, npar, NegVolPars,
				    Deferred, OKAxis);

  // get the material corresponding to the tracking medium

  G4Material* material=0;
  if (nmed>0) {
    material = G3Med.get(nmed);
  } else {
    G4cerr << "No material found for tracking medium index " << nmed << endl;
    assert(material != 0);
  }

  if (G3Vol.GetVTE(vname) == 0) {
    // external store is needed to handle deferred LV creation 
    VolTableEntry* VTE = 
      new VolTableEntry(vname, shape, Rpar, npar, nmed, material, solid,
			Deferred, NegVolPars, OKAxis);
    G3Vol.PutVTE(VTE);
  }
  //  G4cout << "Listing VTD" << endl;
  //  G3Vol.ListVTE();
}












































