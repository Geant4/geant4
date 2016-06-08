// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3toG4MakeSolid.hh,v 1.3 1999/12/05 17:50:04 gcosmo Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
#ifndef _G3TOG4MAKESOLID_
#define _G3TOG4MAKESOLID_ 1
G4VSolid* G3toG4MakeSolid(const G4String& vname, const G4String& shape, 
			  const G4double* Rpar, const G4int npar, 
			  G4bool& NegVolPars, G4bool& Deferred, 
			  G4bool* OKAxis);
#endif

    
