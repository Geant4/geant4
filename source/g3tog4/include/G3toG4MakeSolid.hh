// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3toG4MakeSolid.hh,v 1.4 1999-12-09 01:27:47 lockman Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G3TOG4MAKESOLID_HH
#define G3TOG4MAKESOLID_HH 1

G4VSolid* G3toG4MakeSolid(const G4String& vname, const G4String& shape, 
			  const G4double* Rpar, const G4int npar, 
			  G4bool& NegVolPars, G4bool& Deferred, 
			  G4bool* OKAxis);
#endif

    
