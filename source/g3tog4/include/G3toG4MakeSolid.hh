// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3toG4MakeSolid.hh,v 1.5 2000-11-24 09:50:11 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------
// Class Description:
//
// Definition of a global method:
//
//   G4VSolid* G3toG4MakeSolid(const G4String& vname,
//                             const G4String& shape, 
//                             const G4double* Rpar,
//                             const G4int npar, 
//                                   G4bool& NegVolPars,
//                                   G4bool& Deferred, 
//                                   G4bool* OKAxis);
//
// which checks the volume parameters and creates the G4VSolid
// subclass object corresponding to the specified shape. 
// If volume parameters are incomplete (negative or none)
// it returns 0.

// ----------------------

#ifndef G3TOG4MAKESOLID_HH
#define G3TOG4MAKESOLID_HH 1

G4VSolid* G3toG4MakeSolid(const G4String& vname, const G4String& shape, 
			  const G4double* Rpar, const G4int npar, 
			  G4bool& NegVolPars, G4bool& Deferred, 
			  G4bool* OKAxis);
#endif

    
