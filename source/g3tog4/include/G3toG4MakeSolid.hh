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
// $Id: G3toG4MakeSolid.hh,v 1.6 2001-07-11 09:58:58 gunter Exp $
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

    
