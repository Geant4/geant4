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
// $Id: G3toG4MakeSolid.hh 67982 2013-03-13 10:36:03Z gcosmo $
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

    
