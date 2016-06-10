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
// $Id: G3toG4MANY.hh 67982 2013-03-13 10:36:03Z gcosmo $
//
// ----------------------
// Class Description:
//
// Definition of global methods:
//
// void G3toG4MANY(G3VolTableEntry* curVTE);
//
// void MakeBooleanSolids(G3VolTableEntry* curVTE, 
//                        G3VolTableEntryVector* overlaps,
//                        const G4Transform3D& transform);
//
// void SubstractSolids(G3VolTableEntry* vte1, 
//                      G3VolTableEntry* vte2,
//                      G4int copy, 
//                      const G4Transform3D& transform);
//
// G4Transform3D GetTransform3D(G3Pos*);
//
// MANY positions are resolved in G3toG4MANY() function, which has to be
// processed before G3toG4BuildTree() (it is not called by default).
// In order to resolve MANY, the user code has to provide additional info
// using G4gsbool(G4String volName, G4String manyVolName) function, for
// all overlapping volumes. Daughters of overlapping volumes are then
// resolved automatically and should not be specified via Gsbool.
//
// Limitation: a volume with a MANY position can have only this one
// ----------  position; if more than one position is needed a new volume
//             has to be defined (gsvolu) for each position.

// ----------------------
//
// By I.Hrivnacova, 22.10.01

#ifndef G3TOG4_MANY_HH
#define G3TOG4_MANY_HH 1

#include <vector>
#include "G3VolTableEntry.hh"
#include "G4Transform3D.hh"

typedef std::vector<G3VolTableEntry*>  G3VolTableEntryVector;

void G3toG4MANY(G3VolTableEntry* curVTE);

void MakeBooleanSolids(G3VolTableEntry* curVTE, 
                       G3VolTableEntryVector* overlaps,
		       const G4Transform3D& transform);

void SubstractSolids(G3VolTableEntry* vte1, 
                     G3VolTableEntry* vte2,
	             G4int copy, 
		     const G4Transform3D& transform);

G4Transform3D GetTransform3D(G3Pos*);

#endif  
