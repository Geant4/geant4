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
// $Id: G3toG4MANY.hh,v 1.1 2001-11-08 16:07:59 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

#include "g4std/vector"
#include "G3VolTableEntry.hh"
#include "G4Transform3D.hh"

typedef G4std::vector<G3VolTableEntry*>  G3VolTableEntryVector;

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
