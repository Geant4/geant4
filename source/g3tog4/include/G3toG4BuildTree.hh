// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3toG4BuildTree.hh,v 1.5 2000-11-24 09:50:10 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------
// Class Description:
//
// Definition of a global method:
//
//         void G3toG4BuildTree(G3VolTableEntry* curVTE,
//                              G3VolTableEntry* motherVTE)
//
// which processes the G3 volumes table (G3VolTable) and creates
// remaining G4 objects that could not be created during the phase of
// filling the G3 tables (defining G3 geometry, eg. by parsing the G3
// input via clparse.cc):
// G4LogicalVolume, G4PVPlacement and G4PVReplica objects.
// After processing of this method the G4 geometry is completely
// created and the G3 tables can be deleted.

// ----------------------
//
// modified by I.Hrivnacova, 13.10.99

#ifndef G3TOG4BUILDTREE_HH
#define G3TOG4BUILDTREE_HH 1

#include "G3VolTableEntry.hh"

void G3toG4BuildTree(G3VolTableEntry* curVTE, G3VolTableEntry* motherVTE);

#endif  
