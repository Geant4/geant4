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
// $Id: G3toG4BuildTree.hh,v 1.7 2001-07-11 09:58:58 gunter Exp $
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
void G3toG4BuildLVTree(G3VolTableEntry* curVTE, G3VolTableEntry* motherVTE);
void G3toG4BuildPVTree(G3VolTableEntry* curVTE);

#endif  
