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
// $Id: G3toG4BuildTree.hh 67982 2013-03-13 10:36:03Z gcosmo $
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
