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
// $Id$
//

// testG4SmartVoxelNodeProxy
//             Ensure asserts are compiled in
//
// Test file for G4SmartVoxelProxy
//
// o Simplistic checks for
//
//   IsNode
//   IsHeader
//   GetNode
//   GetHeader

#include <assert.h>

// Global defs
#include "globals.hh"

// Tested entities
#include "G4SmartVoxelProxy.hh"

// Required
#include "G4SmartVoxelNode.hh"
#include "G4SmartVoxelHeader.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"

G4bool testG4SmartVoxelNodeProxy()
{
    G4SmartVoxelNode *tNode;
    G4SmartVoxelProxy *tNodeProxy;
    tNode=new G4SmartVoxelNode(1);
    assert(tNode != 0);		// Sanity check
    tNodeProxy=new G4SmartVoxelProxy(tNode);
    assert(tNodeProxy != 0);		// Sanity check
   
    assert(tNodeProxy->IsNode());
    assert(!tNodeProxy->IsHeader());
    assert(tNodeProxy->GetNode()==tNode);
    delete tNodeProxy;
    delete tNode;
    return true;
}

G4bool testG4SmartVoxelHeaderProxy()
{
    G4Box tBox("dummyBox",1,1,1);
    G4LogicalVolume tVol(&tBox,0,"dummyLogical",0,0,0);

    G4SmartVoxelHeader *tHeader;
    G4SmartVoxelProxy *tHeaderProxy;
    tHeader=new G4SmartVoxelHeader(&tVol);
    assert(tHeader != 0);		// Sanity check
    tHeaderProxy=new G4SmartVoxelProxy(tHeader);
    assert(tHeaderProxy != 0);	// Sanity check
    
    assert(!tHeaderProxy->IsNode());
    assert(tHeaderProxy->IsHeader());
    assert(tHeaderProxy->GetHeader()==tHeader);
    delete tHeaderProxy;
    delete tHeader;
    return true;
}

int main()
{
#ifdef NDEBUG
    G4Exception("main","000",FatalException,"FAIL: *** Assertions must be compiled in! ***");
#endif
    assert(testG4SmartVoxelNodeProxy());
    assert(testG4SmartVoxelHeaderProxy());
    return 0;
}


