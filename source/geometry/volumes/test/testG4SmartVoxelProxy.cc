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
// $Id: testG4SmartVoxelProxy.cc,v 1.3 2001-07-11 10:00:36 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif
    assert(testG4SmartVoxelNodeProxy());
    assert(testG4SmartVoxelHeaderProxy());
    return 0;
}


