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

// testBlockingList 
//             Ensure asserts are compiled in
//
// History:
//
// 24.7.96 P.Kent Verify all functions

#include <assert.h>

#include "globals.hh"
#include "G4BlockingList.hh"

G4bool testG4BlockingList()
{
	G4BlockingList bList(10,3);	// Size 10, resize 3
	assert (bList.Length()==10);
	assert (bList.IsBlocked(0)==false);
	assert (bList.IsBlocked(1)==false);
	bList.BlockVolume(0);
	assert (bList.IsBlocked(0)==true);
	assert (bList.IsBlocked(1)==false);
	bList.BlockVolume(1);
	assert (bList.IsBlocked(0)==true);
	assert (bList.IsBlocked(1)==true);

	bList.Reset();
	assert (bList.Length()==10);
	assert (bList.IsBlocked(0)==false);
	assert (bList.IsBlocked(1)==false);
	
	bList.Enlarge(16);
	assert (bList.Length()==18);

	bList.BlockVolume(15);
	assert (bList.IsBlocked(15)==true);
	bList.Reset();
	assert (bList.IsBlocked(15)==false);

	return true;
}

int main()
{
#ifdef NDEBUG
  G4Exception("main","000",FatalException,"FAIL: *** Assertions must be compiled in! ***");
#endif
    assert(testG4BlockingList());
    return 0;
}

