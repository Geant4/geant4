// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testG4BlockingList.cc,v 1.1 1999-01-08 16:31:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif
    assert(testG4BlockingList());
    return 0;
}

