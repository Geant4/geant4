// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BlockingList.cc,v 1.2 1999-12-15 14:50:26 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4BlockingList Implementation
//

#include "G4BlockingList.hh"

// Clear List and reset tag
// Fix: Out of line for HP-CC
void G4BlockingList::FullyReset()
{
	fBlockTagNo=1;
	for (G4int i=fBlockingList.length()-1;i>=0;i--)
		{
			fBlockingList(i)=0;
		}	
	
}
