// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MPVEntry.hh,v 1.1 1999-01-07 16:09:41 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
////////////////////////////////////////////////////////////////////////
// G4MPVEntry Class Definition 
////////////////////////////////////////////////////////////////////////
//
// File:	G4MPVEntry.hh	
// Description: A G4MPVEntry is an MaterialPropertyVector Entry.  
//		One Material Property Vector contains many MPVEntries  
// Version:	1.0
// Created:	1996-02-08	
// Author:	Juliet Armstrong
// Updated:     1997-03-25 by Peter Gumplinger
//		> cosmetics (only)
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#ifndef G4MPVEntry_h
#define G4MPVEntry_h 1

/////////////
// Includes
/////////////

#include "G4ios.hh"
#include "globals.hh"

/////////////////////
// Class Definition
/////////////////////

class G4MPVEntry {

public:

	//////////////
	// Operators
	//////////////G4MPVEntry.hh
		
	// Well defined semantics for these operators
	// required by RWTPtrSortedVector

	G4bool operator <(const G4MPVEntry &right) const;	
	G4bool operator ==(const G4MPVEntry &right) const;
	G4MPVEntry& operator =(const G4MPVEntry &right);

	///////////////////////////////
	// Constructor and Destructor
	///////////////////////////////

	G4MPVEntry(G4double aPhotonMomentum, G4double aPropertyValue); 

	G4MPVEntry(const G4MPVEntry &right);

	~G4MPVEntry();

	////////////
	// Methods
	////////////

	G4double GetPhotonMomentum();

	G4double GetProperty();
	
	//////////
	// Tests
	//////////

	void DumpEntry();

private:

	/////////////////////////
	// Private Data members 
	/////////////////////////

	G4double thePhotonMomentum;
	G4double theProperty;
};

////////////////////
// Inline methods
////////////////////

// GetPhotonMomentum
// -----------------
//

inline 
G4double G4MPVEntry::GetPhotonMomentum()
{
	return thePhotonMomentum;
}

// GetProperty
// -----------
//

inline 
G4double G4MPVEntry::GetProperty()
{
	return theProperty;
}

#endif /* G4MPVEntry_h */
