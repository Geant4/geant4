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
// $Id: G4MPVEntry.hh,v 1.5 2001-07-11 10:01:26 gunter Exp $
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
// Updated:     1999-10-29 add method and class descriptors
//              1997-03-25 by Peter Gumplinger
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

// Class Description:
// A G4MPVEntry is an MaterialPropertyVector Entry.
// One Material Property Vector contains many MPVEntries.
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////

class G4MPVEntry {

public: // Without description

	//////////////
	// Operators
	//////////////G4MPVEntry.hh
		
	// Well defined semantics for these operators
	// required by G4RWTPtrSortedVector

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
