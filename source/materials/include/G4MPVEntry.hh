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
// $Id: G4MPVEntry.hh,v 1.6 2006/06/29 19:11:11 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
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
