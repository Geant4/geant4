// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MaterialPropertiesTable.hh,v 1.1 1999-01-07 16:09:41 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
////////////////////////////////////////////////////////////////////////
// G4MaterialPropertiesTable Definition
////////////////////////////////////////////////////////////////////////
//
// File:        G4MaterialPropertiesTable.hh
// Description: An Material properties table is a hash table, with 
// 		key = property name, and value = G4MaterialPropertyVector 
// Version:     1.0
// Created:     1996-02-08
// Author:      Juliet Armstrong
// Updated:     1997-03-25 by Peter Gumplinger
//              > cosmetics (only)
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#ifndef G4MaterialPropertiesTable_h
#define G4MaterialPropertiesTable_h 1

/////////////
// Includes
/////////////

#include <math.h>
#include <rw/tphdict.h>
#include <rw/cstring.h>
#include "globals.hh"
#include "G4MaterialPropertyVector.hh"

#define RefractionIndex "RINDEX"

// unsigned hashString(const G4String &str);

/////////////////////
// Class Definition
/////////////////////

class G4MaterialPropertiesTable {

public:

	//////////////
        // Operators
        //////////////

	G4MaterialPropertiesTable&
		operator =(const G4MaterialPropertiesTable &right);
		
	////////////////
	// Constructor
	////////////////

	G4MaterialPropertiesTable(); 
		  
	G4MaterialPropertiesTable(const G4MaterialPropertiesTable &right);

	///////////////
	// Destructor
	///////////////

	~G4MaterialPropertiesTable();

	////////////
	// Methods
	////////////

	void AddProperty(char     *key,
		         G4double *PhotonMomenta,
		         G4double *PropertyValues,
		         G4int     NumEntries);

	void AddProperty(char *key, G4MaterialPropertyVector *opv);

	void RemoveProperty(char *key);

	G4MaterialPropertyVector* GetProperty(char *key); 

	void AddEntry(char *key, G4double aPhotonMomentum,
                                 G4double  aPropertyValue);

	void RemoveEntry(char *key, G4double  aPhotonMomentum);

	void DumpTable();	

private:

	/////////////////////////
	// Private Data members
	/////////////////////////

	RWTPtrHashDictionary<G4String, G4MaterialPropertyVector> MPT;  
	RWTPtrHashDictionaryIterator<G4String, G4MaterialPropertyVector> 
						MPTiterator; 
};

#endif /* G4MaterialPropertiesTable_h */
