// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MaterialPropertiesTable.hh,v 1.2 1999-10-30 01:42:49 gum Exp $
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
// Updated:     1999-10-29 add method and class decriptors
//              1997-03-25 by Peter Gumplinger
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

// Class Description:
// An Material properties table is a hash table, with key = property 
// name, and value = G4MaterialPropertyVector.
// Class Description - End:

#define RefractionIndex "RINDEX"

// unsigned hashString(const G4String &str);

/////////////////////
// Class Definition
/////////////////////

class G4MaterialPropertiesTable {

public: // Without description

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

public: // With description

	void AddProperty(char     *key,
		         G4double *PhotonMomenta,
		         G4double *PropertyValues,
		         G4int     NumEntries);
        // Add a new property to the table by giving a key-name and the 
        // arrays x and y of size NumEntries.

	void AddProperty(char *key, G4MaterialPropertyVector *opv);
        // Add a new property to the table by giving a key-name and an 
        // already constructed G4MaterialPropertyVector.

	void RemoveProperty(char *key);
        // Remove a property from the table.

	G4MaterialPropertyVector* GetProperty(char *key);
        // Get the property from the table corresponding to the key-name. 

	void AddEntry(char *key, G4double aPhotonMomentum,
                                 G4double  aPropertyValue);
        // Add a new entry (pair of numbers) to the table for a given key.

	void RemoveEntry(char *key, G4double  aPhotonMomentum);
        // Remove an entry from the table for a given key and x-value.

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
