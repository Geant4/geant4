// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MaterialPropertiesTable.hh,v 1.8 1999-12-15 14:50:49 gunter Exp $
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
// Updated:     1999-11-05 Migration from G4RWTPtrHashDictionary to STL
//                         by John Allison
//              1999-10-29 add method and class descriptors
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
#include "g4std/map"
#include "globals.hh"
#include "G4MaterialPropertyVector.hh"

// Class Description:
// An Material properties table is a hash table, with key = property
// name, and value = G4MaterialPropertyVector.
// Class Description - End:

#define RefractionIndex "RINDEX"

/////////////////////
// Class Definition
/////////////////////

class G4MaterialPropertiesTable {

	//////////////
        // Operators
        //////////////

private:

	G4MaterialPropertiesTable&
		operator =(const G4MaterialPropertiesTable &right);
		
	////////////////
	// Constructor
	////////////////

public: // Without description

	G4MaterialPropertiesTable(); 
	
private:

	G4MaterialPropertiesTable(const G4MaterialPropertiesTable &right);

public: // Without description

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

	G4std::map<G4String, G4MaterialPropertyVector*, G4std::less<G4String> > MPT;
	typedef G4std::map<G4String, G4MaterialPropertyVector*,
	  G4std::less<G4String> >::iterator MPTiterator;
};

#endif /* G4MaterialPropertiesTable_h */
