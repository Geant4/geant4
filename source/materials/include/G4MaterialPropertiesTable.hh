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
// $Id: G4MaterialPropertiesTable.hh,v 1.11 2002-11-07 02:29:58 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
////////////////////////////////////////////////////////////////////////
// G4MaterialPropertiesTable Definition
////////////////////////////////////////////////////////////////////////
//
// File:        G4MaterialPropertiesTable.hh
// Description: An Material properties table is a hash table, with 
// 		key = property name, and value either G4double or
//              G4MaterialPropertyVector 
// Version:     1.0
// Created:     1996-02-08
// Author:      Juliet Armstrong
// Updated:     2002-11-05 add named material constants by P. Gumplinger
//              1999-11-05 Migration from G4RWTPtrHashDictionary to STL
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
// name, and value either G4double or G4MaterialPropertyVector.
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////

class G4MaterialPropertiesTable {

        ////////////////
        // Constructor
        ////////////////

public: // Without description

	G4MaterialPropertiesTable(); 
	
	///////////////
	// Destructor
	///////////////

public: // Without description

	~G4MaterialPropertiesTable();

	////////////
	// Methods
	////////////

public: // With description

        void AddConstProperty(const char     *key,
                              G4double PropertyValue);
        // Add a new property to the table by giving a key-name and value 

	void AddProperty(const char     *key,
		         G4double *PhotonMomenta,
		         G4double *PropertyValues,
		         G4int     NumEntries);
        // Add a new property to the table by giving a key-name and the
        // arrays x and y of size NumEntries.

	void AddProperty(const char *key, G4MaterialPropertyVector *opv);
        // Add a new property to the table by giving a key-name and an
        // already constructed G4MaterialPropertyVector.

        void RemoveConstProperty(const char *key);
        // Remove a constant property from the table.

	void RemoveProperty(const char *key);
        // Remove a property from the table.

        G4double GetConstProperty(const char *key);
        // Get the constant property from the table corresponding to the key-name

	G4MaterialPropertyVector* GetProperty(const char *key);
        // Get the property from the table corresponding to the key-name.

	void AddEntry(const char *key, G4double aPhotonMomentum,
                                 G4double  aPropertyValue);
        // Add a new entry (pair of numbers) to the table for a given key.

	void RemoveEntry(const char *key, G4double  aPhotonMomentum);
        // Remove an entry from the table for a given key and x-value.

	void DumpTable();

private:

	/////////////////////////
	// Private Data members
	/////////////////////////

	G4std::map<G4String, G4MaterialPropertyVector*, G4std::less<G4String> > MPT;
	typedef G4std::map<G4String, G4MaterialPropertyVector*,
	  G4std::less<G4String> >::iterator MPTiterator;

        G4std::map< G4String, G4double, G4std::less<G4String> > MPTC;
        typedef G4std::map< G4String, G4double,
          G4std::less<G4String> >::iterator MPTCiterator;

};

#endif /* G4MaterialPropertiesTable_h */
