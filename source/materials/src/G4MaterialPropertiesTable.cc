// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4MaterialPropertiesTable.cc,v 1.1 1999-01-07 16:09:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
////////////////////////////////////////////////////////////////////////
// G4MaterialPropertiesTable Implementation
////////////////////////////////////////////////////////////////////////
//
// File: G4MaterialPropertiesTable.cc 
// Version:     1.0
// Created:     1996-02-08
// Author:      Juliet Armstrong
// Updated:     1997-03-26 by Peter Gumplinger
//              > cosmetics (only)
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#include "G4MaterialPropertiesTable.hh"

unsigned hashString(const G4String &str) { return str.hash(); }

        //////////////
        // Operators
        //////////////

G4MaterialPropertiesTable&
G4MaterialPropertiesTable::operator =(const G4MaterialPropertiesTable &right)
{
	// clear any current contents of MPT

        MPT.clearAndDestroy();

        // want to make an actual copy -- not a shallow copy which is
	// the default for RWTPrtHashDictionary's assignment operator
 
        RWTPtrHashDictionary<G4String, G4MaterialPropertyVector> 
						rightMPT(right.MPT);
        RWTPtrHashDictionaryIterator<G4String, G4MaterialPropertyVector> 
						rightIterator(rightMPT); 
        rightIterator.reset();
        while (++rightIterator) {
		G4MaterialPropertyVector *newProp =
                        new G4MaterialPropertyVector(*(rightIterator.value()));
                RWCString *newKey =
                        new RWCString(*(rightIterator.key()));
                MPT.insertKeyAndValue(newKey, newProp);
        }
        return *this;
}

        /////////////////
        // Constructors
        /////////////////

G4MaterialPropertiesTable::G4MaterialPropertiesTable() : MPT(hashString),
                                                         MPTiterator(MPT) {}

G4MaterialPropertiesTable::G4MaterialPropertiesTable
			   (const G4MaterialPropertiesTable &right) : 
			   MPT(hashString), MPTiterator(MPT)
{
        // want to make an actual copy -- not a shallow copy which is
	// the default for RWTPrtHashDictionary's assignment operator

        RWTPtrHashDictionary<G4String, G4MaterialPropertyVector> 
						rightMPT(right.MPT);
        RWTPtrHashDictionaryIterator<G4String, G4MaterialPropertyVector> 
						rightIterator(rightMPT); 

        rightIterator.reset();

        while (++rightIterator) {
		G4MaterialPropertyVector *newProp =
                        new G4MaterialPropertyVector(*(rightIterator.value()));
                RWCString *newKey =
                        new RWCString(*(rightIterator.key()));
                MPT.insertKeyAndValue(newKey, newProp);
        }
}

        ////////////////
        // Destructors
        ////////////////

G4MaterialPropertiesTable::~G4MaterialPropertiesTable()
{
	MPT.clearAndDestroy();
}

        ////////////
        // Methods
        ////////////

void G4MaterialPropertiesTable::AddProperty(char     *key,
					    G4double *PhotonMomenta,
					    G4double *PropertyValues,
					    G4int     NumEntries)
{
	G4MaterialPropertyVector *mpv = 
			new G4MaterialPropertyVector(PhotonMomenta, 
					  	     PropertyValues, 
						     NumEntries);
	G4String *newKey = new G4String(key);
	MPT.insertKeyAndValue(newKey, mpv);
}

void G4MaterialPropertiesTable::AddProperty(char *key,
					    G4MaterialPropertyVector *mpv)
{
//	Provides a way of adding a property to the Material Properties
//	Table given an G4MaterialPropertyVector Reference and a key 

	G4String *theKey = new G4String(key);
	MPT.insertKeyAndValue(theKey, mpv);	
} 

void G4MaterialPropertiesTable::RemoveProperty(char *key)
{
	G4String target(key);
	MPT.remove(&target);
}

G4MaterialPropertyVector* G4MaterialPropertiesTable::GetProperty(char *key)
{
	G4String target(key);

	return MPT.findValue(&target);
}

void G4MaterialPropertiesTable::AddEntry(char     *key,
					 G4double  aPhotonMomentum,
					 G4double  aPropertyValue)
{
	G4String target(key);
	G4MaterialPropertyVector *targetVector=MPT.findValue(&target);
	if (targetVector != NULL) {
		targetVector->AddElement(aPhotonMomentum, aPropertyValue);
	}
	else {
		G4Exception("G4MaterialPropertiesTable::AddEntry ==> "
			    "Material Property Vector not found.");
	}
}

void G4MaterialPropertiesTable::RemoveEntry(char *key, 
					    G4double  aPhotonMomentum)
{
        G4String target(key);
        G4MaterialPropertyVector *targetVector=MPT.findValue(&target);
	if (targetVector) {
		targetVector->RemoveElement(aPhotonMomentum);
 	}
        else {
                G4Exception("G4MaterialPropertiesTable::AddEntry ==> "
			    "Material Property Vector not found.");
        }
}
void G4MaterialPropertiesTable::DumpTable()
{
	MPTiterator.reset();
	while(++MPTiterator) {
		G4cout << *MPTiterator.key() << endl;
		MPTiterator.value()->DumpVector();
	}
}
