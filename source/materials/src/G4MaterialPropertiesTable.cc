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
// $Id: G4MaterialPropertiesTable.cc,v 1.13 2002-11-07 02:30:29 gum Exp $
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
// Updated:     2002-11-05 add named material constants by P. Gumplinger
//              1999-11-05 Migration from G4RWTPtrHashDictionary to STL
//                         by John Allison
//              1997-03-26 by Peter Gumplinger
//              > cosmetics (only)
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#include "globals.hh"
#include "G4MaterialPropertiesTable.hh"

        /////////////////
        // Constructors
        /////////////////

G4MaterialPropertiesTable::G4MaterialPropertiesTable() {}

        ////////////////
        // Destructors
        ////////////////

G4MaterialPropertiesTable::~G4MaterialPropertiesTable()
{
        MPTiterator i;
        for (i = MPT.begin(); i != MPT.end(); ++i) {
            delete (*i).second;
        }
        MPT.clear();
        MPTC.clear();
}

        ////////////
        // Methods
        ////////////

void G4MaterialPropertiesTable::AddConstProperty(const char     *key,
                                                 G4double PropertyValue)
{
//      Provides a way of adding a constant property to the Mataerial Properties
//      Table given a key

        MPTC [G4String(key)] = PropertyValue;
}

void G4MaterialPropertiesTable::AddProperty(const char     *key,
					    G4double *PhotonMomenta,
					    G4double *PropertyValues,
					    G4int     NumEntries)
{
//      Privides a way of adding a property to the Material Properties
//      Table given a pair of numbers and a key

	G4MaterialPropertyVector *mpv = 
			new G4MaterialPropertyVector(PhotonMomenta, 
					  	     PropertyValues, 
						     NumEntries);
	MPT [G4String(key)] = mpv;
}

void G4MaterialPropertiesTable::AddProperty(const char *key,
					    G4MaterialPropertyVector *mpv)
{
//	Provides a way of adding a property to the Material Properties
//	Table given an G4MaterialPropertyVector Reference and a key 

	MPT [G4String(key)] = mpv;
} 

void G4MaterialPropertiesTable::RemoveConstProperty(const char *key)
{
        MPTC.erase(G4String(key));
}

void G4MaterialPropertiesTable::RemoveProperty(const char *key)
{
	MPT.erase(G4String(key));
}

G4double G4MaterialPropertiesTable::GetConstProperty(const char *key)
{
//      Returns the constant material property corresponding to a key

        MPTCiterator j;
        j = MPTC.find(G4String(key));
        if ( j != MPTC.end() ) {
           return j->second;
        } 
        else {
           G4Exception("G4MaterialPropertiesTable::GetConstProperty ==> "
                       "Constant Material Property not found.");
           return G4double(0.0);
        }
}

G4MaterialPropertyVector* G4MaterialPropertiesTable::GetProperty(const char *key)
{
//      Returns a Material Property Vector corresponding to a key

        return MPT [G4String(key)];
}

void G4MaterialPropertiesTable::AddEntry(const char     *key,
					 G4double  aPhotonMomentum,
					 G4double  aPropertyValue)

//      Allows to add an entry pair directly to the Material Property Vector
//      given a key

{
	G4MaterialPropertyVector *targetVector=MPT [G4String(key)];
	if (targetVector != 0) {
		targetVector->AddElement(aPhotonMomentum, aPropertyValue);
	}
	else {
		G4Exception("G4MaterialPropertiesTable::AddEntry ==> "
			    "Material Property Vector not found.");
	}
}

void G4MaterialPropertiesTable::RemoveEntry(const char *key,  
					    G4double  aPhotonMomentum)
{
//      Allows to remove an entry pair directly from the Material Property Vector
//      given a key

        G4MaterialPropertyVector *targetVector=MPT [G4String(key)];
	if (targetVector) {
		targetVector->RemoveElement(aPhotonMomentum);
 	}
        else {
                G4Exception("G4MaterialPropertiesTable::RemoveEntry ==> "
			    "Material Property Vector not found.");
        }
}
void G4MaterialPropertiesTable::DumpTable()
{
  MPTiterator i;
  for (i = MPT.begin(); i != MPT.end(); ++i) {
		G4cout << (*i).first << G4endl;
                if ( (*i).second != 0 ) {
		  (*i).second->DumpVector();
                }
                else {
                  G4cout << "NULL Material Property Vector Pointer." << G4endl;
                }
  }
  MPTCiterator j;
  for (j = MPTC.begin(); j != MPTC.end(); ++j) {
                 G4cout << j->first << G4endl;
                 if ( j->second != 0 ) {
                   G4cout << j->second << G4endl;
                 }
                 else {
                   G4cout << "No Material Constant Property." << G4endl;
                 }
  }

}
