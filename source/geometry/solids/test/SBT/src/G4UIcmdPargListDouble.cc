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
// G4UIcmdPargListDouble.cc
//
// Implementation of an argument consisting of
// a list of double values for G4UIcmdWithParg
//

#include "G4UIcmdPargListDouble.hh"

#include "g4std/strstream"

//
// Constructor
//
G4UIcmdPargListDouble::G4UIcmdPargListDouble( const G4String &theName, G4int theMaxItem, 
					      const G4double theUnits )
			: G4UIcmdPargList( theName, theMaxItem )
{
	units = theUnits;
	//
	// Allocate memory for maximum number of items
	//
	storage = new G4double [theMaxItem];
}


//
// Destructor
//
G4UIcmdPargListDouble::~G4UIcmdPargListDouble()
{
	delete [] storage;
}


//
// FetchItem
//
G4bool G4UIcmdPargListDouble::FetchItem( const char *string, const G4int item )
{
	G4std::istrstream is( (char *)string );
	is >> storage[item];
	if (is.fail()) {
		G4cerr << "Error interpreting '" << string << "' as a double value" << G4endl;
		return false;
	}
	storage[item] *= units;
	
	return true;
}


//
// WriteItem
//
G4std::ostream &G4UIcmdPargListDouble::WriteItem( G4std::ostream &ios, const G4int item )
{
	ios << storage[item];
	return ios;
}
