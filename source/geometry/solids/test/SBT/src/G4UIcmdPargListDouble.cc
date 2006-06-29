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
// G4UIcmdPargListDouble.cc
//
// Implementation of an argument consisting of
// a list of double values for G4UIcmdWithParg
//

#include "G4UIcmdPargListDouble.hh"

#include <sstream>

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
	std::istringstream is( string );
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
std::ostream &G4UIcmdPargListDouble::WriteItem( std::ostream &ios, const G4int item )
{
	ios << storage[item];
	return ios;
}
