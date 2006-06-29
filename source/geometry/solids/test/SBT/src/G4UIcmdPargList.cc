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
// G4UIcmdPargList.cc
//
// Implementation of a (generic) list argument for G4UIcmdWithPargs
//

#include "G4UIcmdPargList.hh"

#include <sstream>
#include <iomanip>

//
// Constructor
//
G4UIcmdPargList::G4UIcmdPargList( const  G4String &theName, G4int theMaxItem ) 
		: G4UIcmdParg( theName )
{
	nItem = 0;
	maxItem = theMaxItem;
}


//
// Fetch argument value from input stream
//
std::istream &G4UIcmdPargList::FetchValue( std::istream &istr )
{
	//
	// Because of limitations of G4UIcommand, we must avoid
	// requiring any spaces.
	//
	// Require '(' to be the leading character, and read in
	// character strings until ',' or ')' are found. If the
	// latter is found, input is terminated.
	//
	// In actual usage, we expect leading blanks to be removed
	// before this method is invoked,
	// but there is no need to assume that.
	//
	// Notice that blank items are allowed. It is up to decendent
	// classes to decide what to do in this case. To allow an
	// empty list (i.e. no items, nItem==0), the special symbol
	// "-" is reserved.
	//
	char buffer[256];	// Buffer for each item
	char a;			// Single character input
	
	//
	// Remove leading blanks
	//
	do {
		istr.get(a);
		if (istr.fail()) return istr;
	} while( a == ' ' );
	
	//
	// Check for empty list ( '-', followed by EOF or blank)
	//
	if (a == '-') {
		istr.get(a);
		if (istr.fail() || a == ' ') {
			nItem = 0;
			return istr;
		}
	}
		
	
	//
	// Otherwise: Insist argument begin with a '('
	//
	if (a != '(') {
		G4cerr << "Listed argument must begin with '(' or consist of only '-'" << G4endl;
		istr.clear(std::ios::failbit|istr.rdstate());	// Is there a better way to do this???
		return istr;
	}

	//
	// Loop over items in argument
	//
	G4int nItemFetched = 0;
	do {
		//
		// Next item
		//
		char *b = buffer;
		for(;;) {
			//
			// Next character
			//
			istr.get(*b);
			if (istr.fail()) return istr;
			
			//
			// Check for delimiter
			//
			if (*b == ',' || *b == ')') {
				a = *b;		// Remember 
				*b = '\0';
				break;
			}
			
			//
			// Go to next character
			//
			if (++b > buffer+255) {
				G4cerr << "Listed argument item must be less than 255 characters" << G4endl;
				istr.clear(std::ios::failbit|istr.rdstate());
				return istr;
			}
		} 
		
		//
		// Interpret
		//
		if (nItemFetched >= maxItem) {
			G4cerr << "Maximum of " << maxItem << " items exceeded in listed argument" << G4endl;
			istr.clear(std::ios::failbit|istr.rdstate());
			return istr;
		}
		
		if (!FetchItem( buffer, nItemFetched++ )) {
			istr.clear(std::ios::failbit|istr.rdstate());
			return istr;
		}
	} while( a == ',' );
	
	//
	// Set nItem only on success
	//
	nItem = nItemFetched;
	
	return istr;
}


//
// ConvertToString
//
G4String G4UIcmdPargList::ConvertToString()
{
	//
	// Empty case
	//
	if (nItem == 0) return "-";
	
	//
	// Allocate a buffer, big enough to hold all items plus
	// a command and '()'. Each item is allowed to take 255 
	// characters.
	//
	// buffSize = 255*nItem+4;
	
	std::ostringstream os;
	
	//
	// Write out everything in turn
	//
	os << '(';
	
	G4int item;
	for( item=0; item<nItem; item++) {
		if (item) os << ',';
		if (!WriteItem(os,item)) return "(error)";
	}
	
	os << ')';
	
	//
	// Convert to G4String
	//
	G4String answer = os.str();
		
	//
	// Return G4String
	//
	return answer;
}
