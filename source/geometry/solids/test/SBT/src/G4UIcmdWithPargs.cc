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
// G4UIcmdWithPargs.cc
//
// Implementation of a G4UIcommand with generic arguments
//

#include "G4UIcmdWithPargs.hh"
#include "G4UIcmdParg.hh"

#include <sstream>

//
// Constructor
//
G4UIcmdWithPargs::G4UIcmdWithPargs( const char *path, G4UImessenger *theMessenger,
			  	    G4UIcmdParg **theArgs, const int theNumParg )
		: G4UIcommand( path, theMessenger )
{
	args = theArgs;
	numArgs = theNumParg;
	
	G4UIcmdParg **thisArg;
	G4int index;
	for( index=0, thisArg=args; index < numArgs; index++, thisArg++ ) {
		//
		// Boy, sure looks like a memory leak, but it isn't.
		// See G4UIcommand::~G4UIcommand
		//
		G4UIparameter *param = new G4UIparameter( (*thisArg)->GetTypeCode() );
		param->SetParameterName( (*thisArg)->GetName() );
		param->SetOmittable(false);
		SetParameter( param );
	}
}


//
// Destructor
//
G4UIcmdWithPargs::~G4UIcmdWithPargs() {;}


//
// GetArguments
//
// Unlike the standard G4UIcommand, the values are returned in the
// in the array of arguments given to the constructor
//
G4bool G4UIcmdWithPargs::GetArguments( G4String argumentString )
{
	const char *buff = argumentString;
	std::istringstream is( buff );
	
	G4UIcmdParg **thisArg = args;
	while( thisArg < args + numArgs ) {
		(*thisArg++)->FetchValue(is);
		if (is.fail()) return false;
	}
	
	return true;
}
