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
// G4UIcmdWithPargs.cc
//
// Implementation of a G4UIcommand with generic arguments
//

#include "G4UIcmdWithPargs.hh"
#include "G4UIcmdParg.hh"

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
	G4std::istrstream is( (char *)buff );
	
	G4UIcmdParg **thisArg = args;
	while( thisArg < args + numArgs ) {
		(*thisArg++)->FetchValue(is);
		if (is.fail()) return false;
	}
	
	return true;
}
