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
// G4UIcmdParg.hh
//
// Specification of a pure virtual class of generic command
// arguments, to be used by G4UIcmdWithPargs.
//

#ifndef G4UIcmdParg_hh
#define G4UIcmdParg_hh

#include "globals.hh"
#include "g4std/iostream"

class G4UIcmdParg {
	public: 
	G4UIcmdParg( const G4String &theName ) { name = theName; }
	virtual ~G4UIcmdParg() {;}
	
	virtual G4String ConvertToString() = 0;
	
	virtual char GetTypeCode() const = 0;
	inline G4String GetName() const { return name; }
	
	virtual G4std::istream &FetchValue( G4std::istream &ios ) = 0;
	
	protected:
	G4String name;
};


#endif
