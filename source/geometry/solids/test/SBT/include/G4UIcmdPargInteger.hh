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
// G4UIcmdPargInteger.hh
//
// Specification of a integer argument for G4UIcmdWithPargs
//

#ifndef G4UIcmdPargInteger_hh
#define G4UIcmdPargInteger_hh

#include "G4UIcmdParg.hh"

class G4UIcmdPargInteger : public G4UIcmdParg {
	public:
	G4UIcmdPargInteger( const  G4String &theName, const G4int def ) : G4UIcmdParg(theName) { value = def; }
	virtual ~G4UIcmdPargInteger() {;}
	
	G4String ConvertToString();
	
	inline G4int GetValue() const { return value; }
	inline void SetValue( G4int newValue ) { value = newValue; }
	
	inline char GetTypeCode() const { return 'i'; }
	
	G4std::istream &FetchValue( G4std::istream &ios );
	
	protected:
	G4int value;
};

#endif
