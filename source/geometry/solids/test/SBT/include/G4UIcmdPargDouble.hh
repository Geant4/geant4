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
// G4UIcmdPargDouble.hh
//
// Specification of a double argument for G4UIcmdWithPargs
//

#ifndef G4UIcmdPargDouble_hh
#define G4UIcmdPargDouble_hh

#include "G4UIcmdParg.hh"

class G4UIcmdPargDouble : public G4UIcmdParg {
	public:
	G4UIcmdPargDouble( const  G4String &theName, const G4double def, const G4double theUnits );
	virtual ~G4UIcmdPargDouble() {;}
	
	G4String ConvertToString();
	
	inline G4double GetValue() const { return value; }
	inline void SetValue( G4double newValue ) { value = newValue; }
	
	inline char GetTypeCode() const { return 'd'; }
	
	G4std::istream &FetchValue( G4std::istream &ios );
	
	protected:
	G4double units, value;
};

#endif
