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
	
	std::istream &FetchValue( std::istream &ios );
	
	protected:
	G4double units, value;
};

#endif
