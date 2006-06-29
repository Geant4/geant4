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
	
	std::istream &FetchValue( std::istream &ios );
	
	protected:
	G4int value;
};

#endif
