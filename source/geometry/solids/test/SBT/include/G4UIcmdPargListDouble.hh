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
// G4UIcmdPargListDouble.hh
//
// Declaration of an argument consisting of
// a list of double values for G4UIcmdWithParg
//

#ifndef G4UIcmdPargListDouble_hh
#define G4UIcmdPargListDouble_hh

#include "G4UIcmdPargList.hh"

class G4UIcmdPargListDouble : public G4UIcmdPargList {
	public:
	G4UIcmdPargListDouble( const G4String &theName, G4int maxItem, const G4double units );
	~G4UIcmdPargListDouble();
	
	inline G4double *GetValues() const { return storage; }
	
	protected:
	virtual G4bool FetchItem( const char *string, const G4int item );
	virtual std::ostream &WriteItem( std::ostream &ios, const G4int item );
	
	G4double units, *storage;
};

#endif
