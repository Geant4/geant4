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
	virtual G4std::ostream &WriteItem( G4std::ostream &ios, const G4int item );
	
	G4double units, *storage;
};

#endif
