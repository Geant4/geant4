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
// G4UIcmdPargList.hh
//
// Specification of an abstract class of generic command argument
// that is a list with the following syntax:
//
//          (xxxxx,xxxxx,xxxxx,....,xxxxx)
//
// with no spaces and with xxxxxx as an object interpreted by
// pure virual method FetchNextItem.
//

#ifndef G4UIcmdPargList_hh
#define G4UIcmdPargList_hh

#include "G4UIcmdParg.hh"

class G4UIcmdPargList : public G4UIcmdParg {
	public:
	G4UIcmdPargList( const G4String &theName, G4int maxItem );
	virtual ~G4UIcmdPargList() {;}
	
	virtual G4String ConvertToString();

	virtual char GetTypeCode() const { return 's'; }
	virtual std::istream &FetchValue( std::istream &ios );
	
	inline G4int GetNItem() const { return nItem; }
	
	protected:
	virtual G4bool FetchItem( const char *string, const G4int item ) = 0;
	virtual std::ostream &WriteItem( std::ostream &ios, const G4int item ) = 0;
	
	G4int	nItem, maxItem;
};	


#endif
