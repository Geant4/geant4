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
	virtual G4std::istream &FetchValue( G4std::istream &ios );
	
	inline G4int GetNItem() const { return nItem; }
	
	protected:
	virtual G4bool FetchItem( const char *string, const G4int item ) = 0;
	virtual G4std::ostream &WriteItem( G4std::ostream &ios, const G4int item ) = 0;
	
	G4int	nItem, maxItem;
};	


#endif
