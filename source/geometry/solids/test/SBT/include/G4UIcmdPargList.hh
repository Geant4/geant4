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
	virtual istream &FetchValue( istream &ios );
	
	inline G4int GetNItem() const { return nItem; }
	
	protected:
	virtual G4bool FetchItem( const char *string, const G4int item ) = 0;
	virtual ostream &WriteItem( ostream &ios, const G4int item ) = 0;
	
	G4int	nItem, maxItem;
};	


#endif
