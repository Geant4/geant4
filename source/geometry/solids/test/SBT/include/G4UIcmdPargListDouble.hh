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
	virtual ostream &WriteItem( ostream &ios, const G4int item );
	
	G4double units, *storage;
};

#endif
