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
	
	istream &FetchValue( istream &ios );
	
	protected:
	G4int value;
};

#endif
