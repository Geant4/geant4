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
