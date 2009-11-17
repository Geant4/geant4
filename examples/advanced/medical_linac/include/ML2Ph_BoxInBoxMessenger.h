#ifndef ML2Ph_BoxInBoxMessengerH
#define ML2Ph_BoxInBoxMessengerH


#include "globals.hh"
#include "G4UImessenger.hh"

class CML2Ph_BoxInBox;
class G4UImessenger;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

class CML2Ph_BoxInBoxMessenger : public G4UImessenger 
{
public:
	CML2Ph_BoxInBoxMessenger(CML2Ph_BoxInBox *Ph_BoxInBox);
	~CML2Ph_BoxInBoxMessenger(void);
	void SetNewValue(G4UIcommand* cmd, G4String newValue);
private:
	CML2Ph_BoxInBox *pPh_BoxInBox;
};

#endif

