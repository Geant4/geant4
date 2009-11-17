#ifndef CML2PhantomConstructionMessengerH
#define CML2PhantomConstructionMessengerH


#include "globals.hh"
#include "G4UImessenger.hh"

class CML2PhantomConstruction;
class G4UImessenger;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

class CML2PhantomConstructionMessenger : public G4UImessenger 
{
public:
	CML2PhantomConstructionMessenger(CML2PhantomConstruction *phantomConstructor);
	~CML2PhantomConstructionMessenger(void);
	void SetNewValue(G4UIcommand* cmd, G4String newValue);
private:
	CML2PhantomConstruction *pPhantomConstructor;

	G4UIcmdWithAnInteger *Phantom_nVoxelsX, *Phantom_nVoxelsY, *Phantom_nVoxelsZ;
	G4UIcmdWithAString *PhantomName, *PhantomSpecficationsFileName;
	G4UIcmdWithADoubleAndUnit *rotationX,  *rotationY,  *rotationZ;
};

#endif

