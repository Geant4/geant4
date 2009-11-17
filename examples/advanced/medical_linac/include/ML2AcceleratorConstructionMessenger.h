#ifndef CML2AcceleratorConstructionMessengerH
#define CML2AcceleratorConstructionMessengerH


#include "globals.hh"
#include "G4UImessenger.hh"

class CML2AcceleratorConstruction;
class G4UImessenger;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

class CML2AcceleratorConstructionMessenger : public G4UImessenger 
{
public:
	CML2AcceleratorConstructionMessenger(CML2AcceleratorConstruction *acceleratorConstructor);
	~CML2AcceleratorConstructionMessenger(void);
	void SetNewValue(G4UIcommand* cmd, G4String newValue);
private:
	CML2AcceleratorConstruction *pAcceleratorConstructor;

	G4UIcmdWithAString *AcceleratorName, *acceleratorSpecficationsFile;
	G4UIcmdWithADoubleAndUnit *rotationX,  *rotationY,  *rotationZ;
};

#endif

