#ifndef CML2Acc1MessengerH
#define CML2Acc1MessengerH


#include "globals.hh"
#include "G4UImessenger.hh"

class CML2Acc1;
class G4UImessenger;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

class CML2Acc1Messenger : public G4UImessenger 
{
public:
	CML2Acc1Messenger(CML2Acc1 *acc1);
	~CML2Acc1Messenger(void);
	void SetNewValue(G4UIcommand* cmd, G4String newValue);
private:
	CML2Acc1 *pAcc1;

	G4UIcmdWithAnInteger *idEnergy;
	G4UIcmdWithADoubleAndUnit *aperture1X,*aperture2X, *aperture1Y, *aperture2Y, *SSD, *leavesA, *leavesB;
};

#endif

