#include "ML2Ph_BoxInBoxMessenger.h"
#include "ML2Ph_BoxInBox.h"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

CML2Ph_BoxInBoxMessenger::CML2Ph_BoxInBoxMessenger(CML2Ph_BoxInBox *Ph_BoxInBox) : pPh_BoxInBox(Ph_BoxInBox)
{
}

CML2Ph_BoxInBoxMessenger::~CML2Ph_BoxInBoxMessenger(void)
{
}
void CML2Ph_BoxInBoxMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue)
{

}
