
#include "RadioactiveDecayGunmessenger.hh"

#include "g4std/iostream"
////////////////////////////////////////////////////////////////////////////////
//
RadioactiveDecayGunmessenger::RadioactiveDecayGunmessenger
(RadioactiveDecayGun* theRadioactiveDecayGun1) :
theRadioactiveDecayGun(theRadioactiveDecayGun1)
{
  ionCmd = new UIcmdWithNucleusAndUnit("/grdm/ion",this);
  ionCmd->SetGuidance("define the primary ion (a,z,e)");
  ionCmd->SetParameterName("A","Z","E",true);

  ionCmd->SetDefaultUnit("keV");
  ionCmd->SetUnitCandidates("keV MeV");
}
////////////////////////////////////////////////////////////////////////////////
//
RadioactiveDecayGunmessenger::~RadioactiveDecayGunmessenger ()
{
  delete ionCmd;
}
////////////////////////////////////////////////////////////////////////////////
//
void RadioactiveDecayGunmessenger::SetNewValue
  (G4UIcommand *command, G4String newValues)
{

  if (command==ionCmd) {theRadioactiveDecayGun->
    SetNucleus(ionCmd->GetNewNucleusValue(newValues));
 }
}
////////////////////////////////////////////////////////////////////////////////
//
//G4String G4RadioactiveDecayGunmessenger::GetCurrentValue (G4UIcommand * command)
//{

//  else if (command==NpIntvCmd)
//    {cv = NpIntvCmd->ConvertToString(theRDMnnun->GetNumberOfPhiIntv());}
//  return cv;
//}
////////////////////////////////////////////////////////////////////////////////










