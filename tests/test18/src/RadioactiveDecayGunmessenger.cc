//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//

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










