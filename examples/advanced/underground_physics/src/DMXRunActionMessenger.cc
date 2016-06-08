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
//
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: a.s.howard@ic.ac.uk
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// RunActionMessenger program
// --------------------------------------------------------------

#include "DMXRunActionMessenger.hh"
#include "DMXRunAction.hh"

#include "G4UIcmdWithAString.hh"
#include "G4ios.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DMXRunActionMessenger::DMXRunActionMessenger(DMXRunAction* run)
:DMXRun(run)
{ 
  FileCmd = new G4UIcmdWithAString("/run/filename",this);
  FileCmd->SetGuidance(" The log file name for the run");
  FileCmd->SetGuidance("   default = rdmex2.log ");
  FileCmd->SetParameterName(" Input ",true);
  FileCmd->SetDefaultValue("rdmex2.log");
  FileCmd->AvailableForStates(Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DMXRunActionMessenger::~DMXRunActionMessenger()
{
  delete FileCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXRunActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
    if(command == FileCmd)
    {DMXRun->SetFilename(newValue);}

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....





