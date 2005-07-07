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
// $Id: XrayFluoRunActionMessenger.cc
// GEANT4 tag $Name: xray_fluo-V04-01-03
//
//
// History:
// -----------
//
// -------------------------------------------------------------------

#include "Tst51RunActionMessenger.hh"
#include "Tst51RunAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

Tst51RunActionMessenger::Tst51RunActionMessenger(Tst51RunAction * run)
:fTst51RunAction(run)
{
  analysisDir = new G4UIdirectory("/analysis/");
  analysisDir->SetGuidance("analysis commands");

  fileNameCmd = new G4UIcmdWithAString("/analysis/filename",this);  
  fileNameCmd->SetGuidance("Change the name of the hbook file");
}

Tst51RunActionMessenger::~Tst51RunActionMessenger()
{  
  delete analysisDir;
  delete fileNameCmd;
}
  
void Tst51RunActionMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{
  if (command == fileNameCmd)
   { fTst51RunAction->SetFileName(newValue); }
}







