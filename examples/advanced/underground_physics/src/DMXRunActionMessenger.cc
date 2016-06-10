//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: alexander.howard@cern.ch
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

#include <sstream>

#include "DMXRunAction.hh"

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4ios.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DMXRunActionMessenger::DMXRunActionMessenger(DMXRunAction* run)
:DMXRun(run)
{ 
  SaveHitsCmd = new G4UIcmdWithAString("/dmx/hitsfile",this);
  SaveHitsCmd->SetGuidance("output file for hits collection (txt)");
  SaveHitsCmd->SetGuidance("Default = hits.out");
  SaveHitsCmd->SetParameterName("savehitsFile", false);
  SaveHitsCmd->SetDefaultValue("hits.out");

  SavePmtCmd = new G4UIcmdWithAString("/dmx/pmtfile",this);
  SavePmtCmd->SetGuidance("output file for pmt hits (txt)");
  SavePmtCmd->SetGuidance("Default = pmt.out");
  SavePmtCmd->SetParameterName("savepmtFile", false);
  SavePmtCmd->SetDefaultValue("pmt.out");

  SaveHistFileCmd = new G4UIcmdWithAString("/dmx/histogramfile",this);
  SaveHistFileCmd->SetGuidance("output file for histograms");
  SaveHistFileCmd->SetGuidance("Default = dmx.his");
  //  SaveHistFileCmd->SetParameterName("savehistFile", false);
  SaveHistFileCmd->SetParameterName("histFile", false);
  SaveHistFileCmd->SetDefaultValue("dmx.his");


  //  FileCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DMXRunActionMessenger::~DMXRunActionMessenger()
{
  delete SaveHitsCmd;  
  delete SavePmtCmd;  
  delete SaveHistFileCmd;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXRunActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if(command == SaveHitsCmd)
    DMXRun->SetsavehitsFile(newValue);

  if(command == SavePmtCmd)
    DMXRun->SetsavepmtFile(newValue);

  if(command == SaveHistFileCmd)
    DMXRun->SetsavehistFile(newValue);


}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....





