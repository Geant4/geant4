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

  InteractPlotCmd = new G4UIcmdWithABool("/dmx/PlotAtEnd",this);
  InteractPlotCmd->SetGuidance("Flag for Interactive plotting at end of run");
  InteractPlotCmd->SetGuidance("Default = false");
  InteractPlotCmd->SetParameterName("interactplot", false);
  InteractPlotCmd->SetDefaultValue(false);

  //  FileCmd->AvailableForStates(Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DMXRunActionMessenger::~DMXRunActionMessenger()
{
  delete SaveHitsCmd;  
  delete SavePmtCmd;  
  delete SaveHistFileCmd;  
  delete InteractPlotCmd;
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

  if(command == InteractPlotCmd) {
    G4int vl;
    const char* t = newValue;
    G4std::istrstream is((char*)t);
    is >> vl;
    DMXRun->Setinteractplot(vl!=0);
  }

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....





