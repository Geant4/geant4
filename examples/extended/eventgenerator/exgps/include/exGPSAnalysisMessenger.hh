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
#ifndef exGPSAnalysisMessenger_h
#define exGPSAnalysisMessenger_h 1

#ifdef G4ANALYSIS_USE

#include "globals.hh"
#include "G4UImessenger.hh"

class exGPSAnalysisManager;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class exGPSAnalysisMessenger: public G4UImessenger
{
public:
  exGPSAnalysisMessenger(exGPSAnalysisManager* );
  ~exGPSAnalysisMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  exGPSAnalysisManager* exGPSAnalysis;
  G4UIdirectory*              exGPSAnalysisDir;
  
  G4UIcmdWithAString*               FileNameCmd;
  G4UIcmdWithAString*               FileTypeCmd;
  G4UIcmdWithADoubleAndUnit*        MaxEngCmd;
  G4UIcmdWithADoubleAndUnit*        MinEngCmd;
  G4UIcmdWithADoubleAndUnit*        MaxPosCmd;
  G4UIcmdWithADoubleAndUnit*        MinPosCmd;
};

#endif // G4ANALYSIS_USE

#endif










