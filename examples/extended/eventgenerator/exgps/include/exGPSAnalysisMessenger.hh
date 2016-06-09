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










