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
// $Id: RemSimAnalysisMessenger.hh
// GEANT4 tag $Name: radioprotection-V07-00-01
//
// Author: Susanna Guatelli (Susanna.Guatelli@ge.infn.it)
//
// History:
// -----------
//  23 Nov 2004  Susanna Guatelli   Created
//
// -------------------------------------------------------------------




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4ANALYSIS_USE
#ifndef RemSimAnalysisMessenger_h
#define RemSimAnalysisMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "RemSimAnalysisManager.hh"
#include "G4UIdirectory.hh"
class RemSimAnalysisManager;
//class G4UIdirectory;
//class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


class RemSimAnalysisMessenger: public G4UImessenger

{
public:
  RemSimAnalysisMessenger(RemSimAnalysisManager*);
  ~RemSimAnalysisMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  RemSimAnalysisManager* remsimAnalysis;
  G4UIdirectory* analysisDir;
  G4UIcmdWithAString* outputFileType;
};
#endif
#endif


