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
// Author: Patricia Mendez (patricia.mendez@cern.ch)
//
// History:
// -----------
//  14 Feb 2003  Patricia Mendez   Created
//
// -------------------------------------------------------------------




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4ANALYSIS_USE
#ifndef FCALAnalysisMessenger_h
#define FCALAnalysisMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"


class FCALAnalysisManager;
class G4UIdirectory;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


class G4UIcmdWithAString;
class FCALAnalysisMessenger: public G4UImessenger

{
public:
  FCALAnalysisMessenger(FCALAnalysisManager* );
  ~FCALAnalysisMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:

  //pointer to FCALAnalysisManager
  FCALAnalysisManager* xrayFluoAnalysis;
  G4UIdirectory* FCALAnalysisDir;
  G4UIcmdWithAString* ouputFileCommand;

};
#endif
#endif


