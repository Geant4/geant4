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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef XrayFluoSteppingAction_h
#define  XrayFluoSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#ifdef G4ANALYSIS_USE
#include "XrayFluoAnalysisManager.hh"
#endif

class  XrayFluoDetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class  XrayFluoSteppingAction : public G4UserSteppingAction
{
public:
#ifdef G4ANALYSIS_USE    
  XrayFluoSteppingAction( XrayFluoDetectorConstruction*,XrayFluoAnalysisManager*);
#else
  XrayFluoSteppingAction();   
#endif
  ~ XrayFluoSteppingAction();
  
  void UserSteppingAction(const G4Step*);
private:
  XrayFluoDetectorConstruction* detector;
#ifdef G4ANALYSIS_USE   
  XrayFluoAnalysisManager* analysisManager;
#endif
};

#endif
