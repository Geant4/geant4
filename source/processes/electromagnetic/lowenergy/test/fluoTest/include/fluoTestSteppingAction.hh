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

#ifndef fluoTestSteppingAction_h
#define fluoTestSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#ifdef G4ANALYSIS_USE
#include "fluoTestAnalysisManager.hh"
#endif

class fluoTestDetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class fluoTestSteppingAction : public G4UserSteppingAction
{
  public:
#ifdef G4ANALYSIS_USE    
  fluoTestSteppingAction(fluoTestDetectorConstruction*,fluoTestAnalysisManager*);
#else
  fluoTestSteppingAction();   
#endif
   ~fluoTestSteppingAction();

    void UserSteppingAction(const G4Step*);
 private:
    fluoTestDetectorConstruction* detector;
#ifdef G4ANALYSIS_USE   
    fluoTestAnalysisManager* analysisManager;
#endif
};

#endif
