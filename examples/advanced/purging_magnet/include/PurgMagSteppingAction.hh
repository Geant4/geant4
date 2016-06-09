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
// Code developed by:
//  S.Larsson
//
//    ***********************************
//    *                                 *
//    *    PurgMagSteppingAction.hh     *
//    *                                 *
//    ***********************************
//
// $Id: PurgMagSteppingAction.hh,v 1.2 2004/06/18 09:17:50 gunter Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef PurgMagSteppingAction_h
#define PurgMagSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4VPVParameterisation.hh"
#include "G4PVParameterised.hh"
#include "G4Tubs.hh"


class PurgMagRunAction;
class PurgMagDetectorConstruction;
class PurgMagAnalysisManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class PurgMagSteppingAction : public G4UserSteppingAction
{
public:
  PurgMagSteppingAction(PurgMagRunAction*,PurgMagDetectorConstruction*);
  ~PurgMagSteppingAction();
  
  void UserSteppingAction(const G4Step*);
  
private:
  PurgMagRunAction*            PurgMagRun;
  PurgMagDetectorConstruction* Detector; 
  
};

#endif




