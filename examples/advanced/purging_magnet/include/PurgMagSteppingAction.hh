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
// Code developed by:
//  S.Larsson
//
//    ***********************************
//    *                                 *
//    *    PurgMagSteppingAction.hh     *
//    *                                 *
//    ***********************************
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef PurgMagSteppingAction_h
#define PurgMagSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4VPVParameterisation.hh"
#include "G4PVParameterised.hh"
#include "G4Tubs.hh"

class PurgMagDetectorConstruction;
class PurgMagAnalysisManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class PurgMagSteppingAction : public G4UserSteppingAction
{
public:
  PurgMagSteppingAction(const PurgMagDetectorConstruction*);
  ~PurgMagSteppingAction();
  
  void UserSteppingAction(const G4Step*);
  
private:
  const PurgMagDetectorConstruction* Detector; 
  
};

#endif




