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

#ifndef fluoTestPrimaryGeneratorMessenger_h
#define fluoTestPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class fluoTestPrimaryGeneratorAction;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class fluoTestPrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    fluoTestPrimaryGeneratorMessenger(fluoTestPrimaryGeneratorAction*);
   ~fluoTestPrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
   
  private:
    fluoTestPrimaryGeneratorAction* Action; 
    G4UIcmdWithAString*          RndmCmd;
    G4UIcmdWithAString*          RndmCmmd;
    G4UIcmdWithADoubleAndUnit*     SigmAngleCmd;
    G4UIcmdWithADoubleAndUnit*  SigmaMomentumCmd;
};

#endif








