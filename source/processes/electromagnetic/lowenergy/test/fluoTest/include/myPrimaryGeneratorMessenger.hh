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

#ifndef myPrimaryGeneratorMessenger_h
#define myPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class myPrimaryGeneratorAction;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class myPrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    myPrimaryGeneratorMessenger(myPrimaryGeneratorAction*);
   ~myPrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
   
  private:
    myPrimaryGeneratorAction* Action; 
    G4UIcmdWithAString*          RndmCmd;
    G4UIcmdWithAString*          RndmCmmd;
    G4UIcmdWithADoubleAndUnit*     SigmAngleCmd;
    G4UIcmdWithADoubleAndUnit*  SigmaMomentumCmd;
};

#endif








