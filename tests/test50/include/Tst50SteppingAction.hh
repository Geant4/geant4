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
// $Id: Tst50SteppingAction.hh,v 1.10 2003-02-06 14:42:37 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Tst50SteppingAction_h
#define Tst50SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class G4Step;
class Tst50AnalysisManager;
class Tst50EventAction;
class Tst50RunAction;
class Tst50PrimaryGeneratorAction;
class Tst50DetectorConstruction;
class Tst50SteppingAction : public G4UserSteppingAction
{
  public:
    Tst50SteppingAction(Tst50EventAction*,Tst50PrimaryGeneratorAction*,Tst50RunAction*, Tst50DetectorConstruction*,G4String,G4bool,G4bool,G4bool, G4bool);
   ~Tst50SteppingAction();

    void UserSteppingAction(const G4Step* Step);

private:
 G4int IDnow;
 G4int IDold;
 Tst50EventAction*          eventaction;
  Tst50PrimaryGeneratorAction* p_Primary;
  Tst50RunAction* runaction; 
  Tst50DetectorConstruction* detector;     
  G4double initial_energy;
  G4double range;
  G4double  KinE_stepBeginning;
  G4String filename;
  G4bool StoppingPower;
  G4bool Range;
  G4bool RadiationY;
  G4bool Foil;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
