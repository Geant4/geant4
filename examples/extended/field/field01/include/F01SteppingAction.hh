// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: F01SteppingAction.hh,v 1.1 2001-03-27 16:21:29 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef F01SteppingAction_h
#define F01SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "globals.hh"

class F01DetectorConstruction;
class F01RunAction;
class F01EventAction;
class F01SteppingMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class F01SteppingAction : public G4UserSteppingAction
{
  public:
    F01SteppingAction(F01DetectorConstruction*, F01EventAction*,
                      F01RunAction* );
   ~F01SteppingAction();

    void UserSteppingAction(const G4Step*);

  private:
    F01DetectorConstruction* detector;
    F01EventAction*          eventaction;
    F01RunAction*            runaction;
    F01SteppingMessenger*    steppingMessenger;

    G4int IDnow,IDold;
    G4int evnoold ;

};

#endif
