// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: F02SteppingAction.hh,v 1.1 2001-03-27 16:26:21 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef F02SteppingAction_h
#define F02SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "globals.hh"

class F02DetectorConstruction;
class F02RunAction;
class F02EventAction;
class F02SteppingMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class F02SteppingAction : public G4UserSteppingAction
{
  public:
    F02SteppingAction(F02DetectorConstruction*, F02EventAction*,
                      F02RunAction* );
   ~F02SteppingAction();

    void UserSteppingAction(const G4Step*);

  private:
    F02DetectorConstruction* detector;
    F02EventAction*          eventaction;
    F02RunAction*            runaction;
    F02SteppingMessenger*    steppingMessenger;

    G4int IDnow,IDold;
    G4int evnoold ;

};

#endif
