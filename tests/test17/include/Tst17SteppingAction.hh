// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst17SteppingAction.hh,v 1.2 1999-12-15 14:54:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Tst17SteppingAction_h
#define Tst17SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "globals.hh"

class Tst17DetectorConstruction;
class Tst17RunAction;
class Tst17EventAction;
class Tst17SteppingMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst17SteppingAction : public G4UserSteppingAction
{
  public:
    Tst17SteppingAction(Tst17DetectorConstruction* DET,Tst17EventAction* EA,
                        Tst17RunAction* RA);
   ~Tst17SteppingAction();

    void UserSteppingAction(const G4Step*);

  private:
    Tst17DetectorConstruction* detector;
    Tst17EventAction* eventaction;
    Tst17RunAction* runaction;
    Tst17SteppingMessenger* steppingMessenger;

    G4int IDnow,IDold;
    G4int evnoold ;

};

#endif
