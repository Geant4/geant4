// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:		SSAEventAction.hh
//
// Version:		0.b.4
// Date:		16/08/99
// Author:		F Lei
// Organisation:	DERA UK
// Customer:		ESA/ESTEC, NOORDWIJK
// Contract:		12115/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// DESCRIPTION
// -----------
//
// The function of the SSAEventAction is to instantiate the SSAEventMessenger
// class and (at the EndOfEventAction) to draw the event trajectories according 
// to the drawFlag status. 
// 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// PUBLIC MEMBER FUNCTIONS
// -----------------------
//
// SSAEventAction ()
//    Constructor:  Instatiates the SSAEventMessenger class.
//
// ~SSAEventAction ()
//    Destructor:  Deletes the SSAEventMessenger class.
//
// void BeginOfEventAction (const G4Run *aRun)
//    Do nothing    
//
// void EndOfRunAction (const G4Run *aRun)
//    Draw the event trajectories.
//
// void SetDrawFlag(G4String val) 
//    to set the drawFlag using val.
// 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 16 August 1999, F Lei, DERA UK
// Adapted from a verson by Bill Lockman, SLAC, to whom all credits go:
//
// $Id: Tst18EventAction.hh,v 1.1 2000-05-23 06:30:17 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#ifndef Tst18EventAction_h
#define Tst18EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class Tst18EventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst18EventAction : public G4UserEventAction
{
  public:
    Tst18EventAction();
   ~Tst18EventAction();

  public:
    void BeginOfEventAction(const G4Event* anEvent);
    void EndOfEventAction(const G4Event* anEvent);
    
    void SetDrawFlag(G4String val)  {drawFlag = val;};
    
  private:
    G4String drawFlag;                         // control the drawing of event
    Tst18EventActionMessenger*  eventMessenger;
};

#endif

    




