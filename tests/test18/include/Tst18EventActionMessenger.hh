// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:		SSAEventActionMessenger.hh
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
// The SSAEventActionMessenger is instatiated by the SSARunManager and introduces 
// into the UI event menu additional command to control the drawing of event
// trajectory. User can choose one from
//       1) none: no particle trajectory will be drawn.
//       2) charged: only for charged particles.
//       3) all: for all particles
// The default option is for all particles.  
//        
// The SSAEventActionMessenger modifies the state of the event Drawing flag 
// according to UI menu command issued by the user.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// PUBLIC MEMBER FUNCTIONS
// -----------------------
//
// SSAEventActionMessenger (SSAEventAction*)
//    Constructor:  Defines the commands available to change the DrawFlag
//    status. 
//
// ~SSAEventActionMessenger ()
//    Destructor deletes G4UIdirectory and G4UIcommand objects.
//
// void SetNewValue (G4UIcommand *command, G4String newValues)
//    Identifies the command which has been invoked by the user, extracts the
//    parameters associated with that command (held in newValues, and uses
//    these values with the appropriate member function.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
// CHANGE HISTORY
// --------------
//
// 16 August 1999, F Lei, DERA UK
// Adapted from a verson by Bill Lockman, SLAC, to whom all credits go:
//
// $Id: Tst18EventActionMessenger.hh,v 1.1 2000-05-23 06:30:17 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 

#ifndef Tst18EventActionMessenger_h
#define Tst18EventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst18EventAction;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst18EventActionMessenger: public G4UImessenger
{
  public:
    Tst18EventActionMessenger(Tst18EventAction*);
   ~Tst18EventActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Tst18EventAction*   eventAction;   
    G4UIcmdWithAString* DrawCmd;
};

#endif






