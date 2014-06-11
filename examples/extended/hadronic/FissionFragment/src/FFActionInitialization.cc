// -------------------------------------------------------------
//  =============== Begin Documentation Comments ===============
//!
//! \file       FFActionInitialization.cc
//! \author     B. Wendt (brycen.linn.wendt@cern.ch)
//! \date       June 06, 2014
//!
//! \brief      Implementation of the FFActionInitialization class
//!
//  ================ End Documentation Comments ================
//
//  Modified: 
//
// -------------------------------------------------------------

#include "globals.hh"

#include "FFActionInitialization.hh"
#include "FFPrimaryGeneratorAction.hh"
#include "FFRunAction.hh"
//#include "FFEventAction.hh"
//#include "FFSteppingAction.hh"

FFActionInitialization::
FFActionInitialization()
:   G4VUserActionInitialization(),
    masterRunAction(new FFRunAction())
{
    // Nothing here
}

void FFActionInitialization::
Build(void) const
{
    FFRunAction* runAction;
#ifdef G4MULTITHREADED
    runAction = new FFRunAction();
#else
    runAction = masterRunAction;
#endif // G4MULTITHREADED

    SetUserAction(runAction);
    SetUserAction(new FFPrimaryGeneratorAction());
  
    //FFEventAction* eventAction = new FFEventAction();
    //SetUserAction(eventAction);
  
    //SetUserAction(new FFSteppingAction(eventAction));
}

void FFActionInitialization::
BuildForMaster(void) const
{
    SetUserAction(masterRunAction);
}

FFActionInitialization::
~FFActionInitialization()
{
    // Nothing here
}


