// -------------------------------------------------------------
//  =============== Begin Documentation Comments ===============
//!
//! \file       FFRunAction.cc
//! \author     B. Wendt (brycen.linn.wendt@cern.ch)
//! \date       June 06, 2014
//!
//! \brief      Implementation of the FFRunAction class
//!
//  ================ End Documentation Comments ================
//
//  Modified: 
//
// -------------------------------------------------------------

#include "globals.hh"

#include "G4RunManager.hh"
#include "G4UserRunAction.hh"

#include "FFRunAction.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
FFRunAction::
FFRunAction()
:   G4UserRunAction()
{
    // Nothing here
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void FFRunAction::
BeginOfRunAction(const G4Run*)
{
    G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void FFRunAction::
EndOfRunAction(const G4Run*)
{
    // TODO Check location of fission fragments here
    // TODO Implement detector tally here
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
FFRunAction::
~FFRunAction()
{
    // Nothing here
}


