//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
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


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
FFActionInitialization::
FFActionInitialization()
:   G4VUserActionInitialization(),
    fMasterRunAction(new FFRunAction())
{
    // Nothing here
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void FFActionInitialization::
Build(void) const
{
    FFRunAction* runAction;
#ifdef G4MULTITHREADED
    runAction = new FFRunAction();
#else
    runAction = fMasterRunAction;
#endif // G4MULTITHREADED

    SetUserAction(runAction);
    SetUserAction(new FFPrimaryGeneratorAction());
  
    //FFEventAction* eventAction = new FFEventAction();
    //SetUserAction(eventAction);
  
    //SetUserAction(new FFSteppingAction(eventAction));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void FFActionInitialization::
BuildForMaster(void) const
{
    SetUserAction(fMasterRunAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
FFActionInitialization::
~FFActionInitialization()
{
    // Nothing here
}


