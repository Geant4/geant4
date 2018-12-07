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
//
/// \file src/RE04ActionInitialization.cc
/// \brief Implementation of the RE04ActionInitialization class
//

#include "RE04ActionInitialization.hh"
#include "RE04PrimaryGeneratorAction.hh"
#include "RE04EventAction.hh"
//#include "RE04TrackingAction.hh"
//#include "RE04SteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE04ActionInitialization::RE04ActionInitialization()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE04ActionInitialization::~RE04ActionInitialization()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE04ActionInitialization::Build() const
{
  //
  SetUserAction(new RE04PrimaryGeneratorAction);
  //
  SetUserAction(new RE04EventAction);
  //
  //SetUserAction(new RE04TrackingAction);
  //SetUserAction(new RE04SteppingAction);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE04ActionInitialization::BuildForMaster() const
{
  //
  //G4UserRunAction* run_action = new RE04RunAction;
  //SetUserAction(run_action);
}

