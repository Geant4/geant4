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
/// \file scavenger/src/TimeStepAction.cc
/// \brief Implementation of the scavenger::TimeStepAction class

#include "TimeStepAction.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

namespace scavenger
{

TimeStepAction::TimeStepAction() : G4UserTimeStepAction()
{
  /**
   * Give to G4ITTimeStepper the user defined time steps
   * eg : from 1 picosecond to 10 picosecond, the minimum time
   * step that the TimeStepper can returned is 0.1 picosecond.
   * Those time steps are used for the chemistry of G4DNA
   * This method is not recommended for IRT method
   */
  /*
  AddTimeStep(1*picosecond, 0.3*picosecond);
  AddTimeStep(10*picosecond, 1*picosecond);
  AddTimeStep(100*picosecond, 3*picosecond);
  AddTimeStep(1000*picosecond, 10*picosecond);
  AddTimeStep(10000*picosecond, 100*picosecond);
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

TimeStepAction::~TimeStepAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

TimeStepAction::TimeStepAction(const TimeStepAction& other) :
        G4UserTimeStepAction(other)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

TimeStepAction&
TimeStepAction::operator=(const TimeStepAction& rhs)
{
  if (this == &rhs) return *this; // handle self assignment
  //assignment operator
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
void TimeStepAction::UserPreTimeStepAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void TimeStepAction::UserPostTimeStepAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void TimeStepAction::Clear()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void TimeStepAction::UserReactionAction(const G4Track& /*trackA*/,
    const G4Track& /*trackB*/,
    const std::vector<G4Track*>* /*products*/)
{
//  G4cout<<trackA.GetTrackID()<<" + "<<trackB.GetTrackID()<<'\n';
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

}