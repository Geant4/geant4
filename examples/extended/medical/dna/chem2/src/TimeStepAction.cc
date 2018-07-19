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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// J. Comput. Phys. 274 (2014) 841-882
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// $Id$
//
/// \file TimeStepAction.cc
/// \brief Implementation of the TimeStepAction class

#include "TimeStepAction.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
//#include "G4ITScheduler.hh"
//#include "G4Molecule.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepAction::TimeStepAction() :
    G4UserTimeStepAction()
{
  /**
   * Inform G4ITTimeStepper of the selected minimum time steps
   * eg : from 1 picosecond to 10 picosecond, the minimum time
   * step that the TimeStepper can returned is 0.1 picosecond.
   *
   * Case 1) If the rection model calculates a minimum reaction time
   * bigger than the user defined time step, the reaction model wins
   *
   * Case 2) If an interaction process with the continuous medium
   * calculates a time step less than the selected minimum time step,
   * the interaction process wins
   */

  AddTimeStep(1 * picosecond, 0.1 * picosecond);
  AddTimeStep(10 * picosecond, 1 * picosecond);
  AddTimeStep(100 * picosecond, 10 * picosecond);
  AddTimeStep(1000 * picosecond, 100 * picosecond);
  AddTimeStep(10000 * picosecond, 1000 * picosecond);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepAction::~TimeStepAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepAction::TimeStepAction(const TimeStepAction& other) :
    G4UserTimeStepAction(other)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TimeStepAction&
TimeStepAction::operator=(const TimeStepAction& rhs)
{
  if (this == &rhs) return *this; // handle self assignment
  //assignment operator
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TimeStepAction::StartProcessing()
{
// You want to know why the simulation stopped ?
// G4ITScheduler::Instance()->WhyDoYouStop();
// At the end of the simulation, information will be printed
// It is better to place this command before the simulation starts
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TimeStepAction::UserPreTimeStepAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TimeStepAction::UserPostTimeStepAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// Here you can retrieve information related to reactions
void
TimeStepAction::UserReactionAction(const G4Track&,
                                   const G4Track&,
                                   const std::vector<G4Track*>* /*products*/)
{
  /*
   for (int i = 0 ; i < nbProducts ; i ++)
   {
   G4cout << "Product[" << i << "] : "
   << GetMolecule(products[i])->GetName()
   << G4endl ;
   }
   */
}

void TimeStepAction::EndProcessing()
{
}
