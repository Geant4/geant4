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
// and papers
// M. Batmunkh et al. J Radiat Res Appl Sci 8 (2015) 498-507
// O. Belov et al. Physica Medica 32 (2016) 1510-1520
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// -------------------------------------------------------------------
// November 2016
// -------------------------------------------------------------------
// $Id$
//
/// \file TimeStepAction.hh
/// \brief Implementation of the TimeStepAction class

#include "TimeStepAction.hh"

#include <G4Scheduler.hh>
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
//#include "G4Molecule.hh"

TimeStepAction::TimeStepAction() : G4UserTimeStepAction()
{
  /**
   * Give to G4ITTimeStepper the user defined time steps
   * eg : from 1 picosecond to 10 picosecond, the minimum time
   * step that the TimeStepper can returned is 0.1 picosecond.
   * Those time steps are used for the chemistry of G4DNA
   */

  AddTimeStep(1*picosecond, 0.1*picosecond);
  AddTimeStep(10*picosecond, 1*picosecond);
  AddTimeStep(100*picosecond, 3*picosecond);
  AddTimeStep(1000*picosecond, 10*picosecond);
  AddTimeStep(10000*picosecond, 100*picosecond);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepAction::~TimeStepAction()
{
  //dtor
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TimeStepAction::TimeStepAction(const TimeStepAction& other) :
        G4UserTimeStepAction(other)
{
  //copy ctor
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

void TimeStepAction::UserPostTimeStepAction()
{
  //    G4cout << "_________________" << G4endl;
  /*
  G4cout << "Time Step : "
      << G4BestUnit(G4ITScheduler::Instance()->GetTimeStep(),
          "Time")
          << G4endl;

  G4cout <<  "End of step: "
      << G4BestUnit(G4ITScheduler::Instance()->GetGlobalTime(),
                            "Time")
      << G4endl;
   */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TimeStepAction::UserReactionAction(const G4Track&,
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
