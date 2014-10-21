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
#include "ReactionAction.hh"
#include "G4IT.hh"
#include "G4ITStepManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

ReactionAction::ReactionAction() :
    G4UserTimeStepAction()
{
  /**
   * Give to G4ITStepManager the user defined time steps
   * eg : from 1 picosecond to 10 picosecond, the minimum time
   * step that the TimeStepper can returned is 0.1 picosecond.
   * Those time steps are used for the chemistry of G4DNA
   */

  AddTimeStep(1 * picosecond, 0.1 * picosecond);
  AddTimeStep(10 * picosecond, 1 * picosecond);
  AddTimeStep(100 * picosecond, 3 * picosecond);
  AddTimeStep(1000 * picosecond, 10 * picosecond);
  AddTimeStep(10000 * picosecond, 100 * picosecond);

}

ReactionAction::~ReactionAction()
{
  //dtor
}

ReactionAction::ReactionAction(const ReactionAction& other) :
    G4UserTimeStepAction(other)
{
  //copy ctor
}

ReactionAction& ReactionAction::operator=(const ReactionAction& rhs)
{
  if (this == &rhs) return *this; // handle self assignment
  //assignment operator
  return *this;
}

void ReactionAction::TimeStepAction()
{
//    G4cout << "_________________" << G4endl;
//    G4cout << "Time Step : "
//  << G4BestUnit(G4ITStepManager::Instance()->GetTimeStep(), "Time") << G4endl;
//    G4cout <<  "End of step : "
//  << G4BestUnit(G4ITStepManager::Instance()->GetGlobalTime(), "Time") << G4endl;
}

void ReactionAction::UserReactionAction(const G4Track&,
                                        const G4Track&,
                                        const std::vector<G4Track*>& /*products*/)
{
//    for (int i = 0 ; i < nbProducts ; i ++)
//    {
//        G4cout << "Product[" << i << "] : "
//  << GetIT(products[i])->GetName() << G4endl ;
//    }
}
