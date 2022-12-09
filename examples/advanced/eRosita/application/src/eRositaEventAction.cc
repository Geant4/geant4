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

#include "eRositaEventAction.hh"

#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

eRositaEventAction::eRositaEventAction()
{
}

eRositaEventAction::~eRositaEventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void eRositaEventAction::BeginOfEventAction(const G4Event* event)
{
    auto numberOfEvents = event->GetEventID() + 1;

    constexpr auto INITIAL_FREQUENCY{100};
    constexpr auto FIRST_FREQUENCY_THRESHOLD{1000};
    constexpr auto SECOND_FREQUENCY_THRESHOLD{10000};
    constexpr auto THIRD_FREQUENCY_THRESHOLD{100000};

    auto frequency = INITIAL_FREQUENCY;
    if (numberOfEvents > FIRST_FREQUENCY_THRESHOLD) {
        frequency = FIRST_FREQUENCY_THRESHOLD;
    }
    if (numberOfEvents > SECOND_FREQUENCY_THRESHOLD) {
        frequency = SECOND_FREQUENCY_THRESHOLD;
    }
    if (numberOfEvents > THIRD_FREQUENCY_THRESHOLD) {
        frequency = THIRD_FREQUENCY_THRESHOLD;
    }

    auto remainder = numberOfEvents % frequency;
    if (remainder == 0) {
        G4cout << "---- eRosita number of events: " << numberOfEvents << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void eRositaEventAction::EndOfEventAction(const G4Event* event)
{
    auto eventIdentifier = event->GetEventID();

    // get the number of stored trajectories

    auto *trajectoryContainer = event->GetTrajectoryContainer();
    auto numberOfTrajectories = 0;
    if (trajectoryContainer != nullptr) {
        numberOfTrajectories = trajectoryContainer->entries();
    }

    // periodic printing

    constexpr auto THRESHOLD{100};

    if (eventIdentifier < THRESHOLD || eventIdentifier % THRESHOLD == 0) {
        G4cout << ">>> Event: " << event->GetEventID() << G4endl;
        G4cout << "    " << numberOfTrajectories << " trajectories stored in this event." << G4endl;
    }
}
