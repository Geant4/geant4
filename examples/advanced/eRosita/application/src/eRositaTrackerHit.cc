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

#include "AnalysisManager.hh"
#include "eRositaTrackerHit.hh"

#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreadLocal G4Allocator<eRositaTrackerHit> *trackerHitAllocator{nullptr};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

eRositaTrackerHit::eRositaTrackerHit()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

eRositaTrackerHit::~eRositaTrackerHit()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

eRositaTrackerHit::eRositaTrackerHit(const eRositaTrackerHit& right) :
    trackIdentifier(right.trackIdentifier),
    chamberNumber(right.chamberNumber),
    depositedEnergy(right.depositedEnergy),
    position(right.position)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

auto eRositaTrackerHit::operator=(const eRositaTrackerHit& right) -> const eRositaTrackerHit&
{
    trackIdentifier = right.trackIdentifier;
    chamberNumber = right.chamberNumber;
    depositedEnergy = right.depositedEnergy;
    position = right.position;

    return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

auto eRositaTrackerHit::operator==(const eRositaTrackerHit& right) const -> G4bool
{
    return this == &right;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void eRositaTrackerHit::Draw()
{
    auto *visualizationManager = G4VVisManager::GetConcreteInstance();
    if (visualizationManager != nullptr) {
        constexpr auto CIRCLE_SCREEN_SIZE = 2.;

        G4Circle circle(position);
        circle.SetScreenSize(CIRCLE_SCREEN_SIZE);
        circle.SetFillStyle(G4Circle::filled);
        G4Colour colour(1., 0., 0.);
        G4VisAttributes attributes(colour);
        circle.SetVisAttributes(attributes);
        visualizationManager->Draw(circle);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void eRositaTrackerHit::Print()
{
    G4cout << "  track identifier: " << trackIdentifier
           << "  deposited energy: " << G4BestUnit(depositedEnergy, "Energy")
           << "  position: " << G4BestUnit(position, "Length") << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void eRositaTrackerHit::PrintToFile() const
{
    // dataFile << trackIdentifier
        // << " " << depositedEnergy
        // << " " << position.x()
        // << " " << position.y()
        // << " " << position.z()
        // << std::endl;
    AnalysisManager::Instance()->Score(depositedEnergy);
}
