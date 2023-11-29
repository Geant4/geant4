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
/// \file B5/include/EventAction.hh
/// \brief Definition of the B5::EventAction class

#ifndef B5EventAction_h
#define B5EventAction_h 1

#include "Constants.hh"

#include "G4UserEventAction.hh"
#include "globals.hh"

#include <vector>
#include <array>

// named constants
const G4int kEm = 0;
const G4int kHad = 1;
const G4int kH1 = 0;
const G4int kH2 = 1;
const G4int kDim = 2;

namespace B5
{

/// Event action

class EventAction : public G4UserEventAction
{
public:
    EventAction();
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction(const G4Event*) override;

    std::vector<G4double>& GetEmCalEdep() { return fCalEdep[kEm]; }
    std::vector<G4double>& GetHadCalEdep() { return fCalEdep[kHad]; }

private:
    // hit collections Ids
    std::array<G4int, kDim> fHodHCID = { -1, -1 };
    std::array<G4int, kDim> fDriftHCID = { -1, -1 };
    std::array<G4int, kDim> fCalHCID = { -1, -1 };
    // histograms Ids
    std::array<std::array<G4int, kDim>, kDim> fDriftHistoID
      {{ {{ -1, -1 }}, {{ -1, -1 }} }};
        // std::array<T, N> is an aggregate that contains a C array.
        // To initialize it, we need outer braces for the class itself
        // and inner braces for the C array
    // energy deposit in calorimeters cells
    std::array<std::vector<G4double>, kDim> fCalEdep
      {{ std::vector<G4double>(kNofEmCells, 0.), std::vector<G4double>(kNofHadCells, 0.) }};
};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
