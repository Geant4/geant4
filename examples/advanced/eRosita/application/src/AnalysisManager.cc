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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "AnalysisManager.hh"
#include <sstream>

AnalysisManager *AnalysisManager::instance{nullptr};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AnalysisManager::AnalysisManager()
{
    dataFile1.open("TrackerPhotonEnergy.out"); // open the file
//    dataFile2.open("TotalEnergy.out"); // open the file
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

AnalysisManager::~AnalysisManager()
{
    dataFile1.close(); // close the file
//    dataFile2.close(); // close the file
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

auto AnalysisManager::Instance() -> AnalysisManager*
{
    // A new instance of AnalysisManager is created, if it does not exist:
    if (instance == nullptr) {
        instance = new AnalysisManager();
    }

    // The instance of AnalysisManager is returned:
    return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AnalysisManager::Destroy()
{
    // The AnalysisManager instance is deleted, if it exists:
    if (instance != nullptr) {
        delete instance;
        instance = nullptr;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AnalysisManager::Score(G4double depositedEnergy)
{
    dataFile1 << depositedEnergy << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void AnalysisManager::ScoreTotalEnergy(G4double totalDepositedEnergy)
{
    dataFile2 << totalDepositedEnergy << std::endl;
}
