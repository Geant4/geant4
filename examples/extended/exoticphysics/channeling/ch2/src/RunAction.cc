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
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "G4AnalysisManager.hh"
#include "G4RegionStore.hh"
#include "G4FastSimulationManager.hh"
#include "G4ChannelingFastSimModel.hh"
#include "G4Threading.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
    : G4UserRunAction()
{
    //using analysis manager for output
    auto analysisManager = G4AnalysisManager::Instance();

#ifdef G4MULTITHREADED
    analysisManager->SetNtupleMerging(true);
#else
    analysisManager->SetNtupleMerging(false);
#endif

    //Creating the ntuple to score the particles and the emitted radiation

    //ALL "detector_primaries" ENTER THE CRYSTAL; "detector_photons", "detector_secondaries"
    //are their daughters.

    //CAUTION: if a primary did not cross the crystal,
    //this primary and its daughters are written ONLY in missed_crystal

    G4String nTupleName[5] =
        {"crystal",
                              "detector_primaries",
                              "detector_photons",
                              "detector_secondaries",
                              "missed_crystal"};
    for(G4int i=0; i<5; i++)
    {
        analysisManager->CreateNtuple(nTupleName[i],nTupleName[i]);
        analysisManager->CreateNtupleIColumn("eventID");
        analysisManager->CreateNtupleSColumn("volume");
        analysisManager->CreateNtupleDColumn("x");
        analysisManager->CreateNtupleDColumn("y");
        analysisManager->CreateNtupleDColumn("angle_x");
        analysisManager->CreateNtupleDColumn("angle_y");
        analysisManager->CreateNtupleDColumn("Ekin");
        analysisManager->CreateNtupleSColumn("particle");
        analysisManager->CreateNtupleIColumn("particleID");
        analysisManager->CreateNtupleIColumn("parentID");

        if(i==1)
        {
            analysisManager->CreateNtupleDColumn("incoming_angle_x");
            analysisManager->CreateNtupleDColumn("deflection_angle_x");
            analysisManager->CreateNtupleDColumn("incoming_angle_y");
            analysisManager->CreateNtupleDColumn("deflection_angle_y");
        }

        analysisManager->FinishNtuple();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
    //opening output file
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    G4String fileName = "results.root";
    analysisManager->OpenFile(fileName);

    //delete the spectrum temporary files if they exist
    std::string filename = "Spectrum_"+std::to_string(G4Threading::G4GetThreadId())+".dat";
    std::remove(filename.c_str());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();

    //getting internal data of G4ChannelingFastSimModel
    G4RegionStore* regionStore = G4RegionStore::GetInstance();
    G4Region* regionCh = regionStore->GetRegion("Crystal");
    G4bool someflag=false;
    G4ChannelingFastSimModel* channeling =
        static_cast<G4ChannelingFastSimModel*>
        (regionCh->GetFastSimulationManager()->GetFastSimulationModel("ChannelingModel",
                                                                      0,someflag));

    if (!IsMaster() && channeling->GetIfRadiationModelActive())
    {
        std::vector<G4double> photonEnergyInSpectrum =
            channeling->GetRadiationModel()->GetPhotonEnergyInSpectrum();
        std::vector<G4double> spectrum =
            channeling->GetRadiationModel()->GetTotalSpectrum();

        G4int threadID = G4Threading::G4GetThreadId();

        std::ofstream file1;
        file1.open("Spectrum_"+std::to_string(threadID)+".dat");

        file1 << std::setprecision(16);
        for(std::size_t i = 0; i<spectrum.size(); i++)
            {file1 << photonEnergyInSpectrum[i] << " " << spectrum[i] << G4endl;}

        file1.close();
    }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
