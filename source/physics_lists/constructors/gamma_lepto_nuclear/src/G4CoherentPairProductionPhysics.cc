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
// Author:      Alexei Sytov

#include "G4CoherentPairProductionPhysics.hh"
#include "G4Gamma.hh"
#include "G4ProcessManager.hh"
#include "G4CoherentPairProduction.hh"
#include "G4ChannelingFastSimModel.hh"
#include "G4RegionStore.hh"
#include "G4FastSimulationManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CoherentPairProductionPhysics::G4CoherentPairProductionPhysics(const G4String& name):
    G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CoherentPairProductionPhysics::ConstructParticle()
{
    G4Gamma::GammaDefinition();  // Define the gamma particle
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CoherentPairProductionPhysics::ConstructProcess()
{
    //create the gamma process
    G4CoherentPairProduction* gammaProcess = new G4CoherentPairProduction();

    if(verboseLevel > 0) {
        G4cout << "G4CoherentPairProductionPhysics::ConstructProcess" << G4endl;
    }

    G4RegionStore* regionStore = G4RegionStore::GetInstance();
    G4Region* RegionCh = regionStore->GetRegion(fNameRegion);

    //if the region is not found
        if(RegionCh==0)
        {
            G4Exception("GetRegion",// Origin of the exception
                        "001",                   // Unique error code
                        FatalException,          // Terminate the program
                        "Region is not found! The program will terminate.");
        }
        else
        {
            //get channeling model
            G4bool someflag=false;
            G4ChannelingFastSimModel* Channeling =
                static_cast<G4ChannelingFastSimModel*>
                (RegionCh->GetFastSimulationManager()->
                 GetFastSimulationModel(fNameChannelingModel,0,someflag));

            //if channeling model is not found
            if(Channeling==0)
            {
                G4Exception("GetFastSimulationModel",// Origin of the exception
                            "001",                   // Unique error code
                            FatalException,          // Terminate the program
                            "Input channeling model is not found! The program will terminate."
                            );
            }
            else
            {
                gammaProcess->Input(Channeling->GetCrystalData());
            }
        }

    //set functions
    if(fIncoherentScattering){gammaProcess->ActivateIncoherentScattering();}
    gammaProcess->SetLowEnergyLimit(fLowEnergyLimit);
    gammaProcess->SetHighAngleLimit(fHighAngleLimit);
    gammaProcess->SetPPKineticEnergyCut(fPPKineticEnergyCut);
    gammaProcess->SetSamplingPairsNumber(fNMCPairs);
    gammaProcess->SetChargeParticleAngleFactor(fChargeParticleAngleFactor);
    gammaProcess->SetNTrajectorySteps(fNTrajectorySteps);

    //set the name of G4Region in which the model is applicable
    gammaProcess->SetG4RegionName(fNameRegion);

    //get process manager for gamma
    G4ProcessManager* pManager = G4Gamma::Gamma()->GetProcessManager();

    //register the G4CoherentPairProduction process
    pManager->AddDiscreteProcess(gammaProcess);
}
