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

#include "PhysicsList.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmDNAPhysics.hh"
#include "G4EmDNAPhysics_option1.hh"
#include "G4EmDNAPhysics_option2.hh"
#include "G4EmDNAPhysics_option3.hh"
#include "G4EmDNAPhysics_option4.hh"
#include "G4EmDNAPhysics_option5.hh"
#include "G4EmDNAPhysics_option6.hh"
#include "G4EmDNAPhysics_option7.hh"
#include "G4EmDNAPhysics_option8.hh"
#include "G4PhysicsConstructorRegistry.hh"
#include "G4ProcessTable.hh"
#include "G4ProcessManager.hh"
#include "G4DNAChampionElasticModel.hh"
#include "G4DNAScreenedRutherfordElasticModel.hh"
#include "G4DNAElastic.hh"
#include "G4PhysicsListHelper.hh"
#include "G4DNAVibExcitation.hh"
#include "G4DNASancheExcitationModel.hh"
#include "G4DNAElectronSolvation.hh"
#include "G4DNAChemistryManager.hh"
#include "G4LeptonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "PhysicsMessenger.hh"
#include "UserChoosingDNASolvationModel.hh"
#include "InformationKeeper.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList()
    : G4VModularPhysicsList()
{
    fPhysMsg = std::make_unique<PhysicsMessenger>(this);
    SetDefaultCutValue(1.0*nanometer);
    SetVerboseLevel(1);
    RegisterPhysicsList("G4EmDNAPhysics_option2");
    G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100*eV, 1*GeV);
    G4DNAChemistryManager::Instance()->SetChemistryActivation(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
    if(fDNAPhysicsList != nullptr)    
    {
        fDNAPhysicsList->ConstructParticle(); 
    }
    
    // construct following pqrticles to get rid of warning
    G4LeptonConstructor lConstructor;
    lConstructor.ConstructParticle();
    G4BosonConstructor  pBosonConstructor;
    pBosonConstructor.ConstructParticle();
    G4MesonConstructor pMesonConstructor;
    pMesonConstructor.ConstructParticle();
    G4BaryonConstructor pBaryonConstructor;
    pBaryonConstructor.ConstructParticle();
    G4ShortLivedConstructor pShortLivedConstructor;
    pShortLivedConstructor.ConstructParticle();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
    AddTransportation();
    
    if(fDNAPhysicsList != nullptr)    
    {
        fDNAPhysicsList->ConstructProcess(); 
    }

    G4VProcess* process = G4ProcessTable::GetProcessTable()->FindProcess(
                "e-_G4DNAVibExcitation", "e-");

    if (process)
    {
        G4DNAVibExcitation* vibExcitation = (G4DNAVibExcitation*) process;
        G4VEmModel* model = vibExcitation->EmModel();
        G4DNASancheExcitationModel* sancheExcitationMod =
                dynamic_cast<G4DNASancheExcitationModel*>(model);
        if(sancheExcitationMod)
        {
            sancheExcitationMod->ExtendLowEnergyLimit(0.025 * eV);
        }
    }

    // Modify elastic scattering models to avoid killing electrons
    // at low energy
    //
    process = G4ProcessTable::GetProcessTable()->FindProcess("e-_G4DNAElastic", "e-");

    if (process)
    {
        G4DNAElastic* vibExcitation = (G4DNAElastic*) process;
        G4VEmModel* model = vibExcitation->EmModel();

        if(G4DNAChampionElasticModel* championMod =
            dynamic_cast<G4DNAChampionElasticModel*>(model))
        {
        championMod->SetKillBelowThreshold(-1);
        }
        else if(G4DNAScreenedRutherfordElasticModel* screenRutherfordMod =
            dynamic_cast<G4DNAScreenedRutherfordElasticModel*>(model))
        {
        screenRutherfordMod->SetKillBelowThreshold(-1);
        }
    }

    // force to use UserChoosingDNASolvationModel to 
    //invoke G4EmDNAChemistryManager::CreateSovaltedElectron() 
    // when G4DNAChemistryManager is disable
    if (!G4DNAChemistryManager::IsActivated()) {
        process = G4ProcessTable::GetProcessTable()->FindProcess("e-_G4DNAElectronSolvation", "e-");
        if (!process){
            auto pPhysicsListHelper = G4PhysicsListHelper::GetPhysicsListHelper();
            pPhysicsListHelper->RegisterProcess(new G4DNAElectronSolvation("e-_G4DNAElectronSolvation"), 
                                                G4Electron::Definition());
        } 
        G4DNAElectronSolvation* solvation = dynamic_cast<G4DNAElectronSolvation*>(process);
        G4double hLimitE= 0;
        G4VEmModel* therm = solvation->GetModelByIndex(0,1);
        if (therm) {
            
            hLimitE = therm->HighEnergyLimit();
            therm->SetHighEnergyLimit(0*keV);
        }
        auto thermz = UserChoosingDNASolvationModel::UserGetMacroDefinedModel();
        thermz->SetHighEnergyLimit(hLimitE);
        solvation->AddEmModel(-1,thermz);
    }
    
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::RegisterPhysicsList(const G4String& name)
{
    if(name == fPhysDNAName) 
    {
        return; 
    }
    G4bool changedPhys = true;
    if(verboseLevel > 0) 
    {
        G4cout << "===== Register constructor ==== " << name << G4endl; 
    }
    
    if(name == "G4EmDNAPhysics") 
    {
        fDNAPhysicsList.reset(new G4EmDNAPhysics(verboseLevel));
        fPhysDNAName = name;
    } 
    else if(name == "G4EmDNAPhysics_option1") 
    {
        fDNAPhysicsList.reset(new G4EmDNAPhysics_option1(verboseLevel));
        fPhysDNAName = name;
    } 
    else if(name == "G4EmDNAPhysics_option2") 
    {
        fDNAPhysicsList.reset(new G4EmDNAPhysics_option2(verboseLevel));
        fPhysDNAName = name;
    } 
    else if(name == "G4EmDNAPhysics_option3") 
    {
        fDNAPhysicsList.reset(new G4EmDNAPhysics_option3(verboseLevel));
        fPhysDNAName = name;
    } 
    else if(name == "G4EmDNAPhysics_option4") 
    {
        fDNAPhysicsList.reset(new G4EmDNAPhysics_option4(verboseLevel));
        fPhysDNAName = name;
    } 
    else if(name == "G4EmDNAPhysics_option5") 
    {
        fDNAPhysicsList.reset(new G4EmDNAPhysics_option5(verboseLevel));
        fPhysDNAName = name;
    }
    else if(name == "G4EmDNAPhysics_option6") 
    {
        fDNAPhysicsList.reset(new G4EmDNAPhysics_option6(verboseLevel));
        fPhysDNAName = name;
    } 
    else if(name == "G4EmDNAPhysics_option7") 
    {
        fDNAPhysicsList.reset(new G4EmDNAPhysics_option7(verboseLevel));
        fPhysDNAName = name;
    } 
    else if(name == "G4EmDNAPhysics_option8") 
    {
        fDNAPhysicsList.reset(new G4EmDNAPhysics_option8(verboseLevel));
        fPhysDNAName = name;
    } 
    else 
    {
        G4cout << "PhysicsList::RegisterConstructor: <" << name << ">"
               << " fails - name is not defined"
               << G4endl;    
        changedPhys = false;
    }

    if (changedPhys) InformationKeeper::Instance()->RecordPhysDNAName(name);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
